import os, re, logging, json, gzip
from glob import glob
from datetime import datetime

#=== timing report ================
class TimeReport:
    def __init__(self, name):
        self.mName = name
        self.mStartTime = datetime.now()
        logging.info("Process %s started at %s"
            % (self.mName, self.mStartTime))

    def rate(self, total, dt):
        return total / (dt.seconds + 1.)

    def portion(self, total, at_msg = None):
        dt = datetime.now() - self.mStartTime
        note = self.mName
        if at_msg:
            note += " at:" + str(at_msg)
        logging.info("%s Records: %d Time: %s; Rate: %.2f"
            % (note, total, str(dt), self.rate(total, dt)))

    def done(self, total):
        end_time = datetime.now()
        dt = end_time - self.mStartTime
        logging.info("Process %s done at %s; Records: %d; Rate: %.2f"
            % (self.mName, end_time, total, self.rate(total, dt)))

#=== chromosome detection ================
def detectFileChrom(filename, chrom_loc):
    chrom_patt = re.compile(r"(\b|\W)%s(\w+)(\b|\W)" % chrom_loc, re.I)
    chrom_list = {str(idx) for idx in range(1, 23)} | {"M", "X", "Y"}
    qq = chrom_patt.search(os.path.basename(filename))
    assert qq is not None and qq.group(2).upper() in chrom_list, (
        "Failed to detect chrom in filename: " + filename)
    return qq.group(2).upper()

#=== file list extension ================
def extendFileList(files):
    if isinstance(files, str):
        files = [files]
    result = []
    for fname in files:
        if '*' in fname:
            result += list(glob(fname))
        else:
            result.append(fname)
    return sorted(result)

#==== Joined readers =====================
class JoinedReader:
    def __init__(self, name, readers, max_count = -1):
        self.mName = name
        self.mReaders = readers
        self.mBuffers = None
        self.mDone = False
        self.mMaxCount = max_count
        self.mTotal, self.mCount = 0, 0
        self.mTimeReport = TimeReport(self.mName)
        self.mFixedChrom = None

    def iterReaders(self):
        return iter(self.mReaders)

    def close(self):
        for reader in self.mReaders:
            reader.close()
        self.mTimeReport.done(self.mTotal)

    def isDone(self):
        return self.mDone

    def nextOne(self):
        if self.mDone:
            return None
        if self.mBuffers is None:
            self.mBuffers = [reader.getNext() for reader in self.mReaders]
            for buf in self.mBuffers:
                if buf is not None:
                    if self.mFixedChrom is None:
                        self.mFixedChrom = buf[0][0]
                    else:
                        assert self.mFixedChrom == buf[0][0]
        min_key = None
        for idx, reader in enumerate(self.mReaders):
            buf = self.mBuffers[idx]
            if buf is None:
                self.mBuffers[idx] = reader.getNext()
                buf = self.mBuffers[idx]
            if buf is not None:
                if min_key is None or min_key > buf[0]:
                    min_key = buf[0]
        if min_key is None:
            self.mDone = True
            return None
        assert min_key[0] == self.mFixedChrom
        res_seq = []
        for idx in range(len(self.mReaders)):
            buf = self.mBuffers[idx]
            if buf is not None and min_key == buf[0]:
                res_seq += buf[1]
                self.mBuffers[idx] = None
        self.mTotal += len(res_seq)
        self.mCount += 1
        if self.mMaxCount > 0 and self.mCount >= self.mMaxCount:
            logging.info("%s: stopped by count limit %d"
                % (self.mName, self.mMaxCount))
            self.mDone = True
        if self.mCount % 100000 == 0:
            self.mTimeReport.portion(self.mTotal,  "%s:%s:" % min_key)
        return [min_key, res_seq]

#=====================================
def dumpReader(reader, indent_mode = False):
    for key, record in reader.read():
        print(json.dumps({"key": list(key)},
            ensure_ascii = False, sort_keys = True))
        if indent_mode:
            rep = json.dumps(record, indent = 4,
                ensure_ascii = False, sort_keys = True)
        else:
            rep = json.dumps(record,
                ensure_ascii = False, sort_keys = True)
        print(rep)

#========================================
def writeDirect(reader, file_pattern, out_dir):
    cur_chrom, cur_outp = None, None
    for key, rec in reader.read():
        if key[0] != cur_chrom:
            if cur_outp is not None:
                cur_outp.close()
            cur_chrom = key[0]
            fname = out_dir + "/" + (file_pattern % cur_chrom)
            logging.info("Writing %s..." % fname)
            cur_outp = gzip.open(fname, "wt", encoding = "utf-8")
        print(json.dumps([key, rec], sort_keys = True,
            ensure_ascii = False), file = cur_outp)
    if cur_outp is not None:
        cur_outp.close()
    logging.info("File preparation done")

#========================================
class DirectReader:
    def __init__(self, file_list):
        self.mFNames = sorted(extendFileList(file_list), reverse = True)

    def read(self):
        for fname in self.mFNames:
            time_rep = TimeReport(fname)
            count = 0
            logging.info("Loading: " + fname)
            with gzip.open(fname, "rt", encoding = "utf-8") as inp:
                for line in inp:
                    key, rec = json.loads(line)
                    count += 1
                    if count % 100000 == 0:
                        time_rep.portion(count)
                    yield [tuple(key), rec]
            time_rep.done(count)
