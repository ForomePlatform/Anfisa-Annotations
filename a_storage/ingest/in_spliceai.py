import sys, logging, gzip
import cyvcf2 as pyvcf
from fastnumbers import fast_real

from .in_util import (JoinedReader, extendFileList, dumpReader,
    writeDirect, DirectReader, detectFileChrom, TimeReport)
#========================================
# Schema for AStorage
#========================================
SCHEMA_SPLICE_AI = {
    "name": "SpliceAI",
    "key": "hg38",
    "no-stat-details": True,
    "io": {
        "block-type": "segment",
        "pos-frame": 256
    },
    "filter-list": {"ref": "REF", "alt": "ALT"},
    "top": {
        "tp": "list",
        "item": {
            "tp": "dict",
            "label": "spliceai-rec",
            "items": [
                {"name": "ALT", "tp": "str", "opt": "gene"},
                {"name": "REF", "tp": "str", "opt": "gene"},
                {"name": "ID", "tp": "str", "opt": "gene"},
                {"name": "SYMBOL", "tp": "str", "opt": "repeat"},
                {"name": "DP_AG",  "tp": "num", "format": "%d"},
                {"name": "DP_AL",  "tp": "num", "format": "%d"},
                {"name": "DP_DG",  "tp": "num", "format": "%d"},
                {"name": "DP_DL",  "tp": "num", "format": "%d"},
                {"name": "MAX_DS", "tp": "num"},
                {"name": "DS_AG",  "tp": "num"},
                {"name": "DS_DG",  "tp": "num"},
                {"name": "DS_AL",  "tp": "num"},
                {"name": "DS_DL",  "tp": "num"}
            ]
        }
    }
}

#========================================
# Ingest logic
#========================================
# VCF_INFO_NAMES = (
#    "ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split('|'))

VCF_INFO_TYPES = [
    ('ALT',     str),
    ('SYMBOL',  str),
    ('DS_AG',   fast_real),
    ('DS_AL',   fast_real),
    ('DS_DG',   fast_real),
    ('DS_DL',   fast_real),
    ('DP_AG',   fast_real),
    ('DP_AL',   fast_real),
    ('DP_DG',   fast_real),
    ('DP_DL',   fast_real)]

DB_FIELDS = [
    "ALT", "REF", "ID", "SYMBOL",
    "DP_AG", "DP_AL", "DP_DG", "DP_DL",
    "MAX_DS", "DS_AG", "DS_DG", "DS_AL", "DS_DL"]

LEN_REC = len(DB_FIELDS)
REC_SHEET = []

for idx, fld_name in enumerate(DB_FIELDS):
    rec_info = (None, None)
    for vcf_idx, info in enumerate(VCF_INFO_TYPES):
        nm, tp = info
        if nm == fld_name:
            rec_info = (vcf_idx, tp)
            break
    if rec_info[0] is None:
        assert fld_name in ("REF", "ID", "MAX_DS")
        assert idx in (1, 2, 8)
    else:
        REC_SHEET.append((idx, rec_info[0], rec_info[1]))

#========================================
def new_single_record(in_record):
    global REC_SHEET, LEN_REC

    chrom = in_record.CHROM
    pos = int(in_record.POS)
    record = [None] * LEN_REC
    record[1] = in_record.REF
    record[2] = in_record.ID

    info_fields = in_record.INFO['SpliceAI'].split('|')
    for idx, vcf_idx, tp in REC_SHEET:
        record[idx] = tp(info_fields[vcf_idx])
    record[8] = max(record[9:])
    if record[8] == 0:
        record = record[:8]
    else:
        while record[-1] == 0:
            del record[-1]

    return [("chr" + chrom, pos), record]

#===========================================================
class InputDataReader:
    def __init__(self, chrom, fname):
        self.mFName = fname
        self.mReader = pyvcf.VCF(fname)
        self.mCurRecord = False
        self.mChrom = chrom
        if not self.mChrom.startswith("chr"):
            self.mChrom = "chr" + self.mChrom

    def getFName(self):
        return self.mFName

    def isOver(self):
        return self.mCurRecord is None

    def close(self):
        self.mReader.close()

    def getNext(self):
        if self.mCurRecord is None:
            return None
        if self.mCurRecord is False:
            key, seq = None, None
        else:
            key, seq = self.mCurRecord[0], [self.mCurRecord[1]]
        while True:
            try:
                in_record = next(self.mReader)
            except Exception:
                self.mCurRecord = None
                return [key, seq]
            self.mCurRecord = new_single_record(in_record)
            if key is None:
                key, seq = self.mCurRecord[0], [self.mCurRecord[1]]
                continue
            if self.mCurRecord[0] == key:
                seq.append(self.mCurRecord[1])
                continue
            assert key is None or (key[0] != self.mCurRecord[0][0]
                or key[1] < self.mCurRecord[0][1])
            assert key[0] == self.mChrom, (
                "Chrom conflict: %s/%s" % (key[0], self.mChrom))
            return [key, seq]

#===========================================================
class ReaderSpliceAI:
    def __init__(self, indel_file_list, snv_file_list,
            chrom_loc = ".chr"):
        self.mIndelFiles = dict()
        for fname in extendFileList(indel_file_list):
            chrom = detectFileChrom(fname, chrom_loc)
            assert chrom not in self.mIndelFiles
            self.mIndelFiles[chrom] = fname
        self.mSnvFiles = dict()
        for fname in extendFileList(snv_file_list):
            chrom = detectFileChrom(fname, chrom_loc)
            assert chrom not in self.mSnvFiles
            self.mSnvFiles[chrom] = fname
        self.mChromSeq = sorted(
            set(self.mIndelFiles.keys()) | set(self.mSnvFiles.keys()),
            reverse = True)
        assert len(self.mChromSeq) > 0

    def read(self):
        for chrom in self.mChromSeq:
            readers = []
            fname = self.mIndelFiles.get(chrom)
            if fname:
                readers.append(InputDataReader(chrom, fname))
            fname = self.mSnvFiles.get(chrom)
            if fname:
                readers.append(InputDataReader(chrom, fname))
            logging.info(("Evaluation of %s chromosome:" % chrom)
                + " ".join([reader.getFName() for reader in readers]))
            join_reader = JoinedReader("chr%s" % chrom, readers)
            while not join_reader.isDone():
                ret = join_reader.nextOne()
                if ret is not None:
                    yield ret
            join_reader.close()

#========================================
def reader_SpliceAI(properties, schema_h = None):
    if schema_h is not None:
        schema_h.getCodecByLabel(
            "spliceai-rec").setSerialization(DB_FIELDS)

    if properties is None:
        return None

    if "direct_file_list" in properties:
        return DirectReader(properties["direct_file_list"])

    return ReaderSpliceAI(
        properties["indel_file_list"],
        properties["snv_file_list"])

#========================================
def splitPreparation(in_long_file, out_dir):
    cur_chrom = None
    cur_fname = None
    output = None
    out_count = None
    used_chroms = set()
    headers = []

    if "indel" in in_long_file:
        mode = "indel"
    elif "snv" in in_long_file:
        mode = "snv"
    else:
        assert False, "Mode (indel/snv) not found"
    file_pattern = out_dir + "/" + mode + ".spliceai.chr%s.vcf.gz"
    time_rep = TimeReport(mode)
    last_line_no = 0

    with gzip.open(in_long_file, "rt") as inp:
        for line_no, line in enumerate(inp):
            last_line_no = line_no
            if line.startswith('#'):
                headers.append(line)
                continue
            chrom = line[:line.index('\t')]
            if chrom != cur_chrom:
                if output is not None:
                    output.close()
                    logging.info("Done %s: %d records"
                        % (cur_fname, out_count))
                    output = None
                    out_count = None
                cur_chrom = chrom
                if '_' not in cur_chrom:
                    assert cur_chrom not in used_chroms
                    used_chroms.add(cur_chrom)
                    cur_fname = file_pattern % cur_chrom
                    logging.info("Open %s" % cur_fname)
                    output = gzip.open(cur_fname, "wt")
                    for hline in headers:
                        output.write(hline)
                    out_count = 0
            if output is not None:
                output.write(line)
                out_count += 1
            if (line_no % 100000) == 0:
                pos = line.split('\t')[1]
                time_rep.portion(line_no, "%s:%s:" % (cur_chrom, pos))
        if output is not None:
            output.close()
            logging.info("Done %s: %d records"
                % (cur_fname, out_count))
    time_rep.done(last_line_no + 1)


#========================================
if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)
    if sys.argv[1] == "SPLIT":
        splitPreparation(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "DIR":
        logging.info("Using direct preparation in " + sys.argv[4])
        reader = reader_SpliceAI({
            "indel_file_list": sys.argv[2],
            "snv_file_list": sys.argv[3]})
        writeDirect(reader, "spliceai_dir_%s.js.gz", sys.argv[4])
    else:
        reader = reader_SpliceAI({
            "indel_file_list": sys.argv[1],
            "snv_file_list": sys.argv[2]})
        dumpReader(reader)
