import json, random, logging

from forome_tools.log_err import logException
from .a_blocker import ABlockerIO_Single
#========================================
class AFastaSchema:
    def __init__(self, storage, name, dbname, schema_descr = None):
        self.mStorage = storage
        self.mName = name
        self.mWriteMode = schema_descr is not None

        self.mSchemaDescr = self.mStorage.preLoadSchemaData(
            dbname, self.mName, schema_descr)

        self.mTotals = self.mSchemaDescr.get("total", [])

        self.mDbConnector = self.getStorage().openConnection(
            dbname, self.isWriteMode())

        self.mTypes = dict()
        for idx, type_name in enumerate(self.mSchemaDescr["types"]):
            blocker = ABlockerIO_Single(self, {
                "block-type": "single",
                "fasta-col-options": {"-compress": "bz2"}},
                key_codec_type = type_name, col_type = "fasta",
                col_name = type_name, conv_bytes = True)
            if len(self.mTotals) <= idx:
                self.mTotals.append(0)
            tp_h = _FastaTypeHandler(type_name, blocker, idx)
            self.mTypes[type_name] = tp_h
            if (self.mWriteMode and not self.mStorage.isDummyMode()):
                tp_h.initSamples(self.mStorage.getSamplesCount())

        self.mBlockSize = self.mSchemaDescr["block-size"]
        if self.mWriteMode:
            self.mTypesSet = set()
        logging.info("Fasta starts with types: %s"
            % ", ".join(sorted(self.mTypes.keys())))

    def flush(self):
        for tp_h in self.mTypes.values():
            tp_h.getBlockIO().flush()
        self.keepSchema()

    def getDbConnector(self):
        return self.mDbConnector

    def keepSchema(self):
        if not self.mWriteMode or self.mStorage.isDummyMode():
            return
        self.mSchemaDescr["total"] = self.mTotals
        self.mStorage.saveSchemaData(self)
        logging.info("Schema %s kept, total = %s"
            % (self.mName, ",".join(map(str, self.mTotals))))

    def close(self):
        self.flush()
        if self.mWriteMode:
            assert len(self.mTypesSet) == len(self.mTypes), (
                "No data loaded for "
                + " ".join(sorted(set(self.mTypes.keys()) - self.mTypesSet)))
        self.keepSchema()
        self.keepSamples()
        if self.mWriteMode and not self.mStorage.isDummyMode():
            self.checkSamples()

    def getStorage(self):
        return self.mStorage

    def getName(self):
        return self.mName

    def getDbName(self):
        return self.mDbConnector.getName()

    def isWriteMode(self):
        return self.mWriteMode

    def getProperty(self, name):
        return self.mSchemaDescr.get(name)

    def getTotal(self):
        return self.mTotals

    def isOptionRequired(self, opt):
        return False

    def getFilteringProperties(self):
        return ["type"]

    def useLastPos(self):
        return True

    def getDBKeyType(self):
        return "fasta"

    def getSchemaDescr(self):
        return self.mSchemaDescr

    def iterTypeNames(self):
        return iter(self.mTypes.keys())

    def loadReader(self, reader):
        hg_type = reader.getName()
        tp_h = self.mTypes[hg_type]
        self.mTypesSet.add(hg_type)
        cur_chrom, prev_diap = None, None
        for chrom, diap, letters in reader.readAll(self.mBlockSize):
            if chrom != cur_chrom:
                assert diap[0] == 1
                cur_chrom = chrom
                prev_diap = None
            elif prev_diap is not None:
                assert prev_diap[1] == diap[0]
            assert (diap[0] - 1) % self.mBlockSize == 0
            tp_h.getBlockIO().putRecord(
                ("chr" + chrom, diap[0] - 1), letters)
            self.mTotals[tp_h.getIdx()] += 1
            tp_h.addSample(chrom, diap, letters)
            prev_diap = diap

    def getRecord(self, key, filtering, last_pos = None):
        if (filtering is None or "type" not in filtering
                or filtering["type"] not in self.mTypes):
            raise Exception("Request requires 'type' argument = %s"
                % sorted(self.mTypes.keys()))

        tp_h = self.mTypes[filtering["type"]]
        chrom, pos = key
        base_pos = (pos - 1) - ((pos - 1) % self.mBlockSize)
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        letters = tp_h.getBlockIO().getRecord((chrom, base_pos))
        loc_pos = pos - 1 - base_pos
        if letters is None or len(letters) < loc_pos:
            return None
        if last_pos is None:
            return letters[loc_pos]
        loc_last = last_pos - base_pos - 1
        if loc_last + 1 < len(letters):
            return letters[loc_pos:loc_last + 1]
        ret = [letters[loc_pos:]]
        while len(letters) == self.mBlockSize:
            base_pos += self.mBlockSize
            loc_last -= self.mBlockSize
            letters = tp_h.getBlockIO().getRecord((chrom, base_pos))
            if letters is None:
                break
            if loc_last + 1 < len(letters):
                ret.append(letters[:loc_last + 1])
                break
            ret.append(letters)
        return ''.join(ret)

    def keepSamples(self):
        if not self.mWriteMode or self.mStorage.isDummyMode():
            return
        with open(self.mStorage.getSchemaFilePath(self, "1.samples"),
                "w", encoding = "utf-8") as output:
            cnt = 0
            for tp_h in self.mTypes.values():
                for chrom, diap, letters in tp_h.iterSamples():
                    print(json.dumps({
                        "no": cnt + 1,
                        "tp": tp_h.getName(),
                        "chrom": chrom,
                        "diap": diap,
                        "letters": letters}), file = output)
                    cnt += 1

    def checkSamples(self):
        smp_input = open(self.mStorage.getSchemaFilePath(self, "1.samples"),
            "r", encoding = "utf-8")

        cnt_ok, cnt_bad, cnt_fail = 0, 0, 0
        with open(self.mStorage.getSchemaFilePath(self, "2.samples"),
                "w", encoding = "utf-8") as output:
            while True:
                try:
                    line_rep = smp_input.readline()
                    if not line_rep:
                        break
                    rep_obj = json.loads(line_rep)
                    letters1 = self.getRecord(
                        (rep_obj["chrom"], rep_obj["diap"][0]),
                        {"type": rep_obj["tp"]}, rep_obj["diap"][1])
                    if rep_obj["letters"] != letters1:
                        cnt_bad += 1
                        rep_obj["ok"] = False
                        rep_obj["r-letters"] = letters1
                    else:
                        cnt_ok += 1
                        rep_obj["ok"] = True
                    print(json.dumps(rep_obj, ensure_ascii = False),
                        file = output)
                except Exception:
                    msg_txt = logException("Check samples")
                    cnt_fail += 1
                    print(json.dumps({"exception": msg_txt},
                        ensure_ascii = False), file = output)
            print(json.dumps({"tp": "result",
                "ok": cnt_ok, "bad": cnt_bad, "fail": cnt_fail}),
                file = output)
        smp_input.close()
        if cnt_bad + cnt_fail == 0:
            logging.info("Samples check(%d) for %s: OK" %
                (cnt_ok, self.mName))
        else:
            logging.error("BAD! Samples check for %s: %d of %d (+ %d failures)"
                % (self.mName, cnt_bad, cnt_bad + cnt_ok, cnt_fail))

#========================================
class _FastaTypeHandler:
    def __init__(self, name, block_io, idx):
        self.mName = name
        self.mBlockIO = block_io
        self.mIdx = idx
        self.mSmpCount = None
        self.mTotal = None
        self.mSamples = None
        self.mRH = None
        self.mPrevIdx = None

    def getName(self):
        return self.mName

    def getIdx(self):
        return self.mIdx

    def getBlockIO(self):
        return self.mBlockIO

    def initSamples(self, smp_count):
        self.mRH = random.Random(179)
        self.mSmpCount = smp_count
        self.mSamples = []
        self.mTotal = 0

    def addSample(self, chrom, diap, letters):
        if self.mPrevIdx is not None:
            if self.mSamples[self.mPrevIdx][0] == chrom:
                add_letters = letters[:5]
                self.mSamples[self.mPrevIdx][1][1] += len(add_letters)
                self.mSamples[self.mPrevIdx][2] += add_letters
                self.mPrevIdx = None
        self.mTotal += 1
        if self.mTotal <= self.mSmpCount:
            idx = len(self.mSamples)
            self.mSamples.append(None)
        else:
            idx = self.mRH.randrange(0, self.mTotal - 1)
            if idx >= self.mSmpCount:
                return
        if self.mRH.choice([True, True, False]):
            pos = self.mRH.randrange(0, len(letters) - 1)
        else:
            pos = max(0, len(letters) - 3)
        test_letters = letters[pos:pos + 5]
        self.mSamples[idx] = [chrom,
            [diap[0] + pos, diap[0] + pos + len(test_letters) - 1],
            test_letters]
        if len(test_letters) < 5:
            self.mPrevIdx = idx

    def iterSamples(self):
        return iter(self.mSamples)

    def getSmpCount(self):
        return len(self.mSamples)
