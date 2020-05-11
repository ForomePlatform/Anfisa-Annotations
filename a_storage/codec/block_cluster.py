#===============================================
class BlockerCluster():
    MAX_POS_COUNT = 50

    def __init__(self, master):
        self.mIO = master
        self.mMaxVarCount = self.mIO._getProperty("max-var-count")
        self.mCurWriterKey = None
        self.mIdxColumns = [self.mIO._regColumn("sgidx")]
        self.mMaxPosCount = self.mIO._getProperty(
            "max-loc-count", self.MAX_POS_COUNT)
        if self.mIO.isWriteMode():
            self.mCountBlocks = 0
            self.mCountVariants = 0

    def _addWriteStat(self, var_count):
        self.mCountBlocks += 1
        self.mCountVariants += var_count

    def writeCountsAreGood(self, pos_count, var_count):
        if pos_count + 1 >= self.mMaxPosCount:
            return False
        return self.mMaxVarCount is None or var_count + 1 < self.mMaxVarCount

    def getType(self):
        return "cluster"

    def getIO(self):
        return self.mIO

    def close(self):
        if self.mIO.isWriteMode():
            stat_info = {"cluster-blocks": self.mCountBlocks}
            if self.mMaxVarCount is not None:
                stat_info["cluster-variants"] = self.mCountVariants
            self.mIO._updateProperty("stat", stat_info)

    def addToIndex(self, chrom, pos_seq):
        self.mIO._putColumns((chrom, pos_seq[-1]),
            [pseq2bytes(pos_seq)], self.mIdxColumns, conv_bytes = False)

    def createWriteBlock(self, encode_env, key, codec):
        if (self.mCurWriterKey is not None
                and key[0] == self.mCurWriterKey[0]):
            assert self.mCurWriterKey[1] < key[1]
        self.mCurWriterKey = key
        return _WriteClusterBlock(self, encode_env, key,
            self.mMaxVarCount is None)

    def createReadBlock(self, decode_env_class, key, codec):
        key_base, data_base = self.mIO._seekColumn(key,
            self.mIdxColumns[0], conv_bytes = False)
        pos_seq = None
        if key_base is not None and key_base[0] == key[0]:
            pos_seq = bytes2pseq(key_base[1], data_base)
            data_seq = self.mIO._getColumns(key_base)
            return _ReadClusterBlock(self, key, pos_seq,
                decode_env_class(data_seq))
        return _ReadClusterBlock(self, key)

#===============================================
def pseq2bytes(pos_seq):
    xseq = [0, 0]
    for idx in range(len(pos_seq) - 1, 0, -1):
        delta = pos_seq[idx] - pos_seq[idx - 1] - 1
        if xseq[-1] == delta and xseq[-2] < 255:
            xseq[-2] += 1
        else:
            xseq += [1, delta]
    if len(xseq) > 2 and xseq[0] == 0:
        assert xseq[1] == 0
        del xseq[:2]
    return b''.join(n.to_bytes(1, 'big') for n in xseq)

def bytes2pseq(pos0, xbytes):
    pseq = [pos0]
    for idx in range(0, len(xbytes), 2):
        cnt, delta = xbytes[idx], xbytes[idx + 1]
        for _ in range(cnt):
            pseq.append(pseq[-1] - delta - 1)
    return pseq[-1::-1]

#===============================================
class _WriteClusterBlock:
    def __init__(self, blocker, encode_env, key, no_variants):
        self.mBlocker = blocker
        self.mEncodeEnv = encode_env
        self.mChrom = key[0]
        self.mPosSeq = []
        self.mVarCount = 0
        self.mNoVariants = no_variants

    def goodToWrite(self, key):
        chrom, pos = key
        if self.mChrom == chrom and self.mBlocker.writeCountsAreGood(
                len(self.mPosSeq), self.mVarCount):
            return len(self.mPosSeq) == 0 or self.mPosSeq[-1] + 0x100 >= pos
        return False

    def addRecord(self, key, record, codec):
        chrom, pos = key
        assert self.mChrom == chrom and (len(self.mPosSeq) == 0
            or 0 < pos - self.mPosSeq[-1] <= 0x100)
        self.mPosSeq.append(pos)
        self.mEncodeEnv.put(record, codec)
        if self.mNoVariants:
            self.mVarCount += 1
        else:
            self.mVarCount += len(record)

    def finishUp(self):
        if len(self.mPosSeq) > 0:
            self.mBlocker.getIO()._putColumns(
                (self.mChrom, self.mPosSeq[-1]), self.mEncodeEnv.result())
            self.mBlocker.addToIndex(self.mChrom, self.mPosSeq)
            self.mBlocker._addWriteStat(self.mVarCount)
        del self.mEncodeEnv

#===============================================
class _ReadClusterBlock:
    def __init__(self, blocker, seek_key, pos_seq = None, decode_env = None):
        self.mBlocker = blocker
        self.mChrom, self.mStartPos = seek_key
        self.mPosSeq = pos_seq
        self.mDecodeEnv = decode_env

    def goodToRead(self, key):
        chrom, pos = key
        if chrom != self.mChrom or pos < self.mStartPos:
            return False
        return self.mPosSeq is None or pos <= self.mPosSeq[-1]

    def getRecord(self, key, codec):
        chrom, pos = key
        assert self.mChrom == chrom
        if self.mDecodeEnv is not None and pos in self.mPosSeq:
            idx = self.mPosSeq.index(pos)
            return self.mDecodeEnv.get(idx, codec)
        return None
