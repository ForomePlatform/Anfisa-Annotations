import logging
from .a_blocker import ABlockerIO_Complex
from .a_blk_pager import ABlockerIO_PosPager
#===============================================
DEBUG_MODE = False

#===============================================
class ABlockerIO_PageCluster(ABlockerIO_Complex):
    MAX_POS_COUNT = 50

    def __init__(self, schema, properties, key_codec_type):
        ABlockerIO_Complex.__init__(self, schema, properties, key_codec_type,
            use_cache = True, conv_bytes = False, pre_part = "PosSeq")
        self.mMaxVarCount = self._getProperty("max-var-count")
        self.mLastWriteKey = None
        self.mMaxPosCount = self._getProperty(
            "max-loc-count", self.MAX_POS_COUNT)
        stat_info = self._getProperty("stat", dict())
        if self.isWriteMode():
            self.mCountBlocks = stat_info.get("cluster-blocks", 0)
            self.mCountVariants = stat_info.get("custer-variants", 0)
        pager_properties = self._getProperty("pager", dict())
        if "block-type" not in pager_properties:
            pager_properties["block-type"] = "page-cluster"
        self.mPagerIO = ABlockerIO_PosPager(
            schema, pager_properties, key_codec_type)
        self._onDuty()

    def writeCountsAreGood(self, pos_count, var_count):
        if pos_count + 1 >= self.mMaxPosCount:
            return False
        return self.mMaxVarCount is None or var_count + 1 < self.mMaxVarCount

    def getBlockType(self):
        return "page-cluster"

    def updateWStat(self):
        ABlockerIO_Complex.updateWStat(self)
        if self.isWriteMode():
            stat_info = self._getProperty("stat")
            stat_info["cluster-blocks"] = self.mCountBlocks
            stat_info["cluster-variants"] = self.mCountVariants
        self.mPagerIO.updateWStat()

    def flush(self):
        self.mPagerIO.flush()
        ABlockerIO_Complex.flush(self)

    def close(self):
        self.mPagerIO.close()
        ABlockerIO_Complex.close(self)

    def putBlock(self, chrom, pos_seq, main_data_seq, var_count):
        key = (chrom, pos_seq[-1])
        self.mPagerIO.regPos(key)
        self._putBlock(key, [pos_seq] + main_data_seq)
        self.mCountBlocks += 1
        self.mCountVariants += var_count

    def openWriteBlock(self, key):
        if (self.mLastWriteKey is not None
                and key[0] == self.mLastWriteKey[0]):
            assert self.mLastWriteKey[1] < key[1]
        self.mLastWriteKey = key
        return _WriteBlock_PageCluster(self, key, self.mMaxVarCount is None,
            self.mPagerIO.getFinalPos(key))

    def openReadBlock(self, key, last_pos = None):
        assert last_pos is None
        chrom = key[0]
        key_base, pos_diap = self.mPagerIO.seekPos(key)
        if key_base is not None:
            seq_data = self._getBlock(key_base)
            pos_seq = [pos + key_base[1] for pos in seq_data[0]]
            return _ReadBlock_PageCluster(self, chrom, pos_diap,
                pos_seq, seq_data[1:])
        return _ReadBlock_PageCluster(self, chrom, pos_diap)

#===============================================
class _WriteBlock_PageCluster:
    def __init__(self, blocker, key, no_variants, final_pos):
        self.mBlocker = blocker
        self.mEncodeEnv = self.mBlocker.getSchema().makeDataEncoder()
        self.mChrom = key[0]
        self.mPosSeq = []
        self.mVarCount = 0
        self.mNoVariants = no_variants
        self.mFinalPos = final_pos

    def goodToAdd(self, key):
        chrom, pos = key
        if self.mChrom != chrom or pos >= self.mFinalPos:
            return False
        if self.mBlocker.writeCountsAreGood(len(self.mPosSeq), self.mVarCount):
            return len(self.mPosSeq) == 0 or self.mPosSeq[-1] + 0x100 >= pos
        return False

    def addRecord(self, key, record):
        chrom, pos = key
        assert self.mChrom == chrom and (len(self.mPosSeq) == 0
            or 0 < pos - self.mPosSeq[-1] <= 0x100)
        self.mPosSeq.append(pos)
        self.mEncodeEnv.put(record)
        if self.mNoVariants:
            self.mVarCount += 1
        else:
            self.mVarCount += len(record)

    def finishUp(self):
        if len(self.mPosSeq) > 0:
            self.mBlocker.putBlock(self.mChrom, self.mPosSeq,
                self.mEncodeEnv.result(), self.mVarCount)
        del self.mEncodeEnv

#===============================================
class _ReadBlock_PageCluster:
    def __init__(self, blocker, chrom, pos_diap,
            pos_seq = None, data_seq = None):
        self.mBlocker = blocker
        self.mChrom = chrom
        self.mStartPos, self.mEndPos = pos_diap
        self.mPosSeq = pos_seq
        self.mDataSeq = (self.mBlocker.getSchema().decodeData(data_seq)
            if data_seq is not None else None)

    def goodToRead(self, key, last_pos = None):
        assert last_pos is None
        chrom, pos = key
        return chrom == self.mChrom and self.mStartPos <= pos < self.mEndPos

    def getRecord(self, key, last_pos = None):
        assert last_pos is None
        chrom, pos = key
        assert self.mChrom == chrom
        if DEBUG_MODE:
            logging.info("For key = %s block from %d-%d size = %s/%s" % (
                str(key), self.mStartPos, self.mEndPos,
                str(len(self.mDataSeq))
                if self.mDataSeq is not None else "None",
                str(len(self.mPosSeq)
                if self.mPosSeq is not None else "None")))
        if self.mDataSeq is not None and pos in self.mPosSeq:
            return self.mDataSeq.get(self.mPosSeq.index(pos))
        return None
