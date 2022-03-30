import sys, gzip, logging
from collections import Counter

from .in_util import TimeReport, detectFileChrom, extendFileList, dumpReader
#========================================
# dbNSDP fields
#========================================
class FieldH:
    sRolesLists = {role: [] for role in ("variant", "facet", "transcript")}
    sKeyLists = {role: [] for role in ("variant", "facet", "transcript")}

    @classmethod
    def getRoleList(cls, role):
        return cls.sRolesLists[role]

    @classmethod
    def getRoleKeys(cls, role):
        return cls.sKeyLists[role]

    @classmethod
    def setupIndexes(cls, fields_idxs):
        for fld_seq in cls.sRolesLists.values():
            for fld in fld_seq:
                fld._setDataIdx(fields_idxs)

    def __init__(self, role, name, tp_norm = str,
            opt = None, specific = None):
        self.mRole = role
        self.mName = name
        self.mTpNorm = tp_norm
        self.mOpt = opt
        self.mSpec = specific
        self.mDataIdx = None
        self.sRolesLists[self.mRole].append(self)
        if self.mSpec == "key":
            self.sKeyLists[self.mRole].append(self.mName)

    def getName(self):
        return self.mName

    def getRole(self):
        return self.mRole

    def getData(self, record):
        return record[self.mDataIdx]

    def normData(self, data):
        if data == '.':
            return None
        if data == '-' and self.mTpNorm is float:
            return None
        return self.mTpNorm(data)

    def fetchData(self, record):
        return self.normData(self.getData(record))

    def isSpec(self, spec):
        return self.mSpec == spec

    def _setDataIdx(self, fields_idxs):
        self.mDataIdx = fields_idxs[self.mName]

    def makeDescr(self):
        ret = {"name": self.mName}
        if self.mTpNorm == str:
            ret["tp"] = "str"
        elif self.mTpNorm == float or self.mTpNorm == int:
            ret["tp"] = "num"
        else:
            assert False, ("No correct type for %s" % self.mName)
        if self.mOpt is not None:
            ret["opt"] = self.mOpt
        return ret


#========================================
FieldH("variant", "ALT", opt = "gene", specific = "key")
FieldH("variant", "REF", opt = "gene", specific = "key")
FieldH("variant", "CADD_raw", float)
FieldH("variant", "CADD_phred", float)
FieldH("variant", "DANN_score", float)
FieldH("variant", "DANN_rankscore", float)
FieldH("variant", "Eigen_raw_coding", float)
FieldH("variant", "Eigen_raw_coding_rankscore", float)
FieldH("variant", "Eigen_phred_coding", float)
FieldH("variant", "Eigen_PC_raw_coding", float)
FieldH("variant", "Eigen_PC_raw_coding_rankscore", float)
FieldH("variant", "Eigen_PC_phred_coding", float)
FieldH("variant", "GTEx_V7_gene", opt = "repeat")
FieldH("variant", "GTEx_V7_tissue")
FieldH("variant", "Geuvadis_eQTL_target_gene")

FieldH("facet", "CADD_raw_rankscore", float)
FieldH("facet", "MetaLR_score", float)
FieldH("facet", "MetaLR_rankscore", float)
FieldH("facet", "MetaLR_pred", opt = "dict")
FieldH("facet", "MutPred_score", float)
FieldH("facet", "MutPred_rankscore", float)
FieldH("facet", "MutPred_protID")
FieldH("facet", "MutPred_AAchange")
FieldH("facet", "MutPred_Top5features")
FieldH("facet", "MPC_rankscore", float)
FieldH("facet", "PrimateAI_score", float)
FieldH("facet", "PrimateAI_rankscore", float)
FieldH("facet", "PrimateAI_pred", opt = "dict")
FieldH("facet", "REVEL_score", float)
FieldH("facet", "SIFT_converted_rankscore", float)
FieldH("facet", "SIFT4G_converted_rankscore", float)
FieldH("facet", "MutationTaster_score")
FieldH("facet", "MutationTaster_pred")

FieldH("transcript", "refcodon", opt = "repeat")
FieldH("transcript", "codonpos", int)

FieldH("transcript", "Ensembl_transcriptid", opt = "repeat", specific = "key")
FieldH("transcript", "Ensembl_geneid", opt = "repeat")
FieldH("transcript", "Ensembl_proteinid", opt = "repeat")
FieldH("transcript", "FATHMM_score", float)
FieldH("transcript", "FATHMM_pred", opt = "dict")
FieldH("transcript", "GENCODE_basic", opt = "dict")
FieldH("transcript", "HGVSc_ANNOVAR")
FieldH("transcript", "HGVSp_ANNOVAR")
FieldH("transcript", "HGVSc_snpEff")
FieldH("transcript", "HGVSp_snpEff")
FieldH("transcript", "MPC_score", float)
FieldH("transcript", "MutationAssessor_score", float)
FieldH("transcript", "MutationAssessor_pred", opt = "dict")
FieldH("transcript", "Polyphen2_HDIV_score", float)
FieldH("transcript", "Polyphen2_HDIV_pred", opt = "dict")
FieldH("transcript", "Polyphen2_HVAR_score", float)
FieldH("transcript", "Polyphen2_HVAR_pred", opt = "dict")
FieldH("transcript", "SIFT_score", float)
FieldH("transcript", "SIFT_pred", opt = "dict")
FieldH("transcript", "SIFT4G_score", float)
FieldH("transcript", "SIFT4G_pred", opt = "dict")
FieldH("transcript", "Uniprot_acc")

#========================================
# Schema for AStorage
#========================================
TRANSCRIPT_ITEMS = [fld.makeDescr()
    for fld in FieldH.getRoleList("transcript")]

FACET_ITEMS = [fld.makeDescr()
    for fld in FieldH.getRoleList("facet")] + [{
        "name": "transcripts",
        "tp": "list",
        "item": {
            "tp": "dict",
            "items": TRANSCRIPT_ITEMS
        }
    }]

VARIANT_ITEMS = [fld.makeDescr()
    for fld in FieldH.getRoleList("variant")] + [{
        "name": "facets",
        "tp": "list",
        "item": {
            "tp": "dict",
            "items": FACET_ITEMS
        }
    }]


SCHEMA_DBNSFP_4 = {
    "name": "DBNSFP",
    "key": "hg38",
    "io": {
        "block-type": "page-cluster",
        "max-var-count": 50
    },
    "filter-list": {"ref": "REF", "alt": "ALT"},
    "top": {
        "tp": "list",
        "item": {
            "tp": "dict",
            "items": VARIANT_ITEMS
        }
    }
}

#========================================
FLD_NAME_MAP = {
    "ref": "REF",
    "alt": "ALT",
    "Eigen_pred_coding": "Eigen_phred_coding"
}

def _normFieldName(name):
    global FLD_NAME_MAP
    name = name.replace('-', '_')
    return FLD_NAME_MAP.get(name, name)

def setupFields(field_line):
    global FLD_NAME_MAP
    assert field_line.startswith('#')
    field_names = field_line[1:].split()
    assert field_names[0].startswith("chr")
    assert field_names[1].startswith("pos")
    FieldH.setupIndexes({_normFieldName(name): idx
        for idx, name in enumerate(field_names)})


#========================================
assert len(FieldH.getRoleKeys("transcript")) == 1, (
    "Transcript keys:" + " ".join(FieldH.getRoleKeys("transcript")))

TRANSCRIPT_KEY = FieldH.getRoleKeys("transcript")[0]
#========================================
class DataCollector:
    def __init__(self):
        self.mCounts = [0, 0, 0]
        self.mCurRecord = None

    def getCounts(self):
        return self.mCounts

    def ingestLine(self, line):
        if line.endswith('\n'):
            line = line[:-1]
        fields = line.split('\t')
        chrom = "chr" + str(fields[0])
        pos = int(fields[1])
        new_record = False
        if self.mCurRecord is None or (chrom, pos) != self.mCurRecord[0]:
            new_record = True
        new_variant = new_record

        variant_data = {fld.getName(): fld.fetchData(fields)
            for fld in FieldH.getRoleList("variant")}

        if not new_variant:
            prev_var_data = self.mCurRecord[1][-1]
            new_variant = any(
                prev_var_data[name] != variant_data[name]
                for name in FieldH.getRoleKeys("variant"))
            if not new_variant:
                for fld in FieldH.getRoleList("variant"):
                    name = fld.getName()
                    assert prev_var_data[name] == variant_data[name], (
                        "Inequalty in variant field %s" % name)

        facet_data = {fld.getName(): fld.fetchData(fields)
            for fld in FieldH.getRoleList("facet")}

        transcript_collections = {fld.getName(): fld.getData(fields)
            for fld in FieldH.getRoleList("transcript")}

        transcript_data = [dict()
            for _ in transcript_collections[TRANSCRIPT_KEY].split(';')]

        for fld in FieldH.getRoleList("transcript"):
            name = fld.getName()
            values = transcript_collections[name].split(';')
            if len(values) != len(transcript_data):
                if fld.isSpec("same"):
                    v0 = values[0]
                    if len(values) > 0:
                        if any(vv != v0 for vv in values[1:]):
                            logging.error("Wrong tr field %s: %s/%d" % (name,
                                transcript_collections[name],
                                len(transcript_data)))
                else:
                    assert False, "Missing tr field %s" % name
            for idx, vv in enumerate(values):
                transcript_data[idx][name] = fld.normData(vv)
        facet_data["transcripts"] = transcript_data
        self.mCounts[2] += len(transcript_data)
        self.mCounts[1] += 1

        ret = None
        if new_record:
            self.mCounts[0] += 1
            variant_data["facets"] = [facet_data]
            self.checkRecord()
            ret, self.mCurRecord = self.mCurRecord, [
                (chrom, pos), [variant_data]]
        elif new_variant:
            self.mCounts[0] += 1
            variant_data["facets"] = [facet_data]
            self.mCurRecord[1].append(variant_data)
        else:
            self.mCurRecord[1][-1]["facets"].append(facet_data)

        return ret

    def finishUp(self):
        self.checkRecord()
        return self.mCurRecord

    def checkRecord(self):
        if not self.mCurRecord:
            return
        v_keys = set()
        for v_data in self.mCurRecord[1]:
            key = "|".join([v_data["REF"], v_data["ALT"]])
            assert key not in v_keys, ("Dup key: %s" % key)
            v_keys.add(key)
            marks = set()
            if len(v_data["facets"]) > 2:
                marks.add("multi-facet-%d" % len(v_data["facets"]))
            tr_counts = Counter()
            for f_data in v_data["facets"]:
                for t_data in f_data["transcripts"]:
                    tr_counts[t_data[TRANSCRIPT_KEY]] += 1
            dup_tr_id = []
            for tr_id, cnt in tr_counts.items():
                if cnt > 1:
                    dup_tr_id.append(tr_id)
            if len(dup_tr_id) > 0:
                marks.add("multi-transcript-%d" % max(tr_counts.values()))
                for tr_id in dup_tr_id:
                    self.fixTranscritDup(v_data, tr_id)
                #for idx in range(len(v_data["facets"])-1, -1, -1):
                #    f_data = v_data["facets"][idx]
                #    if len(f_data["transcripts"]) == 0:
                #        del v_data["facets"][idx]
            if len(marks) > 0:
                logging.info("Complications at %s|%s: %s"
                    % (str(self.mCurRecord[0]), key, " ".join(sorted(marks))))

    def fixTranscritDup(self, v_data, tr_id):
        tr_place_seq = []
        for f_idx, f_data in enumerate(v_data["facets"]):
            for t_idx, t_data in enumerate(f_data["transcripts"]):
                if t_data[TRANSCRIPT_KEY] == tr_id:
                    tr_place_seq.append([t_data, f_data, f_idx, t_idx])
        while len(tr_place_seq) > 1:
            t_data1 = tr_place_seq[0][0]
            t_data2 = tr_place_seq[-1][0]
            cnt1, cnt2 = 0, 0
            res_data = dict()
            for key, val1 in t_data1.items():
                val2 = t_data2[key]
                if val1 == val2:
                    res_data[key] = val1
                    continue
                if val1 is None:
                    res_data[key] = val2
                    cnt2 += 1
                    continue
                if val2 is None:
                    res_data[key] = val1
                    cnt1 += 1
                logging.info(
                    "Failed to fix multi-transcript for %s: %d/%d %d/%d"
                    % (tr_id, tr_place_seq[0][-2], tr_place_seq[0][-1],
                        tr_place_seq[-1][-2], tr_place_seq[-1][-1]))
                return
            if cnt2 > cnt1:
                t_data, f_data = tr_place_seq[-1][:2]
                self._removeTranscript(*(tr_place_seq[0][:2]))
                del tr_place_seq[0]
            else:
                t_data, f_data = tr_place_seq[0][:2]
                self._removeTranscript(*(tr_place_seq[-1][:2]))
                del tr_place_seq[-1]
            if any(val != res_data[key] for key, val in t_data.items()):
                logging.info("Joined multi-transcript: %s" % tr_id)
                f_data["transcripts"][
                    f_data["transcripts"].index(t_data)] = res_data

    @staticmethod
    def _removeTranscript(t_data, f_data):
        tr_seq = f_data["transcripts"]
        del tr_seq[tr_seq.index(t_data)]

#========================================
#========================================
class ReaderDBNSFP4:
    def __init__(self, file_list, chrom_loc = "chr"):
        self.mFiles = extendFileList(file_list)
        self.mChromLoc = chrom_loc

    def read(self):
        exceptions = 0
        for chrom_file in self.mFiles:
            chrom = detectFileChrom(chrom_file, self.mChromLoc)
            logging.info("Evaluation of %s in %s" % (chrom, chrom_file))
            with gzip.open(chrom_file, 'rt') as text_inp:
                time_rep = TimeReport("chr" + chrom)
                collector = DataCollector()
                for line_no, line in enumerate(text_inp):
                    if line_no == 0:
                        setupFields(line)
                        continue
                    try:
                        info = collector.ingestLine(line)
                        if info is not None:
                            yield info
                        if (line_no % 10000) == 0:
                            total_var, _, _ = collector.getCounts()
                            time_rep.portion(total_var)
                    except IndexError:
                        exceptions += 1
                info = collector.finishUp()
                if info:
                    yield info
                total_var, total_facets, total_tr = collector.getCounts()
                time_rep.done(total_var)
                logging.info("transcripts: %d, facets: %d, exceptions: %d"
                    % (total_tr, total_facets, exceptions))

#========================================
def reader_dbNSFP4(properties, schema_h = None):
    return ReaderDBNSFP4(
        properties["file_list"],
        properties.get("chrom_loc", "chr"))


#========================================
if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)
    reader = reader_dbNSFP4({"file_list": sys.argv[1]})
    dumpReader(reader, indent_mode = True)
