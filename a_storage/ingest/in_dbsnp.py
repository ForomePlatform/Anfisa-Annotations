import sys, logging, gzip

from .in_util import TimeReport, extendFileList, dumpReader
#========================================
# Schema for AStorage
#========================================
SCHEMA_dbSNP = {
    "name": "dbSNP",
    "key": "hg38",
    "io": {
        "block-type": "page-cluster",
        "max-var-count": 500
    },
    "top": {
        "tp": "list",
        "item": {
            "tp": "dict",
            "label": "db-snp-rec",
            "items": [
                {"name": "ALT", "tp": "str", "opt": "gene"},
                {"name": "REF", "tp": "str", "opt": "gene"},
                {"name": "rs_id",  "tp": "str"}
            ]
        }
    }
}

#========================================
# Ingest logic
#========================================
DB_FIELDS = ["ALT", "REF", "rs_id"]
def prepare_data(fields):
    chrom, pos, rs_id, ref, alt = fields[:5]
    return ("chr" + chrom, int(pos)), [alt, ref, rs_id]

#========================================
class ReaderSNP:
    def __init__(self, file_list):
        self.mFiles = extendFileList(file_list)
        self.mCurRec = None

    def read(self):
        total = 0
        for file_name in self.mFiles:
            logging.info("Evaluation of %s" % file_name)
            assert file_name.endswith(".gz")
            with gzip.open(file_name, 'rt') as reader:
                time_rep = TimeReport(file_name)
                for line in reader:
                    if line.startswith('#'):
                        if not line.startswith('##'):
                            assert line.startswith(
                                "#CHROM\tPOS\tID\tREF\tALT\t")
                        continue
                    key, sub_rec = prepare_data(line.split('\t'))
                    if self.mCurRec and self.mCurRec[0] == key:
                        self.mCurRec[1].append(sub_rec)
                        continue
                    if self.mCurRec:
                        yield self.mCurRec
                        total += 1
                        if total % 100000 == 0:
                            time_rep.portion(total, self.mCurRec[0])
                    self.mCurRec = [key, [sub_rec]]
                if self.mCurRec:
                    yield self.mCurRec
                    self.mCurRec = None
                time_rep.done(total)

#========================================
def reader_SNP(properties, schema_h = None):
    global DB_FIELDS
    if schema_h is not None:
        schema_h.getCodecByLabel("db-snp-rec").setSerialization(DB_FIELDS)
    return ReaderSNP(properties["file_list"])


#========================================
if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)
    reader = reader_SNP({"file_list": sys.argv[1]})
    dumpReader(reader)
