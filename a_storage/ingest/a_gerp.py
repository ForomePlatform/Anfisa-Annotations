import sys, csv, logging
from fastnumbers import fast_real

from .a_util import TimeReport, detectFileChrom, extendFileList, dumpReader
#========================================
DB_FIELDS = ["GerpN", "GerpRS"]
def new_record(chrom, pos, lst):
    record = [fast_real(item) if item != "NA" else None
        for item in lst]
    return [("chr" + chrom, pos), record]

#========================================
class ReaderGerp:
    def __init__(self, file_list, chrom_loc = "chr"):
        self.mFiles = extendFileList(file_list)
        self.mChromLoc = chrom_loc

    def read(self):
        for chrom_file in self.mFiles:
            chrom = detectFileChrom(chrom_file, self.mChromLoc)
            logging.info("Evaluation of %s in %s"
                % (chrom, chrom_file))
            with open(chrom_file, 'r') as header:
                time_rep = TimeReport("chr" + chrom)
                position = 0
                reader = csv.reader(header, delimiter = '\t')
                for record in reader:
                    position += 1
                    yield new_record(chrom, position, record)
                    if position % 100000 == 0:
                        time_rep.portion(position)
                time_rep.done(position)

#========================================
def reader_Gerp(properties, schema_h = None):
    global DB_FIELDS
    if schema_h is not None:
        schema_h.getCodecByLabel("gerp-rec").setSerialization(DB_FIELDS)
    return ReaderGerp(
        properties["file_list"],
        properties.get("chrom_loc", "chr"))


#========================================
if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)
    reader = reader_Gerp({"file_list": sys.argv[1]})
    dumpReader(reader)
