# INSERTING tables ClinVar2Sub_Sig, CV_Submitters,
#           variant_summary, MySQL DATABASE clinvar
#
import csv
import gzip
import mysql.connector
import time
import sys

from util import execute_insert, reportTime

from parser_xml import XML_File

csv.field_size_limit(sys.maxsize)
#============= XML table ================


TABLE_META      = "Metadata"

COLUMN_SUB      = ("SubmitterID", "SubmitterName")
TYPE_SUB        = ('INT', 'TEXT CHARACTER SET utf8')
COLUMN_SIG      = ("SubmitterID", "RCVaccession", 'ClinicalSignificance')
TYPE_SIG        = ('INT', 'varchar(12)', 'varchar(60)')
COLUMN_META     = ("XML_Date", "CSV_upload_Date", "count")
TYPE_META       = ('TEXT DEFAULT NULL', 'TEXT DEFAULT NULL', 'INT(1)')
INDEX_SUB       = [("Index0", ["SubmitterID"], True, False)]
INDEX_SIG       = (("IndexUnic",
    ("RCVaccession", "SubmitterID", "ClinicalSignificance"), True, False),
    ("Index0", ("RCVaccession", "SubmitterID"), False, False),
    ("Index1", ("SubmitterID", "RCVaccession"), False, False))
INDEX_META      = tuple()

REFERENCE_PATHS = [
    ('ClinVarAccession', 'Acc')]  # RCVaccession

ASSERTION_PATHS = [
    ('ClinVarAccession', 'OrgID'),  # SubmitterID
    ('ClinVarSubmissionID', 'submitter'),  # SubmitterName
    'ClinicalSignificance/Description']  # ClinicalSignificance

TABLES_MAP = [((1, 0), (1, 1)), ((1, 0), (0, 0), (1, 2))]

INSTR_CREATE_SIG = """CREATE TABLE IF NOT EXISTS ClinVar2Sub_Sig(
    SubmitterID                 INT,
    RCVaccession                varchar(12),
    ClinicalSignificance        varchar(60),
    INDEX IndexUnic(RCVaccession, SubmitterID, ClinicalSignificance),
    INDEX Index0(RCVaccession, SubmitterID),
    INDEX Index1(SubmitterID, RCVaccession));"""

INSTR_INSERT_SIG = "INSERT INTO ClinVar2Sub_Sig(%s) VALUES(%s)" % (
    ", ".join(COLUMN_SIG),
    ", ".join(["%s" for _ in COLUMN_SIG]))

INSTR_CREATE_SUB = """CREATE TABLE IF NOT EXISTS ClinVar_Submitters(
    SubmitterID                 INT,
    SubmitterName               TEXT CHARACTER SET utf8,
    INDEX Index0(SubmitterID));"""

INSTR_INSERT_SUB = "INSERT INTO ClinVar_Submitters(%s) VALUES(%s)" % (
    ", ".join(COLUMN_SUB),
    ", ".join(["%s" for _ in COLUMN_SUB]))

#========================================
#--- table variant_summary ----------------

# capital letters?
INSTR_CREATE = """CREATE TABLE IF NOT EXISTS ClinVar_variant_summary(
    AlleleID                 INT(11),
    Type                     TEXT,
    Name                     TEXT,
    GeneID                   INT(11),
    GeneSymbol               TEXT,
    HGNC_ID                  TEXT,
    ClinicalSignificance     TEXT,
    ClinSigSimple            INT(11),
    LastEvaluated            TEXT,
    RS                       INT(11),
    nsv_esv                  TEXT,
    RCVaccession             TEXT,
    PhenotypeIDS             TEXT,
    PhenotypeList            TEXT,
    Origin                   TEXT,
    OriginSimple             TEXT,
    Assembly                 TEXT,
    ChromosomeAccession      TEXT,
    Chromosome               VARCHAR(12),
    Start                    INT(11),
    Stop                     INT(11),
    ReferenceAllele          TEXT,
    AlternateAllele          TEXT,
    Cytogenetic              TEXT,
    ReviewStatus             TEXT,
    NumberSubmitters         INT(11),
    Guidelines               TEXT,
    TestedInGTR              TEXT,
    OtherIDs                 TEXT,
    SubmitterCategories      INT(11),
    VariationID              INT(11),
    PositionVCF              INT(11),
    ReferenceAlleleVCF       LONGTEXT,
    AlternateAlleleVCF       LONGTEXT,
    INDEX index1(AlleleID),
    INDEX c_idx(Chromosome),
    INDEX p_idx(Start, Stop),
    INDEX Cs0_idx(ClinSigSimple));"""

COLUMNS = [
    "AlleleID",
    "Type",
    "Name",
    "GeneID",
    "GeneSymbol",
    "HGNC_ID",
    "ClinicalSignificance",
    "ClinSigSimple",
    "LastEvaluated",
    "RS",
    "nsv_esv",
    "RCVaccession",
    "PhenotypeIDS",
    "PhenotypeList",
    "Origin",
    "OriginSimple",
    "Assembly",
    "ChromosomeAccession",
    "Chromosome",
    "Start",
    "Stop",
    "ReferenceAllele",
    "AlternateAllele",
    "Cytogenetic",
    "ReviewStatus",
    "NumberSubmitters",
    "Guidelines",
    "TestedInGTR",
    "OtherIDs",
    "SubmitterCategories",
    "VariationID",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF"]

INSTR_INSERT = "INSERT INTO ClinVar_variant_summary (%s) VALUES (%s)" % (
    ", ".join(COLUMNS),
    ", ".join(["%s" for _ in COLUMNS]))

#========================================
def new_record(lst):
    rec = []
    for item in lst:
        if item == "NA":
            rec.append(None)
        else:
            rec.append(item)
    return rec

#========================================
def ingestCLINVAR(db_host, db_port, user, password, database,
        batch_size, summary_fname, xml_fname):

    conn = mysql.connector.connect(
        host = db_host,
        port = db_port,
        user = user,
        password = password,
        database = database,
        autocommit = True)
    assert conn.is_connected(), "Failed to connect to DB"
    print("Connected to %s..." % database)

    curs = conn.cursor()
    #========= XML tables ===========

    # Create/rename database
    xml = XML_File(file_name = xml_fname, ref_paths = REFERENCE_PATHS,
        asr_paths = ASSERTION_PATHS, tables_map = TABLES_MAP)
    batch_gen = xml.get_batch(batch_size)
    print(INSTR_CREATE_SIG)
    curs.execute(INSTR_CREATE_SIG)
    # Significance insertion, submitter accumulation
    sub_dict = {}  # unique sumitter dict
    ic = 0  # insert counter
    tb = time.time()  # rate time, batch rate time
    count = 0

    for batch_list in batch_gen:
        # Dict and list accumulation
        insert_list = []
        for d in batch_list:
            if d[0][0] not in sub_dict:
                sub_dict.update({d[0][0]: d[0][1]})
            elif sub_dict[d[0][0]] != d[0][1]:
                if sub_dict[d[0][0]] is None:
                    sub_dict[d[0][0]] = d[0][1]

            insert_list += [d[1]]

        # Data base insertion
        ic += execute_insert(conn, INSTR_INSERT_SIG, insert_list)
        count += len(insert_list)
        reportTime("", count, tb)
    if len(insert_list) > 0:
        ic += execute_insert(conn, INSTR_INSERT_SIG, insert_list)
        count += len(insert_list)
    reportTime("Done", count, tb)

    # Init database submitters table
    print(INSTR_CREATE_SUB)
    curs.execute(INSTR_CREATE_SUB)

    # Submitter insertion
    sub_list = []  # unique submitter list
    ic = 0  # insert counter
    tb = time.time()  # rate time, batch rate time
    for sid in sub_dict:
        sub_list += [(sid, sub_dict[sid])]
        if len(sub_list) >= batch_size:
            # Data base insertion
            ic += execute_insert(conn, INSTR_INSERT_SUB, sub_list)
            reportTime("", ic, tb)
            del sub_list[:]
            tb = time.time()
    if len(sub_list) > 0:
        ic += execute_insert(conn, INSTR_INSERT_SUB, sub_list)
    reportTime("Done", ic, tb)

    #========= end of XML ============
    print(INSTR_CREATE)
    curs.execute(INSTR_CREATE)
    print(summary_fname)
    with gzip.open(summary_fname, "rt") as header:
        position = 0
        list_of_records = []
        total, cnt = 0, 0
        reader = csv.reader(header, delimiter = "\t")
        start_time = time.time()
        for record in reader:
            position += 1
            if position == 1:
                continue
            list_of_records.append(
                new_record(record))
            if len(list_of_records) >= batch_size:
                total += execute_insert(conn, INSTR_INSERT,
                    list_of_records)
                list_of_records = []
                cnt += 1
                if cnt >= 100:
                    cnt = 0
                    reportTime("", total, start_time)
        if len(list_of_records) > 0:
            total += execute_insert(conn, INSTR_INSERT, list_of_records)
        reportTime("Done:", total, start_time)


#========================================
if __name__ == "__main__":
    ingestCLINVAR(
        db_host  = "localhost",
        db_port  = 3306,
        user     = "test",
        password = "test",
        database = "clinvar",
        batch_size = 1000,
        summary_fname = "/home/trosman/work/clinvar/variant_summary.txt.gz",
        xml_fname =
        "/home/trosman/work/clinvar/ClinVarFullRelease_2020-01.xml.gz")
