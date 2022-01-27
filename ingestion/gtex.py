#
# Creates tables GENE, TISSUE, GENE2TISSUE
#       in db gtex and inserts data into them
#

import gzip, time
from io import TextIOWrapper
import mysql.connector
from collections import defaultdict

from util import execute_insert, reportTime

#=== table GENE ============

INSTR_CREATE_GENE = """CREATE TABLE IF NOT EXISTS GTex_GENE(
    GeneId         VARCHAR(24),
    SubId          VARCHAR(19),
    Symbol         VARCHAR(19),
    TopT1          INT(2),
    TopT2          INT(2),
    TopT3          INT(2),
    KEY `Symbol` (Symbol),
    PRIMARY KEY(GeneId, SubId));"""

COLUMNS_GENE = [
    "GeneId",
    "SubId",
    "Symbol",
    "TopT1",
    "TopT2",
    "TopT3"]

#=== table TISSUE ============
INSTR_CREATE_TISSUE = """CREATE TABLE IF NOT EXISTS GTex_TISSUE(
    TissueNo        INT(2),
    Name            VARCHAR(41),
    PRIMARY KEY(TissueNo));"""

COLUMNS_TISSUE = [
    "TissueNo",
    "Name"]

#=== table GENE2TISSUE ============
INSTR_CREATE_GENE2TISSUE = """CREATE TABLE IF NOT EXISTS GTex_GENE2TISSUE(
    GeneId           VARCHAR(24),
    SubId            VARCHAR(19),
    TissueNo         INT(2),
    Expression       FLOAT,
    RelExp           FLOAT,
    PRIMARY KEY(GeneId, SubId, TissueNo));"""

COLUMNS_GENE2TISSUE = [
    "GeneId",
    "SubId",
    "TissueNo",
    "Expression",
    "RelExp"]

#================================================
INSTR_INSERT_GENE = "INSERT INTO GTex_GENE (%s) VALUES (%s)" % (
    ", ".join(COLUMNS_GENE),
    ", ".join(['%s' for _ in COLUMNS_GENE]))

INSTR_INSERT_TISSUE = "INSERT INTO GTex_TISSUE (%s) VALUES (%s)" % (
    ", ".join(COLUMNS_TISSUE),
    ", ".join(['%s' for _ in COLUMNS_TISSUE]))

INSTR_INSERT_GENE2TISSUE = "INSERT INTO GTex_GENE2TISSUE (%s) VALUES (%s)" % (
    ", ".join(COLUMNS_GENE2TISSUE),
    ", ".join(['%s' for _ in COLUMNS_GENE2TISSUE]))

#========================================
def parseGeneId(gene_id_ext):
    gene_id, _, sub_id = gene_id_ext.partition('.')
    return gene_id, sub_id if sub_id else ''

#========================================
def fillTissues(conn, fields):
    tissues = [(idx + 1, fields[idx])
        for idx in range(2, len(fields))]
    c = conn.cursor()
    c.executemany(INSTR_INSERT_TISSUE, tissues)
    print('TISSUE is done')

#========================================
def new_gene_record(rec_keys, expr_sheet):
    rec = rec_keys[:]
    for neg_expr, tissue_no in expr_sheet[:3]:
        if neg_expr >= 0:
            break
        rec.append(tissue_no)
    while len(rec) < len(rec_keys) + 3:
        rec.append(None)
    return rec

#========================================
class GeneIdRegistry:
    def __init__(self):
        self.mReg = defaultdict(list)

    def regOne(self, gene_id, sub_id):
        self.mReg[gene_id].append(sub_id)

    def checkAll(self):
        cnt_multiple = 0
        max_multiplication = 0
        for gene_id in sorted(self.mReg.keys()):
            sub_id_seq = self.mReg[gene_id]
            assert sub_id_seq[0].isdigit, (
                "Not a digital primary sub_id %s for %s"
                % (sub_id_seq[0], gene_id))
            if len(sub_id_seq) > 1:
                sub_id_prefix = sub_id_seq[0] + '_'
                assert all(sub_id.startswith(sub_id_prefix)
                    for sub_id in sub_id_seq[1:]), (
                    "Too complex sub_id seq for %s: %s"
                    % (gene_id, str(sub_id_seq)))
                cnt_multiple += 1
                max_multiplication = max(
                    max_multiplication, len(sub_id_seq))
        print("Multiplication for sub_id: count = %d, max = %d"
            % (cnt_multiple, max_multiplication))

#========================================
def ingestGTEX(db_host, db_port, user, password, database,
        batch_size, filename):

    conn = mysql.connector.connect(
        host = db_host,
        port = db_port,
        user = user,
        password = password,
        database = database,
        autocommit = True)
    assert conn.is_connected(), "Failed to connect to DB"
    print('Connected to %s...' % database)

    curs = conn.cursor()
    print(INSTR_CREATE_GENE)
    curs.execute(INSTR_CREATE_GENE)
    print(INSTR_CREATE_TISSUE)
    curs.execute(INSTR_CREATE_TISSUE)
    print(INSTR_CREATE_GENE2TISSUE)
    curs.execute(INSTR_CREATE_GENE2TISSUE)

    id_registry = GeneIdRegistry()

    with gzip.open(filename, 'rb') as inp:
        gene_records, gen2tis_records = [], []
        total_gene, total_gen2tis = 0, 0
        start_time = time.time()
        text_inp = TextIOWrapper(inp,
            encoding = "utf-8", line_buffering = True)
        line_no = 0
        for line in text_inp:
            line_no += 1
            if line_no < 3:
                continue
            fields = line.rstrip().split('\t')
            if line_no == 3:
                fillTissues(conn, fields)
                continue
            expr_sheet = sorted((-float(fields[idx]), idx + 1)
                for idx in range(2, len(fields)))
            top_expr = - expr_sheet[0][0]
            if top_expr <= 0:
                continue

            gene_id, sub_id = parseGeneId(fields[0])
            symbol_name = fields[1]
            id_registry.regOne(gene_id, sub_id)

            gene_records.append(new_gene_record(
                [gene_id, sub_id, symbol_name], expr_sheet))
            if len(gene_records) >= batch_size:
                total_gene += execute_insert(conn,
                    INSTR_INSERT_GENE, gene_records)
                gene_records = []
                reportTime("Records_GENE", total_gene, start_time)

            for neg_expr, tissue_no in expr_sheet:
                if neg_expr >= 0:
                    break
                expr = - neg_expr
                gen2tis_records.append([gene_id, sub_id, tissue_no,
                    expr, expr / top_expr])
            if len(gen2tis_records) >= batch_size:
                total_gen2tis += execute_insert(
                    conn, INSTR_INSERT_GENE2TISSUE, gen2tis_records)
                gen2tis_records = []
                reportTime("Records_GTexGENE2TISSUE",
                    total_gen2tis, start_time)

        if len(gene_records) >= 0:
            total_gene += execute_insert(conn,
                INSTR_INSERT_GENE, gene_records)
            reportTime("Done_GENE", total_gene, start_time)
        if len(gen2tis_records) >= 0:
            total_gen2tis += execute_insert(conn,
                INSTR_INSERT_GENE2TISSUE, gen2tis_records)
            reportTime("Done_GTexGENE2TISSUE", total_gen2tis, start_time)
    curs.close()
    conn.close()
    id_registry.checkAll()


#========================================
if __name__ == '__main__':
    ingestGTEX(
        db_host  = "localhost",
        db_port  = 3306,
        user     = 'test',
        password = 'test',
        database = "gtex",
        batch_size = 10000,
        filename = "/db/data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9"
        + "_gene_median_tpm.gct.gz")
