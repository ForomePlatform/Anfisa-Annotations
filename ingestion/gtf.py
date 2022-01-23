#
# Creates tables GTF, GTF_genes
#       in db ensembl and fill them
#

import gzip, time
from io import TextIOWrapper
import mysql.connector

from util import execute_insert, reportTime

#=== table GTF ============

INSTR_CREATE_GTF = """CREATE TABLE IF NOT EXISTS GTF(
    chromosome  varchar(4) DEFAULT NULL,
    source      varchar(50) DEFAULT NULL,
    feature     varchar(50) DEFAULT NULL,
    start       int DEFAULT NULL,
    end         int DEFAULT NULL,
    score       float DEFAULT NULL,
    strand      varchar(1) DEFAULT NULL,
    frame       int DEFAULT NULL,
    attribute   varchar(2048) DEFAULT NULL,
    gene        varchar(64) DEFAULT NULL,
    biotype     varchar(64) DEFAULT NULL,
    exon        int DEFAULT NULL,
    transcript  varchar(64) DEFAULT NULL,
    KEY `PosIdx` (`chromosome`,`source`,`feature`,`start`,`end`,`strand`),
    KEY `SIdx` (`start`),
    KEY `EIdx` (`end`),
    KEY `FIdx` (`feature`),
    KEY `GIdx` (`gene`),
    KEY `TIdx` (`transcript`));"""

COLUMNS_GTF = [
    "chromosome",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
    "gene",
    "biotype",
    "exon",
    "transcript"]

INSTR_INSERT_GTF = "INSERT INTO GTF (%s) VALUES (%s)" % (
    ", ".join(COLUMNS_GTF),
    ", ".join(['%s'] * len(COLUMNS_GTF)))

#=== table GTF_gene ============

GENE_BUCKET_SIZE = 1000000

INSTR_CREATE_GTF_GENE = """CREATE TABLE IF NOT EXISTS GTF_gene(
    chromosome  varchar(4) DEFAULT NULL,
    start       int DEFAULT NULL,
    end         int DEFAULT NULL,
    gene        varchar(64) DEFAULT NULL,
    bucket      int DEFAULT NULL,
    KEY `MainIdx` (`gene`),
    KEY `SIdx` (`start`),
    KEY `EIdx` (`end`),
    KEY `bIdx` (`bucket`));"""


COLUMNS_GTF_GENE = [
    "chromosome",
    "start",
    "end",
    "gene",
    "bucket"]

INSTR_SELECT_GTF_GENE = """SELECT distinct `chromosome`, `start`, `end`, `gene`
    FROM GTF WHERE feature = 'gene' ORDER by 1, 2"""

INSTR_INSERT_GTF_GENE = "INSERT INTO GTF_gene (%s) VALUES (%s)" % (
    ", ".join(COLUMNS_GTF_GENE),
    ", ".join(['%s' for _ in COLUMNS_GTF_GENE]))

#========================================
def ingestGTF(db_host, db_port, user, password, database,
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

    with gzip.open(filename, 'rb') as inp:
        text_inp = TextIOWrapper(inp,
            encoding = "utf-8", line_buffering = True)
        loadGTF(text_inp, conn, batch_size)

    conn_read = mysql.connector.connect(
        host = db_host,
        port = db_port,
        user = user,
        password = password,
        database = database)

    loadGenes(conn_read, conn, batch_size)
    conn_read.close()
    conn.close()

#========================================
def parseAttributes(attr, keys):
    values = [None] * len(keys)
    for pair in attr.split(';'):
        seq = pair.strip().split(' ')
        if len(seq) < 2:
            continue
        if seq[0] not in keys:
            continue
        idx = keys.index(seq[0])
        values[idx] = seq[1].strip('"')
    return values

#========================================
def loadGTF(text_inp, conn, batch_size):
    curs = conn.cursor()
    print(INSTR_CREATE_GTF)
    curs.execute(INSTR_CREATE_GTF)
    curs.close()

    records = []
    total = 0
    start_time = time.time()
    line_no = 0
    for line in text_inp:
        line_no += 1
        if (line[0] == '#'):
            continue
        fields = [fld.strip() for fld in line.split('\t')]
        if len(fields[0]) > 4:
            continue
        chrom, source, feature = fields[:3]
        p_start, p_end = map(int, fields[3:5])
        score, strand, frame, attributes = fields[5:9]
        if score == ".":
            score = None
        if frame == ".":
            frame = None

        gene, biotype, exon, transcript = parseAttributes(attributes,
            ["gene_name", "gene_biotype", "exon_number", "transcript_id"])
        if exon is not None:
            exon = int(exon)

        records.append([chrom, source, feature, p_start, p_end,
            score, strand, frame, attributes,
            gene, biotype, exon, transcript])

        if len(records) >= batch_size:
            total += execute_insert(conn, INSTR_INSERT_GTF, records)
            records = []
            reportTime("Records_GTF", total, start_time)

    if len(records) >= 0:
        total += execute_insert(conn, INSTR_INSERT_GTF, records)
    reportTime("Done_GTF", total, start_time)

#========================================
def loadGenes(conn_read, conn, batch_size):
    curs = conn.cursor()
    print(INSTR_CREATE_GTF_GENE)
    curs.execute(INSTR_CREATE_GTF_GENE)
    curs.close()

    select_curs = conn_read.cursor()
    select_curs.execute(INSTR_SELECT_GTF_GENE)
    row = select_curs.fetchone()
    records = []
    total = 0
    start_time = time.time()
    while row:
        chrom, pos_start = row[0], int(row[2])
        bucket = int(pos_start / GENE_BUCKET_SIZE) * GENE_BUCKET_SIZE
        records.append([row[idx] for idx in range(4)] + [bucket])
        row = select_curs.fetchone()
        if len(records) >= batch_size:
            total += execute_insert(conn, INSTR_INSERT_GTF_GENE, records)
            records = []
            reportTime(f"Records_GTF_gene ({chrom}:{bucket}): ", total, start_time)
    select_curs.close()
    if len(records) >= batch_size:
        total += execute_insert(conn, INSTR_INSERT_GTF_GENE, records)
    reportTime("Done_GTF_gene", total, start_time)

#========================================
if __name__ == '__main__':
    ingestGTF(
        db_host  = "localhost",
        db_port  = 3306,
        user     = 'test',
        password = 'test',
        database = "ensembl",
        batch_size = 10000,
        filename = "/db/data/ensembl/Homo_sapiens.GRCh38.105.chr.gtf.gz")
