import time, re
import mysql.connector
from binascii import crc32

from util import execute_insert, reportTime

#========================================
class Parsing:
    # Some tools to parse strings

    sPattVar = re.compile(r"^rs\d+$")

    @classmethod
    def checkGoodVar(cls, variant):
        if not variant:
            return False
        if variant.startswith('rs'):
            assert cls.sPattVar.match(variant)
            return True
        return False

    @classmethod
    def splitValues(cls, value):
        if value.startswith('"'):
            assert value.endswith('"')
            return value[1:-1].split('","')
        return value.split(',')

    sPattTitleId = re.compile(r"^(.+)\((\w+)\)\s*$")
    @classmethod
    def parseTitleID(cls, value):
        q = cls.sPattTitleId.match(value)
        if q is not None:
            return [q.group(1).strip(), q.group(2)]
        title = value.strip()
        return [title, "__%d" % crc32(title.encode("utf-8"))]

#========================================
class TabWriter:
    sBaseColumns = ["Variant", "AssocKind", "AssocID"]

    def __init__(self, conn, table, spec_columns, instr_create,
            feed_mode = None):
        self.mConn = conn
        self.mTable = table
        self.mSpecColumns = spec_columns
        self.mFeedMode = feed_mode
        assert 1 <= len(self.mSpecColumns) <= 2

        self.mConn.cursor().execute(instr_create)

        columns = self.sBaseColumns + self.mSpecColumns
        self.mInstrInsert = ("INSERT INTO " + self.mTable + " "
            + "(" + ",".join(columns)+")"
            + " VALUES (" + ", ".join(['%s'] * len(columns)) + ");")
        self.mValues = []
        self.mTotal = 0

    def flush(self):
        if len(self.mValues) > 0:
            self.mTotal += execute_insert(self.mConn,
                self.mInstrInsert, self.mValues)
        self.mValues = []

    def close(self):
        assert len(self.mValues) == 0
        self.mValues = None
        return self.mTotal

    def feed(self, assoc_kind, assoc_id, variant, value):
        if value in ("", None):
            return
        assert Parsing.checkGoodVar(variant)
        if self.mFeedMode == "single":
            assert len(self.mSpecColumns) == 1
            self.mValues.append([variant, assoc_kind, assoc_id, value])
            return
        if self.mFeedMode == "int?":
            assert len(self.mSpecColumns) == 1
            for v in str(value).split(';'):
                self.mValues.append([variant, assoc_kind, assoc_id, v])
            return
        for val in Parsing.splitValues(value):
            val = val.strip()
            res_list = [variant, assoc_kind, assoc_id]
            if len(self.mSpecColumns) == 1:
                res_list.append(val.strip())
            else:
                res_list += Parsing.parseTitleID(val)
            self.mValues.append(res_list)

#========================================
class TabScanner:
    def __init__(self, assoc_kind, table, columns, tab_writers):
        self.mAssocKind = assoc_kind
        self.mTable = table
        self.mColumns = columns
        self.mTabWriters = tab_writers
        self.mCurId = None

    def scan(self, conn):
        while True:
            if self.scanPortion(conn) == 0:
                break

    def scanPortion(self, conn):
        instr = "SELECT " + ", ".join(self.mColumns) + " FROM " + self.mTable
        if self.mCurId is not None:
            instr += " WHERE " + self.mColumns[0] + " > " + str(self.mCurId)
        instr += " ORDER BY " + self.mColumns[0] + " LIMIT 1000;"
        cur = conn.cursor()
        cur.execute(instr)
        count = 0
        for fields in cur.fetchall():
            count += 1
            self.processOne(fields)
        cur.close()
        if count > 0:
            for writer in self.mTabWriters:
                writer.flush()
        return count

    def processOne(self, fields):
        self.mCurId, variant = fields[:2]
        if not Parsing.checkGoodVar(variant):
            return
        for idx, writer in enumerate(self.mTabWriters):
            writer.feed(self.mAssocKind, self.mCurId, variant, fields[2 + idx])

#========================================
def pgkbReTab(db_host, db_port, user, password, database):
    conn = mysql.connector.connect(
        host = db_host,
        port = db_port,
        user = user,
        password = password,
        database = database,
        autocommit = True)
    assert conn.is_connected(), "Failed to connect to DB"
    print('Connected to %s...' % database)

    start_time = time.time()

    tab_writers = [
        TabWriter(conn, "PharmCHEMICALS", ["ChTitle", "ChID"],
          instr_create = """CREATE TABLE IF NOT EXISTS PharmCHEMICALS(
            Variant VARCHAR(20),
            AssocKind VARCHAR(10),
            AssocID int(10) NOT NULL,
            ChTitle TEXT,
            ChID VARCHAR(20),
            PRIMARY KEY (AssocKind, AssocID, ChID),
            KEY (Variant));"""),
        TabWriter(conn, "PharmDISEASES", ["DisTitle", "DisID"],
          instr_create = """CREATE TABLE IF NOT EXISTS PharmDISEASES(
            Variant VARCHAR(20),
            AssocKind VARCHAR(10),
            AssocID int(10) NOT NULL,
            DisTitle TEXT,
            DisID VARCHAR(20),
            PRIMARY KEY (AssocKind, AssocID, DisID),
            KEY (Variant));"""),
        TabWriter(conn, "PharmPMIDS", ["PMID"],
          instr_create = """CREATE TABLE IF NOT EXISTS PharmPMIDS(
            Variant VARCHAR(20),
            AssocKind VARCHAR(10),
            AssocID int(10) NOT NULL,
            PMID int(10),
            PRIMARY KEY (AssocKind, AssocID, PMID),
            KEY (Variant));""", feed_mode = "int?"),
        TabWriter(conn, "PharmNOTES", ["Note"],
          instr_create = """CREATE TABLE IF NOT EXISTS PharmNOTES(
            Variant VARCHAR(20),
            AssocKind VARCHAR(10),
            AssocID int(10) NOT NULL,
            Note BLOB,
            PRIMARY KEY (AssocKind, AssocID),
            KEY (Variant));""", feed_mode = "single")]

    for scanner in (
            TabScanner("clinical", "PharmCAmeta",
                ["CAID", "LOC", "RC", "RD", "PMIDS", "AT"], tab_writers),
            TabScanner("drug", "PharmVDA",
                ["AID", "VAR", "CHEM", "NULL", "PMID", "NOTES"], tab_writers),
            TabScanner("fa", "PharmVFA",
                ["AID", "VAR", "CHEM", "NULL", "PMID", "NOTES"], tab_writers),
            TabScanner("pheno", "PharmVPA",
                ["AID", "VAR", "CHEM", "NULL", "PMID", "NOTES"], tab_writers)):
        scanner.scan(conn)

    total = 0
    for writer in tab_writers:
        total += writer.close()
    reportTime("Done:", total, start_time)


#========================================
if __name__ == '__main__':
    pgkbReTab(
        db_host = "localhost",
        db_port = 3306,
        user   = 'test',
        password = 'test',
        database = "pharmgkb")
