import sys, re, json, csv
from argparse import ArgumentParser

#=====================================
def makeRev(db_name, rev_name):
    rev_map = dict()
    with open(db_name, "r", encoding = "utf-8") as inp:
        line_no = 0
        for line in inp:
            line_no += 1
            if line_no == 1:
                continue
            info = json.loads(line)
            refs = info.get("gtf-refs")
            if refs is not None:
                rev_map[info["_id"]] = refs
    with open(rev_name, "w", encoding = "utf-8") as outp:
        for nm in sorted(rev_map.keys()):
            print("\t".join([nm] + rev_map[nm]), file = outp)
    print(f"File {rev_name} is created with {len(rev_map)} records",
        file = sys.stderr)
    return rev_map

#=====================================
def parseRev(fname):
    rev_map = dict()
    with open(fname, "r", encoding = "utf-8") as inp:
        for line in inp:
            names = line.strip().split()
            rev_map[names[0]] = names[1:]
    return rev_map

#=====================================
sSymPattern = re.compile(r'^[A-Za-z0-9\-]*$')

def symbolIsOK(name):
    return sSymPattern.match(name)

#=====================================
def parseList(fname):
    header = []
    symbols = set()
    bad_names = []
    with open(fname, "r") as inp:
        for line in inp:
            if line.startswith('#'):
                header.append(line.rstrip())
                continue
            for sym in line.partition('#')[0].replace(',', ' ').split():
                if symbolIsOK(sym):
                    symbols.add(sym)
                else:
                    bad_names.append(sym)
    return header, symbols, bad_names

#=====================================
def parseCSV(fname):
    symbols = set()
    bad_names = []
    with open(fname, "r", encoding = "utf-8") as inp:
        rd = csv.DictReader(inp)
        for row in rd:
            sym = row["ensembl_gene_name"]
            if symbolIsOK(sym):
                symbols.add(sym)
            else:
                bad_names.append(sym)
    return [], symbols, bad_names

#=====================================
parser = ArgumentParser()
parser.add_argument("-d", "--db",
    help = "full database file, do not set if review file is ready")
parser.add_argument("-r", "--rev", default = "gene_db.rev",
    help = "review file (file stored if --db is set)")
parser.add_argument("-C", "--calm",   action = "store_true",
    help = "No output result pannel")
parser.add_argument("panel", nargs = 1,
    help= "panel file")

run_args = parser.parse_args()
#=====================================
if run_args.db:
    rev_map = makeRev(run_args.db, run_args.rev)
else:
    rev_map = parseRev(run_args.rev)

panel_name = run_args.panel[0]
print("Processing", panel_name, file = sys.stderr)

if panel_name.endswith(".csv"):
    header, symbols, bad_names = parseCSV(run_args.panel[0])
else:
    header, symbols, bad_names = parseList(run_args.panel[0])

for name in bad_names:
    print("ERROR: bad name", name, file = sys.stderr)

if not run_args.calm:
    for head in header:
        print(head)

genes = set()

for sym in sorted(symbols):
    rec = rev_map.get(sym)
    if rec is None:
        print("ERROR: no record for", sym, file = sys.stderr)
        continue
    if len(rec) == 1 and sym == rec[0]:
        genes.add(sym)
        continue
    if len(rec) == 0:
        print("MSG: no gtf for", sym, file = sys.stderr)
        continue
    print("MAP:", sym, "->", ", ".join(rec), file = sys.stderr)
    genes.update(rec)

if not run_args.calm:
    for gene in sorted(genes):
        print(gene)

print("===================", file = sys.stderr)
