import sys, re, os
from argparse import ArgumentParser

#=====================================
def parseRev(fname):
    rev_map = dict()
    with open(fname, "r") as inp:
        for line in inp:
            names = line.strip().split()
            rev_map[names[0]] = names[1:]
    return rev_map

sSymPattern = re.compile(r'^[A-Za-z0-9\-]*$')

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
                if not sSymPattern.match(sym):
                    bad_names.append(sym)
                else:
                    symbols.add(sym)
    return header, symbols, bad_names

#=====================================
parser = ArgumentParser()
parser.add_argument("-r", "--rev", help = "review file")
parser.add_argument("panel", nargs = 1,
    help= "panel file")
parser.add_argument("-C", "--calm",   action = "store_true",
    help = "No output result pannel")

run_args = parser.parse_args()
#=====================================
rev_fname = run_args.rev
if rev_fname is None:
    rev_fname = os.path.dirname(os.path.abspath(sys.argv[0])) + "/gtf_rev.txt"

rev_map = parseRev(rev_fname)

print("Processing", run_args.panel[0], file = sys.stderr)
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
    print("MAP:", sym, "->", ", ".join(rec))
    genes.update(rec)

for gene in sorted(genes):
    print(gene)

print("===================", file = sys.stderr)
