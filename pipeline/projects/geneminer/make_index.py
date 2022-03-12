import sys, json
from argparse import ArgumentParser
#=====================================
def loadGTFGenes(filename):
    gtf_genes = set()
    with open(filename, "r") as inp:
        for line in inp:
            chrom, gene = line.split()
            if chrom == "chromosome":
                continue
            gtf_genes.add(gene)
    return gtf_genes

#=====================================
def loadHGNC(filename, gtf_genes, all_names):
    with open(filename, "r") as inp:
        hgnc_data = json.loads(inp.read())
    hgnc_records = []
    for record in hgnc_data["response"]["docs"]:
        if record["locus_group"] != "protein-coding gene":
            continue
        sym = record["symbol"]
        alias_seq = record.get("alias_symbol", [])
        prev_seq = record.get("prev_symbol", [])
        names = set(alias_seq + prev_seq) | {sym}
        all_names.update(names)
        if sym in gtf_genes:
            hgnc_records.append([sym, names])
    return hgnc_records

#=====================================
def collectNames(name, gtf_genes, hgnc_records):
    ret = set()
    if name in gtf_genes:
        ret.add(name)
    for sym, names in hgnc_records:
        if name in names:
            ret.add(sym)
    return ret

#=====================================
parser = ArgumentParser()
parser.add_argument("--gtf",  default = "./gtf_genes.txt",
    help = "gtf_genes.txt")
parser.add_argument("--hgnc", default = "./hgnc_complete_set.json",
    help = "hgnc_complete_set.json")
run_args = parser.parse_args()

#=====================================
gtf_genes = loadGTFGenes(run_args.gtf)
all_names = gtf_genes.copy()
hgnc_records = loadHGNC(run_args.hgnc, gtf_genes, all_names)

cnt = 0
cnt_max = len(all_names)
for name in sorted(all_names):
    cnt += 1
    names = collectNames(name, gtf_genes, hgnc_records)
    print("\t".join([name] + sorted(names)))
    if cnt % 1000 == 0:
        print("Eval %d/%d (%s)" % (cnt, cnt_max, name), file = sys.stderr)
print("Done:", cnt_max, file = sys.stderr)
