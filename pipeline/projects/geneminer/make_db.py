#  Copyright (c) 2019. Partners HealthCare and other members of
#  Forome Association
#
#  Developed by Sergey Trifonov based on contributions by Joel Krier,
#  Michael Bouzinier, Shamil Sunyaev and other members of Division of
#  Genetics, Brigham and Women's Hospital
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

import sys, json, os
from argparse import ArgumentParser
from collections import defaultdict
from forome_tools.vcf import VCF_Support
#=====================================
parser = ArgumentParser()
parser.add_argument("-g", "--hgnc",
    help = "Path to hgnc_complete set (hgnc_complete_set_2022-04-01.json)")
parser.add_argument("-e", "--ensembl",
    help = "Path to ensembl archive (Homo_sapiens.GRCh38.105.chr.gtf)")
run_args = parser.parse_args()

#===============================================
def properGTFrecord(gene_name, rec_list):
    return gene_name != "" and (
        any (rec["gene_biotype"] == "protein_coding" for rec in rec_list))

def properHGNCrecord(record):
    return record.get("locus_group") == "protein-coding gene"


def versionOfFile(fpath):
    version = os.path.basename(fpath)
    while True:
        vhead, _, ext = version.rpartition('.')
        if ext in ("chr", "gtf", "gz", "json"):
            version = vhead
        else:
            break
    return version

#===============================================
sUniqueGeneFields = ["chrom"]
sOrdGeneFields = sUniqueGeneFields + [
    "gene_version", "gene_source", "gene_biotype"]
sAllGeneFields = sOrdGeneFields + ["gene_id", "gene_name"]

def prepareEnsembl(ensembl_file):
    symbols = defaultdict(list)
    gtf_records = dict()
    vcf_supp = VCF_Support(multi_fields = {"tag"})
    cnt = 0
    for rec in vcf_supp.readFile(ensembl_file):
        cnt += 1
        if cnt % 100000 == 0:
            print(f"GTF: {cnt} lines / {len(gtf_records)} records "
                + f"/ {len(symbols)} symbols...",
                end = '\r', file = sys.stderr)
        gene_name, gene_id = rec.get("gene_name"), rec["gene_id"]
        if gene_name is None:
            gene_name = ""
        grec = gtf_records.get(gene_id)
        if grec is None:
            grec = {fld: rec.get(fld, "") for fld in sAllGeneFields}
            gtf_records[gene_id] = grec
            symbols[gene_name].append(grec)
        else:
            for fld in sAllGeneFields:
                assert rec.get(fld, "") == grec[fld], (
                    f"GTF: at line {cnt} gene_name={gene_name} "
                    + f"conflict in {fld}: "
                    + f"{rec.get(fld)}|{grec[fld]}")

    print(f"GTF: {cnt} lines / {len(gtf_records)} records "
        + f"/ {len(symbols)} symbols", file = sys.stderr)

    meta_info = {"meta": [versionOfFile(ensembl_file)]}
    for gene_name in sorted(symbols.keys()):
        rec_list = symbols[gene_name]
        rec_list.sort(key = lambda rec: rec["gene_id"])
        if properGTFrecord(gene_name, rec_list):
            for rec in rec_list[1:]:
                for fld in sUniqueGeneFields:
                    assert rec[fld] == rec_list[0][fld], (
                        f"GTF: at line {cnt} "
                        + f"gene_name={rec['gene_name']}: "
                        + f"unique conflict in {fld}: "
                        + f"{rec.get(fld)}|{rec_list[0][fld]}")

    return meta_info, symbols, gtf_records

#===============================================
sHGNC_fields = ["symbol", "ensembl_gene_id", "gene_group",
    "locus_type", "locus_group", "alias_symbol", "prev_symbol"]

def loadHGNC(hgnc_file):
    with open(hgnc_file, "r") as inp:
        hgnc_data = json.loads(inp.read())
    hgnc_symbols = dict()
    print("HGNC load:", hgnc_data["response"].keys(), file = sys.stderr)
    for record in hgnc_data["response"]["docs"]:
        rec = dict()
        for fld in sHGNC_fields:
            if fld in record:
                rec[fld] = record[fld]
        gene_symbol = rec["symbol"]
        assert gene_symbol not in hgnc_symbols, (
            f"HGNC symbol duplication: {gene_symbol}")
        hgnc_symbols[gene_symbol] = rec
    print(f"HGNC loaded from {hgnc_file}: {len(hgnc_symbols)} records",
        file = sys.stderr)
    return hgnc_symbols, versionOfFile(hgnc_file)

#===============================================
def crossData(gtf_symbols, gtf_records, hgnc_symbols):
    gtf_genes = set()
    for gene_name, rec_list in gtf_symbols.items():
        if properGTFrecord(gene_name, rec_list):
            gtf_genes.add(gene_name)
    print(f"GTF filtered: {len(gtf_genes)}/{len(gtf_symbols)}",
        file = sys.stderr)

    gtf_refs = defaultdict(set)
    all_hgnc_names = set()
    cnt_proper_hgnc = 0
    for sym, record in hgnc_symbols.items():
        alias_seq = record.get("alias_symbol", [])
        prev_seq = record.get("prev_symbol", [])
        names = set(alias_seq + prev_seq) | {sym}
        all_hgnc_names.update(names)
        if properHGNCrecord(record):
            cnt_proper_hgnc += 1
            if sym in gtf_genes:
                for nm in names:
                    gtf_refs[nm].add(sym)

    print("HGNC filtered:", cnt_proper_hgnc, "/", len(hgnc_symbols),
        "/", len(gtf_refs), file = sys.stderr)

    extra_hgnc = all_hgnc_names - set(hgnc_symbols.keys())
    if len(extra_hgnc) > 0:
        print(f"HGNC extra names: {len(extra_hgnc)}",
            sorted(extra_hgnc)[:10], "...", file = sys.stderr)

    all_names = all_hgnc_names | set(gtf_symbols.keys())
    print(f"All names: {len(all_names)}", file = sys.stderr)

    records = []
    for gene_name in sorted(all_names):
        if gene_name == "":
            continue
        rec_list = gtf_symbols.get(gene_name)
        hgnc_rec = hgnc_symbols.get(gene_name)
        refs = gtf_refs.get(gene_name)
        if rec_list is None and hgnc_rec is None and gtf_refs is None:
            continue
        rec = {"_id": gene_name}
        if rec_list is not None:
            assert len(rec_list) > 0
            rec["gtf"] = rec_list
        if refs is not None:
            rec["gtf-refs"] = sorted(refs)
        if hgnc_rec is not None:
            rec["hgnc"] = hgnc_rec
        records.append(rec)
    print("Result:", len(records), file = sys.stderr)
    return records

#===============================================
meta_rec, symbols, gtf_records = prepareEnsembl(run_args.ensembl)

hgnc_symbols, hgnc_version = loadHGNC(run_args.hgnc)
meta_rec["meta"].append(hgnc_version)

records = crossData(symbols, gtf_records, hgnc_symbols)

print(json.dumps(meta_rec))
for rec in records:
    print(json.dumps(rec))
