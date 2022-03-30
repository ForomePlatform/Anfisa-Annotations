import argparse
import json
import logging
import os

import colorlog
import pandas as pd
import requests
from Bio import Entrez

from pipeline.projects.geneminer.gene_name_negotiation import get_ensembl, get_hgnc, get_hgnc_symbols, \
    get_hgnc_symbols_map, search_and_negotiate, get_confirmation


# ontologyId = "HP:0003198"

# gene_response = requests.get(f"https://hpo.jax.org/api/hpo/gene/{genes[0]['entrezGeneId']}", verify=False)
# gene = json.loads(gene_response.text)["gene"]


def get_start_stop(gene_info):
    start = int(gene_info['ChrStart'])
    stop = int(gene_info['ChrStop'])

    if start > stop:
        m = start
        start = stop
        stop = m
    start += 1
    stop += 1
    return start, stop


def hpo_csv2ensembl(
        hpo_id, hpo_csv_fp,
        ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, hpo_name2id,
        confirm, do_not_confirm_coords):
    hpo_genes_df = pd.read_csv(hpo_csv_fp, header=0, dtype={0: int, 1: str, 2: str}, sep="\t")
    gene_entrez_ids = hpo_genes_df["GENE_ENTREZ_ID"]
    gene_entrez_names = hpo_genes_df["GENE_SYMBOL"]
    gene_entrez_names_set = set(gene_entrez_names)

    assert len(gene_entrez_names) == len(gene_entrez_names_set)

    return hpo2ensembl(
        hpo_id, gene_entrez_ids, gene_entrez_names,
        ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, hpo_name2id,
        confirm, do_not_confirm_coords)


def hpo_id2ensembl(
        hpo_id, ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, hpo_name2id,
        confirm, do_not_confirm_coords):

    # HPO gene list download
    hpo_genes_response = requests.get(f"https://hpo.jax.org/api/hpo/term/{hpo_id}/genes?max=-1", verify=False)
    hpo_genes = json.loads(hpo_genes_response.text)["genes"]
    gene_entrez_ids = [gene['entrezGeneId'] for gene in hpo_genes]
    gene_entrez_names = [gene['entrezGeneSymbol'] for gene in hpo_genes]
    gene_entrez_names_set = set(gene_entrez_names)

    assert len(gene_entrez_names) == len(gene_entrez_names_set)

    return hpo2ensembl(
        hpo_id, gene_entrez_ids, gene_entrez_names, ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, hpo_name2id, confirm, do_not_confirm_coords)


def hpo2ensembl(
        hpo_id, gene_entrez_ids, gene_entrez_names,
        ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, hpo_name2id,
        confirm, do_not_confirm_coords):

    Entrez.email = "some.email@gmail.com"
    handle = Entrez.efetch(db="gene", id=gene_entrez_ids, rettype="gb", retmode="xml")
    genes_entrez = Entrez.read(handle)

    assert len(gene_entrez_ids) == len(genes_entrez)

    gene_entrez_name_2_ensembl_id = {
        gene["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]:
        gene_ref_db["Dbtag_tag"]["Object-id"]["Object-id_str"]
        for gene in genes_entrez
        for gene_ref_db in gene["Entrezgene_gene"]["Gene-ref"]["Gene-ref_db"]
        if gene_ref_db["Dbtag_db"] == "Ensembl"}

    gene_ensembl_id_2_entrez_name = {}

    # with open(output_fp, "w") as out:

    for gene_entrez_name, gene_ensembl_id in gene_entrez_name_2_ensembl_id.items():
        if gene_ensembl_id in gene_ensembl_id_2_entrez_name:
            raise Exception(
                f"gene_ensembl_id '{gene_ensembl_id}' is duplicated: "
                f"'{gene_entrez_name}' and '{gene_ensembl_id_2_entrez_name[gene_ensembl_id]}'")
        else:
            gene_ensembl_id_2_entrez_name[gene_ensembl_id] = gene_entrez_name

    for gene_entrez_name, gene_entrez_id in zip(gene_entrez_names, gene_entrez_ids):
        if gene_entrez_name in hpo_name2id:
            print(f"Gene name '{gene_entrez_name}' was already considered")
            continue
        gene_handle = Entrez.esummary(db="gene", id=str(gene_entrez_id))
        gene_info = Entrez.read(gene_handle)['DocumentSummarySet']['DocumentSummary'][0]
        print(
            f"\tgene entrez name:\t{gene_entrez_name}\n"
            f"\tgene entrez id:\t{gene_entrez_id}\n"
            f"\tgene entrez name:\t{gene_info['NomenclatureName']}\n"
            f"\tgene entrez synonyms:\t{gene_info['OtherAliases']}\n"
            f"\tgene entrez location:\n\t{gene_info['GenomicInfo'][0]}\n"
        )
        start, stop = get_start_stop(gene_info['GenomicInfo'][0])

        coords = {
            "chr": gene_info['GenomicInfo'][0]['ChrLoc'],
            "start": start,
            "stop": stop
        }

        if gene_entrez_name not in gene_entrez_name_2_ensembl_id:
            # print additional info from entrez
            print(
                f"ensembl id was not specified for gene '{gene_entrez_name}' '{gene_entrez_id}'")

            gene_result = search_and_negotiate(
                gene_entrez_name,
                ensembl,
                hgnc_symbols,
                hgnc_alias2s,
                hgnc_prev2s,
                confirm=True,
                coords=coords,
                do_not_confirm_coords=do_not_confirm_coords
            )
            if gene_result:
                gene_name_id = gene_result["gene_id"]
                ensembl_gene_name = gene_result["gene_name"]
                if gene_name_id in gene_ensembl_id_2_entrez_name:
                    logger.error(
                        f"Ensemble gene id '{gene_name_id}' was already matched "
                        f"to gene name {gene_ensembl_id_2_entrez_name[gene_name_id]} from this panel,"
                        f"skipping '{gene_entrez_name}' from result")
                    continue
            else:
                gene_name_id = None
                ensembl_gene_name = None
        else:
            gene_name_ensembl = ensembl.loc[
                ensembl["gene_id"] == gene_entrez_name_2_ensembl_id[gene_entrez_name]].copy()
            if len(gene_name_ensembl):
                gene_name_ensembl["relationship"] = "entrez_id -> ensembl_id"
                gene_name_ensembl["coords_coincidence"] = False
                gene_name_ensembl.loc[
                                (gene_name_ensembl.iloc[:, 0] == coords["chr"]) &
                                (gene_name_ensembl.iloc[:, 3] == coords["start"]) &
                                (gene_name_ensembl.iloc[:, 4] == coords["stop"]), "coords_coincidence"] = True
                gene_result = get_confirmation(
                    gene_entrez_name, gene_name_ensembl, ensembl, confirm, do_not_confirm_coords)
                gene_name_id = gene_result["gene_id"]
                ensembl_gene_name = gene_result["gene_name"]
            else:
                raise Exception(
                    f"Nothing was found in ensembl db by id '{gene_entrez_name_2_ensembl_id[gene_entrez_name]}'")

        gene_entrez_name_2_ensembl_id[gene_entrez_name] = gene_name_id
        gene_ensembl_id_2_entrez_name[gene_name_id] = gene_entrez_name

        # out.write(f"{gene_entrez_name}\t{gene_name_id}\n")
        # if gene_entrez_name in hpo_name2id:
        #     hpo_name2id[gene_entrez_name]["source"] += f";{hpo_id}"
        # else:
        hpo_name2id[gene_entrez_name] = {
            "gene_name": gene_entrez_name,
            "ensembl_gene_id": gene_name_id,
            "source": hpo_id,
            "ensembl_gene_name": ensembl_gene_name,
            "additional_info": ""
        }

        print("\n---------------------------------------------------------------\n")

    return hpo_name2id


def main(ontology_id, gtf_fp, hgnc_fp, output_fp, confirm, do_not_confirm_coords):
    hpo_name2id = {}

    ensembl = get_ensembl(gtf_fp)
    hgnc = get_hgnc(hgnc_fp)

    # dictionary from HGNC db: symbol -> {"symbol": set(), "previous": set(), "alias": set()}
    hgnc_symbols = get_hgnc_symbols(hgnc)
    hgnc_alias2s, hgnc_prev2s = get_hgnc_symbols_map(hgnc_symbols)

    hpo_name2id = hpo_id2ensembl(
        ontology_id,
        ensembl, hgnc_alias2s, hgnc_prev2s, hgnc_symbols, hpo_name2id,
        confirm=confirm, do_not_confirm_coords=do_not_confirm_coords)

    hpo_name2id_df = pd.DataFrame(hpo_name2id)
    hpo_name2id_df.to_csv(output_fp, index=False)


def parse_options():
    parser = argparse.ArgumentParser(
        description='Process panel from Human Phenotype Ontology')

    parser.add_argument(
        'ontology_id',
        help="ID from HPO",
        type=str
    )

    parser.add_argument(
        'output',
        help="Path to file with resulted ensembl gene ids",
        type=str,
        default=None
    )

    parser.add_argument(
        '-ens',
        '--ensembl-annotation-fp',
        help="Path to tsv file with ensembl annotation",
        type=str,
        default="Homo_sapiens.GRCh38.105.chr.gtf"
    )

    parser.add_argument(
        '-hgnc',
        '--hgnc-tsv-fp',
        help="Path to tsv file with HGNC database",
        type=str,
        default="hgnc_complete_set.txt"
    )

    parser.add_argument(
        '--confirm',
        help="If to confirm genes which were unambiguously matched",
        action="store_true"
    )

    parser.add_argument(
        '--do-not-confirm-coords',
        help="If not to confirm genes which coordinates coincided absolutely",
        action="store_true"
    )

    options = parser.parse_args()

    return options


if __name__ == "__main__":
    opts = parse_options()

    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        '%(log_color)s%(levelname)s:%(message)s'))

    logging.basicConfig(
        level=logging.INFO,
        handlers=[handler])

    logger = logging.getLogger()

    # ontologyId = "HP:0003198"
    # gtf_fp = r"C:\Users\kseniya.petrova\projs\Anfisa\Homo_sapiens.GRCh38.105.chr.gtf"
    # hgnc_fp = r"C:\Users\kseniya.petrova\projs\Anfisa\hgnc_complete_set.txt"
    # output_fp = r"C:\Users\kseniya.petrova\projs\Anfisa\outputs\Myopathy.HP_0003198.txt"
    main(
        opts.ontology_id,
        opts.ensembl_annotation_fp,
        opts.hgnc_tsv_fp,
        opts.output,
        opts.confirm,
        opts.do_not_confirm_coords)


def do_list():
    panel_list = {
        "Myopia": "HP:0000545",
        "Talipes equinovalgus": "HP:0001772"
    }
    outputs_dir = r"C:\Users\kseniya.petrova\projs\Anfisa\outputs"

    for name, hp in panel_list.items():
        print(f"{name}: {hp}")
        main(
            hp,
            r"C:\Users\kseniya.petrova\projs\Anfisa\Homo_sapiens.GRCh38.105.chr.gtf",
            r"C:\Users\kseniya.petrova\projs\Anfisa\hgnc_complete_set.txt",
            os.path.join(outputs_dir, f"{name.replace(' ', '_')}.{hp.replace(':', '_')}.txt"),
            confirm=False,
            do_not_confirm_coords=True)
