import argparse
import logging
import os

import colorlog
import pandas as pd

from pipeline.projects.geneminer.HPO_API import hpo_id2ensembl, hpo_csv2ensembl
from pipeline.projects.geneminer.gene_name_negotiation import get_ensembl, get_hgnc, get_hgnc_symbols, \
    get_hgnc_symbols_map, search_and_negotiate


def checkpoint(panel_res, output_fp):
    panel_df = pd.DataFrame.from_dict(panel_res, orient="index")
    panel_df["gene_name"] = panel_df.index
    panel_df = panel_df[['gene_name', 'ensembl_gene_name', 'ensembl_gene_id', 'source', 'additional_info']]
    panel_df.to_csv(output_fp, index=False)


def main(
        panel_fp, ontology_ids, ontology_fps, gtf_fp, hgnc_fp, output_fp,
        confirm, do_not_confirm_coords, reanalyse):

    if os.path.exists(output_fp) and not reanalyse:
        panel_df = pd.read_csv(output_fp)
        panel_df.set_index("gene_name", inplace=True)
        panel_res = panel_df.to_dict(orient="index")
    else:
        panel_res = {}

    ensembl = get_ensembl(gtf_fp)
    hgnc = get_hgnc(hgnc_fp)

    # dictionary from HGNC db: symbol -> {"symbol": set(), "previous": set(), "alias": set()}
    hgnc_symbols = get_hgnc_symbols(hgnc)
    hgnc_alias2s, hgnc_prev2s = get_hgnc_symbols_map(hgnc_symbols)

    if ontology_ids:
        for ontology_id in ontology_ids:
            panel_res = hpo_id2ensembl(
                ontology_id,
                ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, panel_res,
                confirm=confirm, do_not_confirm_coords=do_not_confirm_coords)
            checkpoint(panel_res, output_fp)

    if ontology_fps:
        for ontology_fp in ontology_fps:
            ontology_id = os.path.basename(ontology_fp).split(".")[-2].replace("_", ":")
            panel_res = hpo_csv2ensembl(
                ontology_id, ontology_fp,
                ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, panel_res,
                confirm=confirm, do_not_confirm_coords=do_not_confirm_coords)
            checkpoint(panel_res, output_fp)

    if panel_fp is not None:
        panel = pd.read_csv(panel_fp, header=0, comment="#", skip_blank_lines=True, dtype={0: str, 1: str, 2: str}, keep_default_na=False)
        # panel.columns = ["gene_name", ]
        for i, gene in panel.iterrows():
            gene_name = gene["gene_name"]
            gene_source = gene["source"] if "source" in gene else ""
            gene_additional_info = gene["additional_info"] if "additional_info" in gene else ""
            print(f"{gene_name}\n\tsource: {gene_source}\n\tadditional_info: {gene_additional_info}")

            if gene_name in panel_res:
                print(f"Gene name '{gene_name}' was already considered")
                if gene_additional_info:
                    panel_res[gene_name]["additional_info"] += ";" + gene["additional_info"]

                if gene_source:
                    panel_res[gene_name]["source"] += ";" + gene["source"]
            else:
                gene_result = search_and_negotiate(
                    gene_name,
                    ensembl,
                    hgnc_symbols,
                    hgnc_alias2s,
                    hgnc_prev2s,
                    confirm
                )
                if gene_result is not None:
                    ensembl_gene_id = gene_result["gene_id"]
                    ensembl_gene_name = ensembl.loc[ensembl["gene_id"] == ensembl_gene_id, "gene_name"].item()

                    panel_res[gene_name] = {
                        "gene_name": gene_name,
                        "ensembl_gene_id": ensembl_gene_id,
                        "source": gene_source,
                        "ensembl_gene_name": ensembl_gene_name,
                        "additional_info": gene_additional_info
                    }

            checkpoint(panel_res, output_fp)

            print("\n---------------------------------------------------------------\n")


def parse_options():
    parser = argparse.ArgumentParser(
        description='Process custom panel from file with hints from Human Phenotype Ontology')

    parser.add_argument(
        '-p', '--panel-fp',
        help="gene name list filepath",
        type=str
    )

    parser.add_argument(
        '-HPO-ids', '--ontology-ids',
        nargs="*",
        help="IDs from HPO",
        type=str
    )

    parser.add_argument(
        '-HPO-fps', '--ontology-file-paths',
        nargs="*",
        help="file-paths of HPO gene names lists",
        type=str
    )

    parser.add_argument(
        '-out', '--output',
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

    parser.add_argument(
        '--reanalyse',
        help="If to reanalyse all panel genes if output already exists",
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
        opts.panel_fp,
        opts.ontology_ids,
        opts.ontology_file_paths,
        opts.ensembl_annotation_fp,
        opts.hgnc_tsv_fp,
        opts.output,
        opts.confirm,
        opts.do_not_confirm_coords,
        opts.reanalyse
    )

