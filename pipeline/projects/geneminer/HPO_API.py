import argparse
import json
import logging
import os
import sys
import time

import colorlog
import pandas as pd
import requests
from Bio import Entrez
from biomart import BiomartServer

from pipeline.projects.geneminer.gene_name_negotiation import get_ensembl, get_hgnc, get_hgnc_symbols, \
    get_hgnc_symbols_map, search_and_negotiate, get_confirmation, choose_from_or_provide_id, provide_id

logger = logging.getLogger()
# ontologyId = "HP:0003198"

# gene_response = requests.get(f"https://hpo.jax.org/api/hpo/gene/{genes[0]['entrezGeneId']}", verify=False)
# gene = json.loads(gene_response.text)["gene"]


class ServerCached(object):
    _server = None
    _ensembl = None
    urls = ['http://uswest.ensembl.org/biomart', 'http://useast.ensembl.org/biomart']

    @property
    def server(self):
        if self._server is None:
            for url in self.urls:
                retries = 5
                while retries:
                    try:
                        self._server = BiomartServer(url)
                        logger.debug("Connection to biomart server succeeded")
                        return self._server
                    except Exception:
                        ex_type, ex_value, ex_trace = sys.exc_info()
                        retries -= 1
                        if not retries:
                            raise ex_type(ex_value).with_traceback(ex_trace)
                        logger.debug(f"Connection to biomart server {retries} retries left")
                        time.sleep(5)
                        continue
            raise Exception(f"Couldn't connect to BioMart server by any of urls from list: {self.urls}")
        return self._server

    @property
    def ensembl(self):
        if self._ensembl is None:
            self._ensembl = self.server.datasets['hsapiens_gene_ensembl']
        return self._ensembl


server_cached = ServerCached()


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

    gene_entrez_id_2_ensembl_id = {}
    gene_ensembl_id_2_entrez_name = {}

    gene_num = len(gene_entrez_ids)
    batch_size = 32
    batch_num = gene_num // batch_size + int(gene_num % batch_size != 0)
    estimated_time_str = "unknown"
    times = []
    first = True

    for i in range(batch_num):
        print(f'Loading batches of genes from BioMart server: {i + 1}/{batch_num}; estimated time left: {estimated_time_str}',
              end='\r', flush=True)
        start = time.time()

        response = server_cached.ensembl.search({
            'filters': {
                "entrezgene_id": gene_entrez_ids[i*batch_size:(i+1)*batch_size],
                # 'ensembl_transcript_id_version': 'ENST00000318602.12'
            },
            'attributes': ["entrezgene_id", "ensembl_gene_id"]
        }, header=1)

        # response format is TSV
        header = True
        for line in response.iter_lines():
            line = line.decode('utf-8')
            if header:
                header = False
            else:
                ids = line.split("\t")
                if len(ids) == 2:
                    entrez_id = int(ids[0])
                    ensembl_id = ids[1]

                    if ensembl_id in gene_ensembl_id_2_entrez_name:
                        gene_ensembl_id_2_entrez_name[ensembl_id].append(entrez_id)
                        # print(f"ensembl_id '{ensembl_id}' duplicated: {gene_ensembl_id_2_entrez_name[ensembl_id]}")
                    else:
                        gene_ensembl_id_2_entrez_name[ensembl_id] = [entrez_id,]

                    if entrez_id in gene_entrez_id_2_ensembl_id:
                        gene_entrez_id_2_ensembl_id[entrez_id].append(ensembl_id)
                        # print(f"entrez_id '{entrez_id}' duplicated: {gene_entrez_id_2_ensembl_id[entrez_id]}")
                    else:
                        gene_entrez_id_2_ensembl_id[entrez_id] = [ensembl_id, ]
                else:
                    print(f"not 2 ids returned: {line}")

        end = time.time()
        if not first:
            times.append(end - start)
            estimated_time = sum(times) * (batch_num - (i + 1)) / (len(times) * 60)
            estimated_time_str = f"{estimated_time: .2f} min"
        first = False

    Entrez.email = "some.email@gmail.com"

    for gene_entrez_name, gene_entrez_id in zip(gene_entrez_names, gene_entrez_ids):
        if gene_entrez_name in hpo_name2id:
            print(f"Gene name '{gene_entrez_name}' was already considered")
            continue
        gene_handle = Entrez.esummary(db="gene", id=str(gene_entrez_id))
        gene_info = Entrez.read(gene_handle)['DocumentSummarySet']['DocumentSummary'][0]
        genomic_info = gene_info['GenomicInfo'][0] if len(gene_info['GenomicInfo']) else ""
        print(
            f"\tgene entrez name:\t{gene_entrez_name}\n"
            f"\tgene entrez id:\t{gene_entrez_id}\n"
            f"\tgene entrez name:\t{gene_info['NomenclatureName']}\n"
            f"\tgene entrez synonyms:\t{gene_info['OtherAliases']}\n"
            f"\tgene entrez location:\n\t{genomic_info}\n"
        )
        if genomic_info:
            start, stop = get_start_stop(gene_info['GenomicInfo'][0])

            coords = {
                "chr": gene_info['GenomicInfo'][0]['ChrLoc'],
                "start": start,
                "stop": stop
            }
        else:
            coords = None

        if gene_entrez_id not in gene_entrez_id_2_ensembl_id:
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
            if gene_result is not None:
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
            ensembl_gene = ensembl.loc[
                ensembl["gene_id"].isin(gene_entrez_id_2_ensembl_id[gene_entrez_id]), :].copy()
            if len(ensembl_gene):
                ensembl_gene["relationship"] = "entrez_id -> ensembl_id"
                ensembl_gene["coords_coincidence"] = False
                if coords is not None:
                    ensembl_gene.loc[
                                    (ensembl_gene.iloc[:, 0] == coords["chr"]) &
                                    (ensembl_gene.iloc[:, 3] == coords["start"]) &
                                    (ensembl_gene.iloc[:, 4] == coords["stop"]), "coords_coincidence"] = True
                if len(ensembl_gene) == 1:
                    gene_result = get_confirmation(
                        gene_entrez_name, ensembl_gene, ensembl, confirm, do_not_confirm_coords)
                elif len(ensembl_gene) > 1:
                    gene_result = choose_from_or_provide_id(
                        gene_entrez_name, ensembl_gene, ensembl, confirm, do_not_confirm_coords)
                else:
                    gene_result = provide_id(gene_entrez_name, ensembl, confirm)

                gene_name_id = gene_result["gene_id"]
                ensembl_gene_name = gene_result["gene_name"]
            else:
                raise Exception(
                    f"Nothing was found in ensembl db by id '{gene_entrez_id_2_ensembl_id[gene_entrez_id]}'")

        if gene_name_id is not None:
            gene_entrez_id_2_ensembl_id[gene_entrez_id] = gene_name_id
            gene_ensembl_id_2_entrez_name[gene_name_id] = gene_entrez_name

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
