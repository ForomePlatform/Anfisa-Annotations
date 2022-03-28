import argparse
import re

import colorlog
import logging

import numpy as np
import pandas as pd


def search_gene_name(
        panel_gene_name,
        ensembl_db,
        hgnc_symbols,
        hgnc_aliases,
        hgnc_previous
):
    """
    1. Collect all possible genes with name panel_gene_name:
    symbol == panel_gene_name
    panel_gene_name in alias names
    panel_gene_name in previous names

    2. Choose from ensembl database all rows with gene_name from
    (symbol, alias, previous_name) of all possible genes

    :param panel_gene_name: gene name as in panel list
    :param ensembl_db: gtf file with ensembl annotation (genes only)
    :param hgnc_symbols: dict symbol -> {"symbol": set(), "previous": set(), "alias": set()}
    :param hgnc_previous: dict previous names prev_name -> [symbols]
    :param hgnc_aliases: dict alias names alias_name -> [symbols]

    :return:
    """

    gene_hgnc_symbols = {
        "alias_symbol":
            hgnc_aliases[panel_gene_name] if panel_gene_name in hgnc_aliases else set(),  # panel_gene_name as an alias_symbol
        "prev_symbol":
            hgnc_previous[panel_gene_name] if panel_gene_name in hgnc_previous else set()  # panel_gene_name as a prev_symbol
    }
    gene_hgnc_all_names = {
        "symbol": [hgnc_symbols[panel_gene_name], ] if panel_gene_name in hgnc_symbols else [],
        "alias_symbol": [hgnc_symbols[alias_gene_name] for alias_gene_name in gene_hgnc_symbols["alias_symbol"]],
        "prev_symbol": [hgnc_symbols[prev_gene_name] for prev_gene_name in gene_hgnc_symbols["prev_symbol"]]}

    ensembl_db["relationship"] = np.NaN

    for gene_symbol_as_type, gene_symbol_as_type_symbol_list in gene_hgnc_all_names.items():
        for gene_symbol_as_type_symbol in gene_symbol_as_type_symbol_list:
            for gene_symbol_as_type_symbol_type, gene_symbol_as_type_symbol_symbols in gene_symbol_as_type_symbol.items():
                rel = \
                    f"the gene_name is a(n) {gene_symbol_as_type_symbol_type} of {gene_symbol_as_type} of '{panel_gene_name}'"
                gene_symbol_as_type_symbol_symbols = gene_symbol_as_type_symbol_symbols.difference(panel_gene_name)
                ensembl_db.loc[ensembl_db["gene_name"].isin(gene_symbol_as_type_symbol_symbols), "relationship"] = rel

    ensembl_db.loc[ensembl_db["gene_name"] == panel_gene_name, "relationship"] = "by gene_name"

    gene_name_ensembl = ensembl_db.loc[~ensembl_db["relationship"].isna()]
    return gene_name_ensembl


def gene_description(gene):
    """
    get a textual description of gene found in ensembl by gene_name
    :param gene:
    :return:
    """
    return f"\tcoords:\tchr{gene.iloc[0]}:{gene.iloc[3]}-{gene.iloc[4]}\n" \
        f"\tdetails:\t{gene.iloc[8]}\n" \
        f"\tmatched:\t{gene.loc['relationship']}"


def get_confirmation(gene_name, gene_name_result, ensembl, confirm):
    """
    Just ask user is this gene was mapped in a right way (if --confirm option was provided)
    :param ensembl:
    :param confirm:
    :param gene_name:
    :param gene_name_result:
    :return:
    """
    if confirm:
        while 1:
            answer = input(
                f"Gene name '{gene_name}' was matched to the gene from GENCODE annotation:\n"
                f"{gene_description(gene_name_result.iloc[0])}\n"
                f"Is it OK?\n"
                f"If yes -- press Enter, otherwise enter valid GENCODE gene id for gene '{gene_name}':\n")
            if answer:
                gene = ensembl.loc[ensembl["gene_id"] == answer, :]
                if len(gene):
                    gene_name_result = gene
                else:
                    print(f"No gene was found by id '{answer}', please, try again")
            else:
                break
    else:
        print(
            f"Gene name '{gene_name}' was matched to the gene from GENCODE annotation:\n"
            f"{gene_description(gene_name_result.iloc[0])}\n")

    return gene_name_result.iloc[0]["gene_id"]


def choose_from_or_provide_id(gene_name, gene_name_result, ensembl, confirm):
    """
    Ask user which gene is to be used from list of possibilities
    :param ensembl:
    :param confirm:
    :param gene_name:
    :param gene_name_result:
    :return:
    """
    while 1:
        genes = '\n'.join([
            f"{i}. {gene_description(gene_name_result.iloc[i])}"
            for i in range(len(gene_name_result))
        ])
        answer = input(
            f"Gene name '{gene_name}' was ambiguously matched to the next genes from GENCODE annotation:\n"
            f"{genes}\n\n"
            f"Please choose one option from the list above by entering its number\n"
            f"\tor enter valid GENCODE gene id for gene '{gene_name}'\n"
            f"\tor press Enter to skip this gene:\n")

        if not answer:
            return None

        if re.match(r"\d+", answer):
            choice = int(answer)
            if 0 <= choice < len(gene_name_result):
                gene = ensembl.iloc[[choice]]
                return get_confirmation(gene_name, gene, ensembl, confirm)
            else:
                print(f"Invalid number '{choice}' was entered, please, try again")
        else:
            gene_name_result2 = ensembl.loc[ensembl["gene_id"] == answer, :]
            if len(gene_name_result2):
                return get_confirmation(gene_name, gene_name_result2, ensembl, confirm)
            else:
                print(f"No gene was found by id '{answer}', please, try again")


def provide_id(gene_name, ensembl, confirm):
    """
    Ask user to provide ensembl_gene_id for missing gene names
    :param confirm:
    :param gene_name:
    :param ensembl:
    :return:
    """
    while 1:
        answer = input(
            f"Gene name '{gene_name}' was not matched to any gene from GENCODE annotation.\n"
            f"Please enter valid ensembl_gene_id for gene '{gene_name}' or press Enter to skip this gene:\n")
        if answer:
            gene_name_result = ensembl.loc[ensembl["gene_id"] == answer, :]
            if len(gene_name_result):
                return get_confirmation(gene_name, gene_name_result, ensembl, confirm)
            else:
                print(f"No gene was found by id '{answer}', please, try again")
        else:
            return None


def get_ensembl(gtf_fp):
    print(f"... loading GENCODE annotation from file '{gtf_fp}'")
    gtf_df = pd.read_csv(gtf_fp, sep="\t", comment="#", header=None, dtype={0: 'str'})
    gtf_genes_only = gtf_df.loc[gtf_df.loc[:, 2] == "gene", :].copy()

    # retrieve gene_name
    gtf_genes_only["gene_name"] = gtf_genes_only.iloc[:, 8].str.extract(r'gene_name "([^"]+)"')
    gtf_genes_only["gene_id"] = gtf_genes_only.iloc[:, 8].str.extract(r'gene_id "([^"]+)"')
    gtf_genes_only["gene_version"] = gtf_genes_only.iloc[:, 8].str.extract(r'gene_version "([^"]+)"')

    return gtf_genes_only


def get_hgnc(hgnc_fp):
    print(f"... loading HGNC from file '{hgnc_fp}'")
    dtype = {32: "str", 34: "str", 38: "str", 40: "str", 50: "str"}
    hgnc_df = pd.read_csv(hgnc_fp, sep="\t", dtype=dtype)

    # check if symbol in HGNC is an unique field
    assert len(hgnc_df["symbol"]) == len(hgnc_df["symbol"].unique())
    hgnc_df = hgnc_df.set_index("symbol", drop=False)

    return hgnc_df


def get_hgnc_symbols(hgnc_df):
    hgnc_symbols = {}

    for index, row in hgnc_df.iterrows():
        hgnc_symbols[index] = {}

        aliases = row["alias_symbol"]
        if not pd.isna(aliases):
            hgnc_symbols[index]["alias_symbol"] = set(aliases.split("|"))
        else:
            hgnc_symbols[index]["alias_symbol"] = set()

        prev_symbols = row["prev_symbol"]
        if not pd.isna(prev_symbols):
            hgnc_symbols[index]["prev_symbol"] = set(prev_symbols.split("|"))
        else:
            hgnc_symbols[index]["prev_symbol"] = set()

        hgnc_symbols[index]["symbol"] = {index}

    return hgnc_symbols


def get_hgnc_symbols_map(hgnc_symbols):
    # dictionary of previous names: prev_name -> [symbols]
    # dictionary of alias names: alias_name -> [symbols]
    hgnc_prev2s = {}
    hgnc_alias2s = {}

    for symbol, other_symbols in hgnc_symbols.items():
        for prev_symbol in other_symbols["prev_symbol"]:
            if prev_symbol in hgnc_prev2s and symbol not in hgnc_prev2s[prev_symbol]:
                hgnc_prev2s[prev_symbol].add(symbol)
            else:
                hgnc_prev2s[prev_symbol] = {symbol}

        for alias_symbol in other_symbols["alias_symbol"]:
            if alias_symbol in hgnc_alias2s and symbol not in hgnc_alias2s[alias_symbol]:
                hgnc_alias2s[alias_symbol].add(symbol)
            else:
                hgnc_alias2s[alias_symbol] = {symbol}
    return hgnc_alias2s, hgnc_prev2s


def negotiate_panel_gene_names(
        gene_panel_fp, ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, confirm):

    print(f"... reading panel gene names from file '{gene_panel_fp}'\n")
    name2id = {}
    id2name = {}
    res_lines = []

    with open(gene_panel_fp, "r") as gene_panel:
        while 1:
            line = gene_panel.readline()
            if line:
                if line.startswith("#"):
                    res_lines.append(line)
                    continue

                gene_name = line.partition('#')[0].strip()
                if gene_name:

                    if gene_name in name2id:
                        logger.error(f"Gene name '{gene_name}' is duplicated in panel gene list")
                        continue

                    gene_name_result = search_gene_name(
                        gene_name,
                        ensembl,
                        hgnc_symbols,
                        hgnc_alias2s,
                        hgnc_prev2s
                    )

                    if len(gene_name_result) == 1:
                        gene_name_id = get_confirmation(
                            gene_name, gene_name_result, ensembl, confirm)
                    elif len(gene_name_result) > 1:
                        gene_name_id = choose_from_or_provide_id(gene_name, gene_name_result, ensembl, confirm)
                    else:
                        gene_name_id = provide_id(gene_name, ensembl, confirm)

                    if gene_name_id:
                        if gene_name_id in id2name:
                            logger.error(
                                f"Ensemble gene id '{gene_name_id}' was already matched "
                                f"to gene name {id2name[gene_name_id]} from this panel,"
                                f"skipping '{gene_name}' from result")
                        else:
                            name2id[gene_name] = gene_name_id
                            id2name[gene_name_id] = gene_name

                            res_lines.append(f"{gene_name}\t{gene_name_id}\n")

                    print("\n---------------------------------------------------------------\n")

                else:
                    res_lines.append(line)
            else:
                break

    return res_lines


def main(gene_panel_fp, gtf_fp, hgnc_fp, confirm, res_fp):

    ensembl = get_ensembl(gtf_fp)
    hgnc = get_hgnc(hgnc_fp)

    # dictionary from HGNC db: symbol -> {"symbol": set(), "previous": set(), "alias": set()}
    hgnc_symbols = get_hgnc_symbols(hgnc)
    hgnc_alias2s, hgnc_prev2s = get_hgnc_symbols_map(hgnc_symbols)

    res_lines = negotiate_panel_gene_names(
        gene_panel_fp, ensembl, hgnc_symbols, hgnc_alias2s, hgnc_prev2s, confirm)

    if res_fp is not None:
        with open(res_fp, "w") as res:
            for line in res_lines:
                res.write(line)


def parse_options():
    parser = argparse.ArgumentParser(
        description='Process list of gene names with unambiguity problems solving by user')

    parser.add_argument(
        'gene_panel_fp',
        help="Path to text file with gene names (one gene name per line)",
        type=str
    )

    parser.add_argument(
        '-res-fp',
        '--result-gene-panel-fp',
        help="Path to file with resulted ensembl gene ids, will be written if argument is provided",
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

    main(
        opts.gene_panel_fp,
        opts.ensembl_annotation_fp,
        opts.hgnc_tsv_fp,
        opts.confirm,
        opts.result_gene_panel_fp
    )
