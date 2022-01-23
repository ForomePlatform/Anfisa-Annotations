import json, sys, os
from argparse import ArgumentParser

#========================================
parser = ArgumentParser(
    description = "Data ingestion utiities")
parser.add_argument("--mode", "-m",
    help = "Mode: gtf/pharmgkb/gtex/clinvar")
parser.add_argument("config", nargs = 1)

args = parser.parse_args()

assert os.path.exists(args.config[0]), (
    "Config file not found: " + str(args.config[0]))

with open(args.config[0], "r", encoding = "utf-8") as inp:
    config = json.loads(inp.read())

#========================================
std_db_host    = config["db.host"]
std_db_port    = config["db.port"]
std_user       = config["db.user"]
std_password   = config["db.password"]
std_database   = config.get("database")

#========================================
if args.mode == "pharmgkb":
    from pharmgkb.ca import ingestCA
    from pharmgkb.ca_meta import ingestCAmeta
    from pharmgkb.ca_meta2ca import ingestCAmeta2CA
    from pharmgkb.spa import ingestSPA
    from pharmgkb.vda import ingestVDA
    from pharmgkb.vda2spa import ingestVDA2SPA
    from pharmgkb.vfa import ingestVFA
    from pharmgkb.vfa2spa import ingestVFA2SPA
    from pharmgkb.vpa import ingestVPA
    from pharmgkb.vpa2spa import ingestVPA2SPA
    from pgkb_retab import pgkbReTab

    ingestCA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/clinical_ann.tsv')

    ingestCAmeta(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/clinical_ann_metadata.tsv')

    ingestCAmeta2CA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/clinical_ann_metadata.tsv')

    ingestSPA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/study_parameters.tsv')

    ingestVDA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/var_drug_ann.tsv')

    ingestVDA2SPA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/var_drug_ann.tsv')

    ingestVFA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/var_fa_ann.tsv')

    ingestVFA2SPA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/var_fa_ann.tsv')

    ingestVPA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/var_pheno_ann.tsv')

    ingestVPA2SPA(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"],
        batch_size = config["pharmgkb.batch_size"],
        filename   = config["pharmgkb.path"] + '/var_pheno_ann.tsv')

    pgkbReTab(
        db_host    = config.get("pharmgkb.db.host", std_db_host),
        db_port    = config.get("pharmgkb.db.port", std_db_port),
        user       = config.get("pharmgkb.db.user", std_user),
        password   = config.get("pharmgkb.db.password", std_password),
        database   = config["pharmgkb.database"])

    sys.exit()

#========================================
if args.mode == "gtf":
    from gtf import ingestGTF
    ingestGTF(
        db_host    = config.get("gtf.db.host", std_db_host),
        db_port    = config.get("gtf.db.port", std_db_port),
        user       = config.get("gtf.db.user", std_user),
        password   = config.get("gtf.db.password", std_password),
        database   = config["gtf.database"],
        batch_size = config["gtf.batch_size"],
        filename  = config["gtf.filename"])
    sys.exit()

#========================================
if args.mode == "gtex":
    from gtex import ingestGTEX
    ingestGTEX(
        db_host    = config.get("gtex.db.host", std_db_host),
        db_port    = config.get("gtex.db.port", std_db_port),
        user       = config.get("gtex.db.user", std_user),
        password   = config.get("gtex.db.password", std_password),
        database   = config["gtex.database"],
        batch_size = config["gtex.batch_size"],
        filename  = config["gtex.filename"])
    sys.exit()

#========================================
if args.mode == "clinvar":
    from clinvar import ingestCLINVAR
    ingestCLINVAR(
        db_host    = config.get("clinvar.db.host", std_db_host),
        db_port    = config.get("clinvar.db.port", std_db_port),
        user       = config.get("clinvar.db.user", std_user),
        password   = config.get("clinvar.db.password", std_password),
        database   = config["clinvar.database"],
        batch_size = config["clinvar.batch_size"],
        summary_fname  = config["clinvar.variant_summary_file"],
        xml_fname = config["clinvar.XML_FILE"])
    sys.exit()

assert False, "Unsupported ingest mode: " + str(args.mode)
