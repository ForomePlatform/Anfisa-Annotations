import sys, codecs, logging
from argparse import ArgumentParser

from utils.json_conf import loadJSonConfig
from a_rocksdb.a_storage import AStorage
from a_rocksdb.a_schema import ASchema
from ingest import getIngestModeSetup
#=====================================
try:
    sys.stderr = codecs.getwriter('utf8')(sys.stderr.detach())
    sys.stdout = codecs.getwriter('utf8')(sys.stdout.detach())
except Exception:
    pass

#=====================================
if __name__ == '__main__':
    logging.root.setLevel(logging.INFO)

    parser = ArgumentParser()
    parser.add_argument("-c", "--config",
        help = "AStrorage configuration file",
        default = "astorage.cfg")
    parser.add_argument("-d", "--dbname", help = "db name, "
        "by default equals to mode name")
    parser.add_argument("--update",
        action = "store_true", help = "Update mode")
    parser.add_argument("--dummy",
        action = "store_true", help = "Dummy DB mode")
    parser.add_argument("-m", "--mode",
        help = "Mode")
    args = parser.parse_args()

    db_config = loadJSonConfig(args.config)
    a_storage = AStorage(db_config, dummy_mode = args.dummy)

    db_name = args.dbname
    if not db_name:
        db_name = args.mode
    schema_cfg, reader_func = getIngestModeSetup(args.mode)
    reader_data = reader_func(db_config["create"][args.mode])

    assert schema_cfg is not None
    a_schema = ASchema(a_storage, args.mode, db_name, schema_cfg,
        update_mode = args.update)

    a_storage.activate()

    for key, record in reader_data.read():
        a_schema.putRecord(key, record)

    a_schema.close()
    a_storage.deactivate()
