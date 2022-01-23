# ingestion subproject

The annotation process provided by Anfisa-Annotation project requires
data from various sources in form of tables in instance of MySQL.

This subproject provides ingestion of these sources into mySQL.
It is written on Python3 and can be used in the following way.

1. Clone or download copy of repository Anfisa-Annotation and
go to ingestion directory

2. Make a copy the file config_proto.js with name config.js:

> cp config_proto.js config.js

Then edit this local file config.js: replace '?'signs with meaningful
values. The file should be formatted in proper JSON form.

3. If the config file is set up properly, one can perform one by one
the modes of ingestion process:

> python3 ingest.py -m _mode_ config.js

Attention: each mode is a long process, so it is recommended to
start it in in a safe way (immune to hangups):

> nohup python3 ingest.py -m _mode_ config.js &> log.txt &

## Modes currently available:

**gtf pharmgkb gtex clinvar**

Mode gtf
--------
Gencode GTF
* Project URL: [https://www.gencodegenes.org/pages/data_format.html](https://www.gencodegenes.org/pages/data_format.html)
* Downloads URL: [https://www.pharmgkb.org/downloads](https://www.pharmgkb.org/downloads)
* Direct download URL: [ftp://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz](ftp://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz)

Mode pharmgkb
--------------
PharmGKB (Pharmacogenomics)

* Project URL: [https://www.pharmgkb.org/](https://www.pharmgkb.org/)
* Download URL: [https://www.pharmgkb.org/downloads](https://www.pharmgkb.org/downloads)
* File: Variant Annotations Help File (annotations.zip)

Mode gtex
---------
GTEx

* Project URL: [https://www.gtexportal.org/home/](https://www.gtexportal.org/home/)
* Download URL: [https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz)

Mode clinvar
-----------
ClinVar

* Project URL: [https://www.ncbi.nlm.nih.gov/clinvar/](https://www.ncbi.nlm.nih.gov/clinvar/)
* CSV File URL (contains data): [https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz ](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz )
* XML File URL (contains data and metadata): [https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/)

**Note**: From newer version of projects Anfisa/Anfisa-Annotations
(v.0.6) the following modes are moved from MySQL to AStorage/RocksDB support, 
the code for these modes is located in ../a_storage/ingest:

**hg19 hg38 gerp gnomad spliceai dbnsfp4**

