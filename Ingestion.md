# Ingestion of annotation sources into Annotation databases

<!-- toc -->

- [Preparation:](#preparation)
  * [Setup MySQL](#setup-mysql)
  * [Setup AStorage service](#setup-astorage-service)
  * [Setup Python virtual environment](#setup-python-virtual-environment)
  * [Configuration](#configuration)
- [Reference DNA sequences](#reference-dna-sequences)
  * [GRCh37 (HG19) Assembly](#grch37-hg19-assembly)
  * [GRCh38 (HG38) Assembly](#grch38-hg38-assembly)
- [The Genome Aggregation Database (gnomAD)](#the-genome-aggregation-database-gnomad)
  * [gnomAD Option 1](#gnomad-option-1)
  * [gnomAD Option 2](#gnomad-option-2)
- [Database for nonsynonymous SNPs' functional predictions (dbNSFP)](#database-for-nonsynonymous-snps-functional-predictions-dbnsfp)
- [Database of Single Nucleotide Polymorphisms (dbSNP) Support Center](#database-of-single-nucleotide-polymorphisms-dbsnp-support-center)
- [ClinVar](#clinvar)
- [Genomic Evolutionary Rate Profiling (GERP) Scores](#genomic-evolutionary-rate-profiling-gerp-scores)
- [Gencode GTF](#gencode-gtf)
  * [GTF in MySQL](#gtf-in-mysql)
  * [GTF in AStorage](#gtf-in-astorage)
- [SpliceAI](#spliceai)
  * [SpliceAI Step 1: Splitting by chromosomes](#spliceai-step-1-splitting-by-chromosomes)
  * [SpliceAI Step 2: Deep compilation](#spliceai-step-2--deep-compilation)
  * [SpliceAI Step 3: Push to RocksDB](#spliceai-step-3--push-to-rocksdb)
- [PharmGKB (Pharmacogenomics)](#pharmgkb-pharmacogenomics)
- [GTEx](#gtex)

<!-- tocstop -->

## Preparation:

### Setup MySQL

### Setup AStorage service

Described in [AStorage subproject](a_storage/README.md)

### Setup Python virtual environment

Make sure that the following directories or corresponding packages 
are added to $PYTHONPATH:

* [Ingestion into MySQL](ingestion)
* [Ingestion into RocksDB aka AStorage](a_storage)

### Configuration

* Make a copy the file [config_proto.js](ingestion/config_proto.js) 
  to a file named `config.js`:

      cp config_proto.js config.js

  Edit this local file `config.js`: replace '?'signs with meaningful
  values. The file should be formatted as valid JSON file.

* Make a copy the file [astorage.cfg.template](a_storage/astorage.cfg.template)
  to a file named `astorage.cfg`:

      cp astorage.cfg.template astorage.cfg
  
  Configure this file

## Reference DNA sequences

Reference DNA sequences for assemblies are ingested in
[FASTA](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format)
format and can be downloaded from US National Library of Medicine
(National Center for Biotechnology Information)
[website](https://www.ncbi.nlm.nih.gov/assembly/).



### GRCh37 (HG19) Assembly

* [Project URL](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)
* Last tested version: patch 13

| Source (downloaded) size | Database size | Time to ingest  | Source file URL                                                                                                                                                                  |
|--------------------------|---------------|-----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|1.6 GB | 1.6 GB | 7 minutes | [p13](http://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz) | 

Command:

    python create_db.py -m fasta

Configuration in `astorage.cfg`:

    "fasta_hg19": "${WORK}/prep/fasta/hg19.fasta.gz",

### GRCh38 (HG38) Assembly



* Last tested version: patch 13
* [Project URL](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39)

| Source (downloaded) size | Database size | Time to ingest  | Source file URL                                                                                                                     |
|--------------------------|---------------|-----------------|-------------------------------------------------------------------------------------------------------------------------------------|
|1.6 GB | 1.6 GB | 7 minutes | FTP: `ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13` | 

Configuration in `astorage.cfg`:

    "fasta_hg38": "${WORK}/prep/fasta/hg38.fasta.gz"

## The Genome Aggregation Database (gnomAD)

* Last tested version: 2.1.1
* [URL](https://gnomad.broadinstitute.org/)

| Source (downloaded) size | Database size | Time to ingest | Source file URL                                                                                                                            |
|--------------------------|---------------|----------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| 21 GB                    | 12 GB         | 16 hours       | [Downloads](https://gnomad.broadinstitute.org/downloads) | 

There two options to ingest the data

### gnomAD Option 1

Command:

    python create_db.py -m gnomAD

Configuration parameters:

    "genome_file_list": "${LOAD}/gnomad.2.1/gnomad.genomes.r2.1.1.sites.*.vcf.bgz",
    "exome_file_list": "${LOAD}/gnomad.2.1/gnomad.exomes.r2.1.1.sites.*.vcf.bgz"

### gnomAD Option 2

Command:

    python -m ingest.a_gnomad211 DIR <genome_file_list> <exome_file_list> <output directory>


Configuration in `astorage.cfg`:

    "direct_file_list": "${WORK}/prep/gnomad/gnomad_dir_*.js.gz"


## Database for nonsynonymous SNPs' functional predictions (dbNSFP)

* A One-Stop Database of Functional Predictions and Annotations for Human Non-synonymous and Splice Site SNVs
* [Project URL](https://sites.google.com/site/jpopgen/dbNSFP)
* [Additional Project URL](http://database.liulab.science/dbNSFP)
* Last tested version: 4.0.a

| Source (downloaded) size | Database size | Time to ingest | Source file URL                                                                                                                                                       |
|--------------------------|---------------|----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 25 GB                    | 13 GB         | 12 hours       | FTP: `ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.0a.zip`, [Google Drive](https://drive.google.com/file/d/1BNLEdIc4CjCeOa7V7Z8n8P8RHqUaF5GZ/view?usp=sharing) | 

Command:

    python create_db.py -m dbNSFP

Configuration in `astorage.cfg`:

    "file_list": "${LOAD}/dbNSFP4/dbNSFP4.0a_variant.chr*.gz"


## Database of Single Nucleotide Polymorphisms (dbSNP) Support Center

* [Project URL](https://www.ncbi.nlm.nih.gov/snp/)
* Download URL: [https://ftp.ncbi.nih.gov/snp/](https://ftp.ncbi.nih.gov/snp/)

Command:

    python create_db.py -m dbSNP


## ClinVar

* [Project URL](https://www.ncbi.nlm.nih.gov/clinvar/)
* Two files are required to create local copy
  * CSV File contains data
  * XML File contains data and metadata


| Source (downloaded) size | Database size | Time to ingest | Source file URL                                                                                                                                                                  |
|--------------------------|---------------|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ? GB (two files)         | ? GB          | ? hours        | [CSV](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz ), [XML](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/) | 


Command:

    python ingest.py -m clinvar config.js

Configuration:

* [config.js](ingestion/config_proto.js)

## Genomic Evolutionary Rate Profiling (GERP) Scores 

* [Project URL](http://mendel.stanford.edu/SidowLab/downloads/gerp/)


| Source (downloaded) size | Database size | Time to ingest | Source file URL |
|--------------------------|---------------|----------------|----------------|
| 5.9 GB        | 5.3 GB        | 26 hours       | [Archive](http://mendel.stanford.edu/SidowLab/downloads/gerp/hg19.GERP_scores.tar.gz)      | 

Command:

    python create_db.py -m Gerp

Configuration in `astorage.cfg`:

    "file_list": "${LOAD}/Gerp/chr*.maf.rates"



## Gencode GTF

* [Project URL](https://www.gencodegenes.org/pages/data_format.html)
* This data is loaded to both AStorage and MySQL (for legacy reasons)


| Source (downloaded) size | Database size | Time to ingest | Source file URL                                                                                                                                                                             |
|--------------------------|---------------|----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 45 MB                    | 30 MB         | 7 minutes      | [General downloads](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/); Direct Download: `ftp://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz` | 


### GTF in MySQL

Command:

    python ingest.py -m gtf config.js

Configuration:

* [config.js](ingestion/config_proto.js)
  

### GTF in AStorage

Command:
  
    python create_db.py -m gtf

Configuration in `astorage.cfg`:

    "file_list":"${LOAD}/Homo_sapiens.GRCh38.99.chr.gtf.gz"

> Note: during the loading process a temporary file `_sort.tmp` is 
> created in the working directory


## SpliceAI

* [GitHub URL](https://github.com/Illumina/SpliceAI)
* Last tested version: v1pre3
* Downloads require free registration
* Files:
  * spliceai_scores.masked.snv.hg38.vcf.gz
  * spliceai_scores.masked.indel.hg38.vcf.gz


| Source (downloaded) size | Database size | Time to ingest | Source file URL                                                                                                                                                                      |
|--------------------------|---------------|----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 92 GB                    | 70 GB         | ~48 hours      | [Downloads](https://basespace.illumina.com/analyses/194103939/files/236418325?projectId=66029966) | 


This process is a three steps process.

### SpliceAI Step 1: Splitting by chromosomes

This step can be done in two processes and takes about 12 hours. It saves 
time required to read VCF format.

* Command:

      python -m ingest.a_spliceai
* Inputs:
  * spliceai_scores.masked.indel.hg38.vcf.gz 
  * spliceai_scores.masked.snv.hg38.vcf.gz
* Outputs:
  * Glob: 
    * indel.spliceai.chr*.vcf.gz
    * snv.spliceai.chr*.vcf.gz
    

### SpliceAI Step 2:  Deep compilation

At this step SpliceAI data is being prepared to be ingested into RocksDB. 
The process can be run in parallel for every chromosome. Can take around 
24 hours if run fully in parallel.

Commands:

    #!/bin/bash
    CHROM=$1
    python -u create_db.py -m SpliceAI  \
    -c deep_spliceai.cfg --deepcomp  \
    --metavar CHROM ${CHROM}  \
    -d ${CHROM} &> c_spl.${CHROM}.log \

Configuration:

Requires a special configuration file `deep_spliceai.cfg`.

    {
	"file-path-def": {
    		"TEST": "/home/trifon/work/MD/data_ex",
    		"WORK": "/home/trifon/work/MD/anno/AStorage"
	},
	"db-dir": "${WORK}/prep/spliceai/deepcomp",
	"schema-dir": "${WORK}/prep/spliceai/schema",
	"samples-count": 300,
	"load-keep-schema-sec": 300,
	"create": {
    	   "SpliceAI": {
        		"indel_file_list": "${TEST}/spliceai/prep/indel.spliceai.${CHROM}.vcf.gz",
        		"snv_file_list": "${TEST}/spliceai/prep/snv.spliceai.${CHROM}.vcf.gz"
    	   }
	}
}

### SpliceAI Step 3:  Push to RocksDB
           
Push to RocksDB takes about 1.5 hours

Command:

    python -u create_db.py -m SpliceAI -d SpliceAI \
    --deepstorage /projects/AStorage/prep/spliceAI/deepcomp \
    --deepschema  /projects/AStorage/prep/spliceAI/schema


## PharmGKB (Pharmacogenomics)


* Project URL: [https://www.pharmgkb.org/](https://www.pharmgkb.org/)
* Download URL: [https://www.pharmgkb.org/downloads](https://www.pharmgkb.org/downloads)
* File: Variant Annotations Help File (annotations.zip)

Command:

    python ingest.py -m pharmgkb config.js

Configuration:

* [config.js](ingestion/config_proto.js)
                

## GTEx

* Project URL: [https://www.gtexportal.org/home/](https://www.gtexportal.org/home/)
* Download URL: [https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz)

Command:

    python ingest.py -m gtex config.js

Configuration:

* [config.js](ingestion/config_proto.js)
