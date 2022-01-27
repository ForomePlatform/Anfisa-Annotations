# Utility geneminer

Gene panel definitions are being preparing by specialists for a long time, and gene nomenclature is being essentially chaning for thisperiod of time. 
The purpose of the utillity is to fix gene panel definition with current version of ensembl/GTF gene nomenclature.


# Sources

Gencode GTF
-----------
* Project URL: [https://www.gencodegenes.org/pages/data_format.html](https://www.gencodegenes.org/pages/data_format.html)
* Direct download URL: [ftp://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz](ftp://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz)

The project logic assumes that GTF information is ingested to MySQL, so the required data is prepared for geneminer purposes as result of such a call:

> mysql -u<?> -p<?> -e 'SELECT distinct chromosome, gene from ensembl_hg38.GTF_gene' > gtf_genes.txt

HGNC
----
* Download [https://www.genenames.org/download/archive/](https://www.genenames.org/download/archive/)

We use the full data downloaded from this site in JSON-form: hgnc_complete_set.json

# Build index

    > python3 make_index.py > gtf_rev.txt

Default options for this utility assume that the files **gtf_genes.txt** and **hghgnc_complete_set.json** (see above) are prepared and llocated in the current directory. The result file gtf_ref.txt is index.

# Gene panel format

We use simple text format for gene panels: each symbol should be located in separated line.
There could be comments in the text file, with leading '#' symbol. There can be inline comments or full line comments. Inline comments are ignored in parsing, full line comments are kept in processing.

There are two special form for coment line, with leading two symbols '##', used for mark ID of pannel and references

    ## ID: ...
    
    ## http...

# Processing panel file

    > python3 up_pannel.py {pannel_input_name} > {panel_output_name}

The script can be run from any directory, the index by default is the file **gtf_rev.txt** that should be located in the same directory as up_pannel.py. 

By default the output of the process is the fixed pannel text format. All errors and remarks are printed to <stderr> channel. With option -C (--calm) the output is empty, only log information is printed. 
