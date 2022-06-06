# Utility geneminer

On the current stage the project contain solutions for two connected problems:

**Task 1**: Prepare database for gene symbols nomenclature used in Anfisa project
    
**Task 2**: Complete gene panels (symbol lists) by symbols actual for current Ensembl/GTF nomenclature that possibly can be fresh names or aliases for symbols in the initial list. 


# Sources
We use two sources for work with gene nomenclature:

**Gencode GTF** or **Ensembl/GTF**
        
* Project URL: [https://www.gencodegenes.org/pages/data_format.html](https://www.gencodegenes.org/pages/data_format.html)

* Direct download URL: [ftp://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz](ftp://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz)

**HGNC**

* Download [https://www.genenames.org/download/archive/](https://www.genenames.org/download/archive/)

* We use the full data downloaded from this site in JSON-form: hgnc_complete_set.json

*Note*. For both sources we collect only proteing-coding symbols 

# Mining procedure: determination of "real" symbols

Gene panel definitions are being preparing by specialists for a long time, and gene nomenclature is being essentially changed for this period of time. So symbol names that are used in existing gene panels can not be interpreted in a single proper way in automatical mode. There can be various errors in interpretation of symbols in the initial list: there are possibilities to add improper symbols or to loose proper ones. The purpose of the current code is to except symbols lost. So it is possible that automatical procedure makes another kind of mistakes, and add extra symbols to the final list.

To solve this problem, we do the following:

* For each symbol from GTF/Ensembl nomenclature we evaluate list of **relevant symbols**: current or previous names of the symbol, joined with alias names

* For each known symbols from HGNC we collect **recommended list**: symbols from GTF/Ensembl nomenclature that are relevant to this symbol

# Prerequirements

Setup python virtual environment for the project Anfisa, activate it.

# Task 1: prepare gene database

    > python make_db.py > gene_db.js

By default, the utility tries to find two preloaded files in the current directory, use the following options to change these settings:
    
 * -g HGNC, --hgnc HGNC  *Path to hgnc_complete set*

     Default: ``hgnc_complete_set_2022-04-01.json:``
 
 * -e ENSEMBL, --ensembl ENSEMBL *Path to ensembl archive*

    Default: ``Homo_sapiens.GRCh38.105.chr.gtf``

The resulting file (``gene_db.js``) is the file in JSON format. The first line of it contains descriptor with versions of the source files. Next lines correspond descriptors of symbols, one symbol per line. 
    
Information on gene symbol currently contains the following blocks (if present):

    - Ensembl/GTF information
    
    - HGNC information
    
    - recommended list from Ensembl/GTF nomenclature (see definition above)
    

# Task 2: 
    
On start we have panel (list of symbols). This list can be in CSV format (only ``ensembl_gene_name`` field is used currently). On finish we have panel that consists only from symbols of Ensembl/GTF nomenclature, this list may be more wide than "perfect", however there is "some" guarantee that no symbols form Ensembl/GTF nomenclature is lost in it.
    
    > python up_pannel.py -d gene_db.js {pannel_input_name} > {panel_output_name}
    > python up_pannel.py {pannel_input_name} > {panel_output_name}

Run the first variant only once, with ``gene_db.js`` file that is result of Task 1 procedure: the file ``gene_db.rev`` will be created in the directory and there is no need to use full database file in next runs.

The script uses option:

 * -r REV, --rev REV
 
    Default: ``gene_db.rev``
    
The script stores resulting symbols list in ``stdout`` and reports problems in ``stderr``
    
The project logic assumes that GTF information is ingested to MySQL, so the required data is prepared for geneminer purposes as result of such a call:

