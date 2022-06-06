# AStorage: technical description

AStorage is a RocksDB support module designed for the Anfisa-Annotations project.

Currently, the module is implemented in Python, but due to technical problems (there is no reliable support for accessing RockDB via Python), the issue of re-implementing the module in Java (support of which is included in the main RockDB project) is being considered in the close future.

The basic principle of use: we use AStorage / RocksDB when we don't want to access via MySQL in following conditions:
- (main) to process large amounts of genome and exome information
- to support specific search features, when we don't want to use MySQL (GTF, where you need to search for “windows of positions”, and not the positions themselves)
- (in the future): we can switch to the same interface the rest of the arrays required for the annotation, but it's not very clear why we should do this: support via RocksDB is relatively expensive.


## Table of contents
* [Service requests](#service-requests)
* [Bulk requests](#bulk-requests)
* [Service structure](#service-structure)
  - [Directories and files](#directories-and-files)
  - [Programs](#programs)
      * [Launches](#launches)
* [Service setup](#service-setup)
  - [Administrative aspects: number of open files](#administrative-aspects-number-of-open-files)
  - [General setup astorage.cfg](#general-setup-astoragecfg)
  - [Schema setup](#schema-setup)
  - [Principles of stacking JSON structures](#principles-of-stacking-json-structures)
  - [Keys and record locking inside RocksDB](#keys-and-record-locking-inside-rocksdb)
* [Supplements](#supplements)
  - [A1. Virtual environment settings](#a1-virtual-environment-settings)
  - [A2. Project structure, brief description for porting to Java](#a2-project-structure-brief-description-for-porting-to-java)
      * [System kernel: /a_rocksdb](#system-kernel-a_rocksdb)
        - [Auxiliary modules for organizing requests](#auxiliary-modules-for-organizing-requests)
        - [Main modules](#main-modules)
        - [Classes that replace ASchema in specific cases](#classes-that-replace-aschema-in-specific-cases)
        - [Implementation of different logics for grouping records](#implementation-of-different-logics-for-grouping-records)
        - [“Deep data compilation” mode](#deep-data-compilation-mode)
        - [Implementation of different logics for grouping records](#implementation-of-different-logics-for-grouping-records)
      * [Encoding/packing algorithms: /codec](#encodingpacking-algorithms-codec)
      * [The logic for generating a “simplified JSON representation”](#the-logic-for-generating-a-simplified-json-representation)
      * [Other directories and files](#other-directories-and-files)

  
> ## Terminology notes: array -> section -> schema
>
>- below **a section** refers to a data component located in its own RocksDB instance
>- **an array** is a group of **sections** that are processed through 
one request; **The array consists of sections**, and not vice versa
>- Data that came from a single source is controlled by **a schema**:
           gnomAD, GTF, SpliceAI
> 
>In the current implementation, each section contains exactly 
one schema, and *the names of the section and its schema 
usually coincide*, but generally, this is not necessary.


## Service requests
Currently, the REST service is running on the underlying server on port 8290. It returns data on two requests:

- **get** - getting data from storage, arguments:
  * **array** - request array
  * **loc** - position on the chromosome in <chrom>:<pos> format
  * (additional arguments)
    * **alt, ref** for main arrays hg19, hg38
    * **feature** for gtf array
- **meta** - meta information descriptor, no arguments

The following arrays are currently configured:
- **hg19** - data on assembly19: sections Gerp, gnomAD
- **hg38** - data on assembly 38: dbNSFP section (later SpliceAI section will be added)
- **gtf** - data from GTF/Ensembl, according to assembly 38, are placed in a separate array as it is possible not to participate in the main procedure of annotation
- **fasta** - genome assembly type=hg38/hg19

>The get method returns a dictionary with the "chrom", "pos", "array" keys, as well as the names of the data sections. Below the section key is a list of entries from the given section that match the given chrom:pos position.
>
>In case of additional arguments in the query, the lists of entries from the sections ***are filtered*** to match the arguments. (The list in the Gerp section is indifferent to alt/ref arguments).

An example of a request to the “standard” **hg19/hg38** arrays:
```shell
curl "localhost:8290/get?array=hg38&loc=12:885081&alt=G"
```

 
Arrays **gtf** and **fasta** process loc argument not only by individual positions, but also by ranges in format
`chrom:pos1-pos2`.
In the case of a range, instead of the “pos” key, the returned structure contains the “start” and “end” keys.

The **gtf** array works with both single positions and ranges, and returns a list of "windows" that the given position/range intersects with. Filter parameter for gtf array: feature.
Example request to **gtf** array:
```shell
curl "localhost:8290/get?array=gtf&loc=12:885081-985081&feature=exon"
```
The **fasta** array returns a nucleotide or a sequence of nucleotides at the position and range specified in loc, the type=hg19/hg38 argument is required.
An example of a request to the **fasta** array:
```shell
curl "localhost:8290/get?array=fasta&type=hg38&loc=18:67760520-67760529"
```

Service **meta-information** is defined manually in the service configuration.

The `indent=1` argument works for all the listed requests; when this argument is activated, the JSON output is formatted in a readable form, with indents.

## Bulk requests
To speed up the work, a method of downloading from the service at once a lot of information has been implemented. The request has the address “/collect”, and is supposedly called according to the POST scheme. The request supports arguments whose values are in JSON format and are contained in a dictionary:
- **variants**: [array of dictionaries with required fields "chrom", "pos", and optional "ref", "alt"]
- **fasta**: “hg19”/”hg38”; may be omitted, default is “hg38”)
- **arrays**: [array of sections]; may be omitted, all of the following are enabled by default: ["Gerp", "SpliceAI", "gnomAD", "dbNSFP", "gtf"]

>To call the query correctly, you need to do one of two following:
>- either put a JSON representation of the request in the request body, specifying content-type="application/json"
>- or simulate a POST request “from the client” specifying these fields, and the values of the variants and arrays fields (the latter is optional) must be presented in JSON format

Thus, the request “/collect”
- does not directly correspond to ordinary requests, but gives out all the information on the options at once, regardless
- options are processed in either of two builds: hg19/hg38, the service uses the liftover functionality to transfer from one system to another, with both sides; this functionality is not used by the service in routin methods
- through a request, you can get almost all information from the service; exceptions are
  1. information from the fasta array about assembly nucleotides;
  2. information from the gtf array by windows, i.e. by position ranges

    
> To connect to bulk request information on fasta array, the following is required:
> - Add to the dictionaries located in **variants** the field “last” - the last position of the selection of the range for fasta (this parameter works only according to the information on fasta)
> - When explicitly enumerating arrays in **arrays**, use the "fasta/hg19" and "fasta/hg38" pseudo-arrays.

## Service structure
### Directories and files
```
${WORK} — Project root directory (/projects/AStorage)   
|    ${WORK}/astorage.json — Service setup      
|
|——— ${WORK}/rdbs — Sections directory
|    |
|    |———${WORK}/rdbs/<section> — For each section - a directory with a deployed structure RocksDB
|    
|——— ${WORK}/schema — Control directory
     |
     |———${WORK}/schema/<section> — For each section - control information directory
              ${WORK}/schema/<section>/<schema>.json — Control description of the schema in JSON format
              ${WORK}/schema/<section>/<schema>.samples — Schema read/write random checks protocol
```
> As noted above, the title of section and schema are usually the same. However, it is not always. In particular, when/if you want to re-create the section, without disturbing the current working copy, you can create a section with a new name with an existing schema name in it.

### Programs
Module code is placed in Anfisa-Annotations repository in `a_storage` top directory. At the moment all this is only available at the brunch `origin/dev-ingestion`. Currently, to run utilities and services you need to specify in PYTHONPATH the path to deployed anfisa project - it uses utils directory. In the short term, this dependence should be transformed into a ‘wheel’, deployed in a virtual environment.

The launch must be carried out in a specially designed virtual environment; a correctly designed and seemingly accessible environment is located on the server in the `/home/trifon/.forome-venv` directory.

#### Launches:
- **create_db.py** - create/recreate a section:

  ```shell
  python -u create_db.py -c <path to astorage.cfg> -m <partition name>
  ```
  
  The main partition settings are contained in `storage.py` and also in the code in the file `a_storage/ingest/s_<section>.py`

- **service.py** - start the service (standalone or via uWSGI), standalone start:

  ```shell
  python service.py <path to astorage.cfg>
  ```

- internal scripts can be used for intermediate data preparation, this is done for two sections, due to the fact that the data in its original form is very difficult to process (DIR is the key call parameter):
  ```shell
  python -u -m ingest.a_gnomad211 DIR “<base genomes pattern>” “<base exomes pattern> <output directory>  
  python -u -m ingest.a_spliceai DIR “<base pattern>” <output directory>
  ```

## Service setup
### Administrative aspects: number of open files
It is critically important that the logged user who launches the database service and utilities, has the ability to keep open a lot of files, much more than ordinary user.
You can find out the user's limit:
```
ulimit -Sn
```
Must be 10,000 or more.


> To fix the problem below you need to be sudo

To set a limit, it must be written in the `/etc/security/limits.conf` file as a line:
```
<user> soft nofile 20000
```
> (If the parameter `ulimit -Hn` is set below the desired value, you must also add a line to "hard/nofile")



To make these settings work without restarting the computer, you must also add the line 
```
session required pam_limits.so
```
to the files `/etc/pam.d/common-session-noninteractive` and `/etc/pam.d/common-session`
And also in `/etc/ssh/sshd_config` change the value to yes:
```
UsePAM yes
```
### General setup astorage.cfg                            
> The astorage.cfg file is organized as a JSON file where the same features as in the anfisa/anfisa.json project are "working":
> - ‘comments’ - // at the beginning of the line
> - ‘macros’ - defined inside the "file-path-def" field.


The general principle of the file: configuring all aspects at once in all programs of the module:

- Module directory location: 
  * **db-dir** section directory
  * **schema-dir** schema directory, section samples directory
- General configuration of RocksDB bases: 
  * **db-options** opening RocksDB databases options
  * **col-options** general settings of different types of columns (details can be redefined within schema settings; works only when creating section)
- Sections creating settings
  * **samples-count** number of entries on random controlling
  * **create** dictionary, data source options are defined for each section
  (direct_file_list is considered if intermediate data preparation has been performed)
- Service logic setting
  * **service**
    - **arrays** grouping arrays from sections
    - **meta** service metadata, request /meta
  * **host, port** offline service start configuration
  * **dir-files** an empty list of directories the service could travers 
  * **logging** setup logging (note: inside the path to log service!)

### Schema setup
> The general principle: settings are given in JSON format, and  initially, 
> for any scheme they are defined in code in files like `ingest/s_<schema>.py`
No additional software solutions are needed, but nesting is desirable.


 The configuration option defined in this file **only works when you create a section**. When filling in the data, the JSON structure is enriched: the default settings are adjusted and the statistics collected during the fill-in are added, and the result of the upload is stored in a file
`${WORK}/schema/<section>/<schema>.json`

And this unmodifiable file is already used when the service is running, along with the `${WORK}/rdbs/<section>` parallel directory, where the section data is deployed in RocksDB.

Over time, a full description of these settings should appear. Now, see examples of configurations in both the specified files and code and the following explanations:

#### Principles of stacking JSON structures
The general considerations on how data is transformed from JSON column format to RocksDB and back:
- The input JSON-structure is converted into a pair of lines. The second line is required only in the structure has/defined (not a dictionary, see below) string fields:
   - the first line is “primitive JSON”, which contains only arrays, numbers, and the null constant; the effect of compression is expected to be the low variety of characters (like 16 or so...) that can be found in such a file.
   - the second “line” is all (non-dictionary) lines encountered in the structure, combined together with the delimiter '\0'
   - these two lines fit into two different RocksDB columns, of types “base” and (if a second one is needed) “str”

- Floating numbers are “ordinary” numbers, for them a conversion takes place, leaving 4 digits in the mantissa, in Python terms '%.3e', and extra characters are trimmed from the representation. If the number is known to be an integer, it is appropriate to specify the “format” option: ”%d” for the field. (Due to the same “format” option, you can also adjust the number of digits in the mantissa)

- String fields have two options for the “opt” option:
  * “opt”: “dict” means that there are few options for the field value, and they are collected in a dictionary (the dictionary is added to the schema description at the end of the load)
  * "opt": "repeat" means that this string may be repeated frequently, and repeating its value makes sense to track and factorize
- Processing lists in a JSON structure considers that lists consist of elements of the same type
- Dictionaries: The schema description must contain all dictionary key values ​​that need to be considered. The exception is “attribute groups” (an example is in gnomAD): this is not a single key, but a list of keys that can be found in the dictionary, and the values ​​for all these keys are described by one type

- Null can be found everywhere as a field value.

> All this above is implemented in the code in the a_storage/codec directory.

#### Keys and record locking inside RocksDB
In the same directory `a_storage/codec` the logic of working with keys and with block packs is collected.

***Keys***: When using arrays, the key role (literally) plays a couple of parameters: chromosome/position. As part of AStorage for the two assemblies - hg19 and hg38, we implement an internal representation of this pair as a 4-byte number. The chromosome ranges are extended and do not fit together. These internal representations are used as keys in RocksDB.

***Block Packs***: Usually the records do not get into RocksDB individually, they are grouped into blocks long enough for the external archive to be able to compress something. Now, this is done in one of two "and a half" ways:

- **Block packing segments**: the range of all position values is divided into fixed-length segments, and all records from one segment are combined into the block
- **Cluster Block Packing**: Entries are grouped into blocks of relatively equal size, grouping creates an additional data object with a table of entries in the block. Additional column is also created with a page-by-page table of all clusters ( to circumvent the seek() call from RocksDB functionality, which has shown itself to be unstable on large volumes)
- **Packaging frames**: this is a pseudo packaging, it is a way to use the RocksDB seek() call to find from the position a block of records related to the "window" positions, which includes the specified item. It is used with the "gtf" scheme (here we use seek() call, but on relatively small volumes, where behavior is fairly stable)

## Supplements
### A1. Virtual environment settings
In the virtual environment, we need to install our “general” 
forome_utils package, as well as the packages specified 
in `requirements.txt` from the Anfisa-Annotations project

And also the plainrocks package:
>During the development, an unusual brunch of the python-rocksdb project was used. In the process of achieving stability, we moved on to the solution described below

There is `a_storage/plainrocks` directory in the project code. It contains everything we currently need from the C++ rocksDB code, which moves to the Python level through Cython. 

To install it in the project environment, you need:

- Pre-install on the system (```sudo apt-get install```) the librocksdb-dev/disco package (version 5.17.2-3 or higher)
- Activate the correct virtual environment
- Go into the `a_storage/plainrocks` code directory and type:
  ```shell
  python setup.py install
  ```
As a result, the plainrocks package will be installed in the virtual environment.It’s all you need to make AStorage work.

### A2. Project structure, brief description for porting to Java 
> Attention! Actual branch is **dev-ingestion.1**

Below **(\*)** denotes the main files that need a port to Java in the first place.

#### System kernel: /a_rocksdb
- **app.py** - organizing the application as an HTTP server, you don't need to rewrite it directly in Java, you need to do your own tweaking

  ##### Auxiliary modules for organizing requests:

- **a_array.py(\*)** - representation of an array, as a collection of sections of AArray

- **a_collect.py(\*)** - implementation of bulk requests

  ##### Main modules:

- **a_storage.py(\*)** - AStorage head logic module

- **a_schema.py(\*)** - the main representation of the ASchema section

- **a_connector.py(\*)** - all direct work with RocksDB, used from under ASchema

- **a_io.py(\*)** - tool-classes supporting read-write logic

- **a_blocker.py(\*)** - general logic for blocking records within a section

  ##### Classes that replace ASchema in specific cases:

- **a_fasta_schema.py** - Schema implementation for the Fasta section

- **a_seg_schema.py** - composite section glued together from several sections, not used

  ##### Implementation of different logics for grouping records:
- **a_blk_segment.py(\*)** - logic for blocking records at intervals of positions of a fixed size (used for “solid” sections)
- **a_blk_cluster.py** - logic for blocking records in batches of more or less the same size (used for "leaky" partitions)
- **a_blk_pager.py** - internal blocking of the table of contents of "cluster" logic in pages, used to organize access to clusters bypassing the seek () method in RocksDB
- **a_blk_frames.py** - special logic for blocking records for gtf
  ##### “Deep data compilation” mode:
- **deep_comp.py** - write/read archives in deep compilation format
- **deep_load.py** - logic for loading deep compilation archives


#### Encoding/packing algorithms: /codec
> The logic of this directory does not need to be copied in Java one to one. You can use some other solutions instead.

- **\_\_init\_\_.py(\*)** - an essential file through which classes from the directory are obtained for use

- **hg_key.py(\*)** - HG19/HG38 system encoding

- **block_support.py(\*)** - packing/formatting algorithms used when blocking records

##### The logic for generating a “simplified JSON representation”:
- **_codec_data.py**
- **codec_num.py**
- **codec_str.py**
- **codec_dict.py**
- **codec_list.py**
- **codec_agroup.py**

#### Other directories and files
- **/ingest** - logic for loading various data sections, a porting to Java is probably not needed at all
- **/plainrocks** - builds a python package for RocksDB, doesn't make sense at all for a Java porting
- **create_db.py(\*)** - general partition creation script
- **service.py** - HTTP application launcher, doesn't make sense for a porting to Java
