#Usage and Examples
=====

## Table of Contents
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
- [Before use](#before-use)
- [taxonkit](#taxonkit)
- [list](#list)
- [lineage](#lineage)
- [reformat](#reformat)
- [name2taxid](#name2taxid)
- [taxid-changelog](#taxid-changelog)
- [genuatocomplete](#genuatocomplete)
<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Before use
Please download and uncompress these files:
- [ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)<br>
<p align="left">And copy "names.dmp", "nodes.dmp", "delnodes.dmp" and "merged.dmp" to data directory: "$HOME/.taxonkit".</p>

## taxonkit
TaxonKit - A Cross-platform and Efficient NCBI Taxonomy Toolkit<br>
Version: 0.5.0<br>
Author: Wei Shen <shenwei356@gmail.com><br>
Source code: [https://github.com/shenwei356/taxonkit](https://github.com/shenwei356/taxonkit)<br>
Documents  : [https://bioinf.shenwei.me/taxonkit](https://bioinf.shenwei.me/taxonkit)<br>
'''
Dataset:

    Please download and decompress "taxdump.tar.gz":
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

    and copy "names.dmp", "nodes.dmp", "delnodes.dmp" and "merged.dmp" to data directory:
    "/home/shenwei/.taxonkit"

    or some other directory, and later you can refer to using flag --data-dir,
    or environment variable TAXONKIT_DB

Usage:
  taxonkit [command]

Available Commands:
  genautocomplete generate shell autocompletion script
  help            Help about any command
  lineage         query lineage of given taxids
  list            list taxon tree of given taxids
  name2taxid      query taxid by taxon scientific name
  reformat        reformat lineage
  taxid-changelog create taxid changelog from dump archives
  version         print version information and check for update

Flags:
      --data-dir string   directory containing nodes.dmp and names.dmp (default "/home/shenwei/.taxonkit")
  -h, --help              help for taxonkit
      --line-buffered     use line buffering on output, i.e., immediately writing to stdin/file for every line of output
  -o, --out-file string   out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
  -j, --threads int       number of CPUs. 2 is enough (default value: 1 for single-CPU PC, 2 for others) (default 2)
      --verbose           print verbose information

Use "taxonkit [command] --help" for more information about a command
'''

