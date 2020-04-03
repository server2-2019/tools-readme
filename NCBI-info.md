# NCBI info

## Table of Contents
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
- [Before use](#before-use)
- [Assembly-genome-refseq](#Assembly-genome-refseq)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Before use - Links
Please download and uncompress these files:
- [Assembly-genome-refseq](https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/)<br>
<p align="left">Assembly Anomalies and Other Reasons a Genome Assembly may be Excluded from RefSeq.</p>

## Assembly-genome-refseq
Assembly Anomalies and Other Reasons a Genome Assembly may be Excluded from RefSeq<br>
Source websie: [https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/](https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/)<br>
### Introduction
```
A genome assembly may be excluded from the NCBI Reference Sequence (RefSeq) project for any of the reasons listed below. 
A few of these reasons for exclusion are likely to be of concern to all users, hence, we report these as "Assembly anomalies".

If the reason for exclusion from RefSeq would also impact your intended use of the genome assembly, you can apply filters to remove 
such assemblies from the search results (these filters appear in the left-hand sidebar under the heading "Exclude", filters for 
more reasons can be exposed using the "Customize …" menu). Filters can be applied for each individual reason, or the "excluded 
from RefSeq" property available under the "Advanced" search menu can be used to filter out all the genome assemblies that were 
excluded from RefSeq.

A filter to remove all assemblies flagged as anomalous from the search results in on by default. This filter must be cleared if 
you want to have assemblies flagged as anomalous included in the search results or if you want to filter out some types of anomaly 
but not others.

"Excluded from RefSeq" reasons can be reported for genome assemblies from any taxonomic group, however, many of the reasons 
for excluding an assembly from RefSeq are only revealed by running the NCBI Prokaryotic Genome Annotation Pipeline (PGAP).
```
### Reasons an assembly is excluded from RefSeq
- *abnormal gene to sequence ratio* - the ratio of the number of predicted genes to the length of the genome divided by 1000 is far outside the usual range for a Complete Genome assembly. The NCBI Prokaryotic Genome Annotation Pipeline typically expects to find an average of one gene for every 1,000 nucleotides in a genome assembly. The typical range is 0.8 to1.2; anything outside the range 0.5 to 1.5 is considered abnormal.
***
