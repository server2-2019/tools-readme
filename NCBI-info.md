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
- **abnormal gene to sequence ratio** - the ratio of the number of predicted genes to the length of the genome divided by 1000 is far outside the usual range for a Complete Genome assembly. The NCBI Prokaryotic Genome Annotation Pipeline typically expects to find an average of one gene for every 1,000 nucleotides in a genome assembly. The typical range is 0.8 to1.2; anything outside the range 0.5 to 1.5 is considered abnormal.
- **derived from environmental source** - the source material for the assembly is from an environmental source rather than a pure culture leading to concern about the accuracy of organism assignment and possible cross-contamination.
- **derived from metagenome** - the genomic sequence was assembled from metagenomic sequencing rather than a pure culture leading to concerns about the accuracy of organism assignment and possible cross-contamination.
- **derived from single cell** - the source material for the assembly was amplified from a single cell leading to concern about the genome sequence accuracy.
- **derived from surveillance project** - assembly generated from a pathogen surveillance project.
- **fragmented assembly** - a prokaryotic assembly with contig L50 above 500, contig N50 below 5000, or more than 2,000 contigs.
- **genome length too large** - total non-gapped sequence length of the assembly is more than 1.5 times that of the average for the genomes in the Assembly resource from the same species, more than 15 Mbp, or is otherwise suspiciously long.
- **genome length too small** - total non-gapped sequence length of the assembly is less than half that of the average for the genomes in the Assembly resource from the same species, less than 300 Kbp, or is otherwise suspiciously short.
- **low gene count** - the number of predicted genes is much lower than expected.
- **low quality sequence** - long stretches of the sequence have a high proportion of ambiguous bases, are low complexity, or have some other indication that the sequence quality is low.
- **many frameshifted proteins** - the CDSs predicted by the NCBI Prokaryotic Genome Annotation Pipeline have a suspiciously high number of frameshifts. For any clade containing at least 10 good quality assemblies the cutoff is more than three standard deviations from average or 5% of annotated CDSs, whichever is larger. For any clade containing less than 10 good quality assemblies the cutoff is more than 30% of total CDSs.
- **metagenome** - a metagenome containing sequences from a mixture of organisms.
- **missing ribosomal protein genes - the NCBI Prokaryotic Genome Annotation Pipeline failed to find at least one copy of each essential ribosomal protein.
- **missing rRNA genes** - the NCBI Prokaryotic Genome Annotation Pipeline failed to find at least one copy each of 5S, 16S, and 23S rRNA gene.
- **missing tRNA genes** - the NCBI Prokaryotic Genome Annotation Pipeline failed to find tRNA genes with anticodons for 2 or more of the expected 20 amino-acids.
- **partial** - the assembly has only partial genome representation.
- **RefSeq annotation failed** - RefSeq failed to produce valid annotation.
- **untrustworthy as type** - does not meet the criteria for a type strain assembly for use in ANI analysis.
***
