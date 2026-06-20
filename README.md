# Adaptive restructuring of *Plasmodium falciparum* populations underlies malaria resurgence in Rwanda

## Code repository


## Overview
```text
This repository contains the computational workflow used to process Oxford Nanopore amplicon sequencing data, perform variant calling, annotate variants, and perform population genomic and antimalarial resistance analyses described in the associated manuscript.
```
---

## Data availability

```text
The raw sequencing data (FASTQ files) and variant files (VCF and gVCF) generated in this study have been deposited in the Gene Expression Omnibus (GEO) repository with accession number GSE318440 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318440.
Because of their size, sequencing data are not included in this repository.
```
---

# Repository structure

```text
scripts/

01_preprocessing/
02_alignment/
03_variant_calling/
04_filtering/
05_downstream_analysis/
```

---

# Computational workflow

1. Basecalling and demultiplexing using MinKNOW.
2. FASTQ concatenation and sample renaming.
3. Quality assessment using NanoPlot and pycoQC.
4. Alignment against the *Plasmodium falciparum* 3D7 reference genome using minimap2.
5. BAM processing using samtools.
6. Variant calling using Clair3.
7. Joint genotyping using GATK GenomicsDBImport and GenotypeGVCFs.
8. Variant filtering using bcftools and vcftools.
9. Variant annotation using SnpEff.
10. Downstream population genomic and drug resistance analyses using Python and R.

---

# Software versions

| Software | Version                          |
| -------- | -------------------------------- |
| MinKNOW  | version used for sequencing runs |
| NanoPlot | 1.42.0                           |
| pycoQC   | version used in analysis         |
| minimap2 | 2.28                             |
| samtools | 1.20                             |
| bedtools | 2.31                             |
| Clair3   | latest Docker release (v1.0.10)  |
| GATK     | 4.6.1.0                          |
| bcftools | 1.20                             |
| vcftools | 0.1.16                           |
| SnpEff   | 5.2                              |
| Python   | 3.x                              |
| R        | 4.x                              |

---

# Reference genome

*Plasmodium falciparum* 3D7 reference genome (PlasmoDB release 68).

---

# Variant annotation

Functional annotation was performed using SnpEff configured with the *Plasmodium falciparum* 3D7 database.



