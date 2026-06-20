############################################################
## Selection scan of Plasmodium falciparum populations
## Rwanda 2021–2025
##
## Generates:
##   - Per-SNP FST
##   - Allele frequency differences
##   - Signed allele frequency changes
##   - Manhattan plots
############################################################

rm(list=ls())

############################################################
## Load libraries
############################################################

library(vcfR)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(viridis)
library(scales)
library(hierfstat)
library(adegenet)

theme_set(theme_bw())

############################################################
## Input files
############################################################

vcf_file <- "final_population.ann.vcf"

metadata_file <- "samples_metadata_final_07-04-2026.tsv"

annotation_file <- "population_variant_table.tsv"
output_dir <- "Selection_scan"

if(!dir.exists(output_dir))
    dir.create(output_dir)
############################################################
## Load annotation table
############################################################

annot_raw <- read.csv(
    annotation_file,
    sep="\t",
    stringsAsFactors=FALSE,
    na.strings=c("",".","NA")
)

cat("Annotation loaded\n")
cat(nrow(annot_raw),"variants\n")

############################################################
## Build SNP annotation table
############################################################

annot_raw <- annot_raw %>%

mutate(

SNP_ID = paste(CHROM, POS, sep="_"),

Gene = GENE_NAME,

Protein = HGVS_PROTEIN

)
############################################################
## Keep unique SNP annotations
############################################################

annot_table <- annot_raw %>%

select(

SNP_ID,
CHROM,
POS,
Gene,
Protein

) %>%

distinct()
############################################################
## Create plotting labels
############################################################

annot_table <- annot_table %>%

mutate(

LABEL = ifelse(

is.na(Protein),

Gene,

paste0(Gene, ":", Protein)

)

)
############################################################
## Check annotation
############################################################

cat("\n")

cat("Annotated SNPs:",
    nrow(annot_table),
    "\n")

cat("\nGenes represented:\n")

print(
    sort(
        table(annot_table$Gene),
        decreasing=TRUE
    )
)
