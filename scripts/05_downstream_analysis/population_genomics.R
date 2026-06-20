############################################################
# Population genomics analyses
# Adaptive restructuring of Plasmodium falciparum populations
# during malaria resurgence in Rwanda
############################################################

############################################################
# Load packages
############################################################

library(vcfR)
library(adegenet)
library(poppr)
library(ape)
library(vegan)
library(hierfstat)

library(dplyr)
library(tidyr)
library(readr)

library(ggplot2)
library(ggtree)
library(scales)

############################################################
# Input files
############################################################

vcf_file <- "final_population.ann.vcf"

metadata_file <- "samples_metadata_final_07-04-2026.tsv"

############################################################
# Import metadata
############################################################

samples_info <- read_tsv(
    metadata_file,
    show_col_types = FALSE
)

############################################################
# Import VCF
############################################################

vcf <- read.vcfR(
    vcf_file,
    verbose = FALSE
)

############################################################
# Convert to genlight object
############################################################

genlight <- vcfR2genlight(vcf)

############################################################
# Remove reference strain
############################################################

genlight <- genlight[
    indNames(genlight) != "3D7",
]

############################################################
# Match metadata
############################################################

samples_info <- samples_info %>%
    filter(
        SAMPLE %in% indNames(genlight)
    )

samples_info <- samples_info[
    match(
        indNames(genlight),
        samples_info$SAMPLE
    ),
]

############################################################
# Check sample order
############################################################

stopifnot(
    all(
        indNames(genlight) ==
            samples_info$SAMPLE
    )
)

############################################################
# Build genotype matrix
############################################################

geno_matrix <- tab(
    genlight,
    NA.method = "mean"
)

############################################################
# Summary
############################################################

cat("\n")

cat("Number of samples:",
    nrow(geno_matrix),
    "\n")

cat("Number of loci:",
    ncol(geno_matrix),
    "\n")

cat("\n")

```r
############################################################
# Neighbor-Joining phylogenetic analysis
############################################################

############################################################
# Select Kigali samples
############################################################

meta_kigali <- samples_info %>%
    filter(region == "Kigali")

gen_kigali <- genlight[
    indNames(genlight) %in% meta_kigali$SAMPLE,
]

############################################################
# Genotype matrix
############################################################

geno_mat <- tab(
    gen_kigali,
    NA.method = "mean"
)

############################################################
# Compute Prevosti genetic distance
############################################################

genind_obj <- df2genind(
    geno_mat,
    ploidy = 1,
    type = "PA"
)

dist_prevosti <- poppr::prevosti.dist(
    genind_obj
)

############################################################
# Neighbor-Joining tree
############################################################

tree <- ape::nj(
    dist_prevosti
)

############################################################
# Match metadata to tree tips
############################################################

meta_kigali <- meta_kigali[
    match(
        tree$tip.label,
        meta_kigali$SAMPLE
    ),
]

############################################################
# Plot tree
############################################################

tree_plot <- ggtree(
    tree,
    layout = "fan"
) %<+% meta_kigali +

geom_tippoint(

    aes(
        colour = period
    ),

    size = 2.5

) +

scale_colour_manual(
    values = period_colors
) +

theme_tree2()

############################################################
# Export tree
############################################################

ggsave(
    filename = "Figure3_NJ_tree.svg",
    plot = tree_plot,
    width = 8,
    height = 8,
    dpi = 600
)

ggsave(
    filename = "Figure3_NJ_tree.pdf",
    plot = tree_plot,
    width = 8,
    height = 8
)

write.tree(
    tree,
    file = "Figure3_NJ_tree.nwk"
)

############################################################
# Export iTOL annotation: period
############################################################

itol_period <- c(

    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tPeriod",
    "COLOR\t#000000",
    "STRIP_WIDTH\t25",
    "MARGIN\t5",
    "SHOW_INTERNAL\t0",
    "DATA",

    paste(
        meta_kigali$SAMPLE,
        period_colors[meta_kigali$period],
        sep = "\t"
    )

)

writeLines(
    itol_period,
    "itol_period.txt"
)

############################################################
# Export iTOL annotation: K13 genotype
############################################################

itol_k13 <- c(

    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tK13",
    "COLOR\t#000000",
    "STRIP_WIDTH\t25",
    "MARGIN\t5",
    "SHOW_INTERNAL\t0",
    "DATA",

    paste(
        meta_kigali$SAMPLE,
        k13_colors[meta_kigali$K13_status],
        sep = "\t"
    )

)

writeLines(
    itol_k13,
    "itol_K13.txt"
)

cat("Neighbor-Joining tree completed.\n")
