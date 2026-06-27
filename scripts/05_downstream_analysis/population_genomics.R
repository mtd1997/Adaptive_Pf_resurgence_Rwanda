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

############################################################
# PCoA analysis
############################################################

period_colors <- c(
  "Pre-Resurgence" = "#1f77b4",
  "Post-Resurgence" = "#d62728"
)

region_colors <- c(
  "Kigali" = "#6B7900",
  "Hotspot" = "#D16100"
)

k13_colors <- c(
  "WT" = "#2ecc71",
  "K189T_only" = "#d62f7d",
  "R561H+K189T" = "#ff7f0e",
  "R561H_only" = "#f4b400",
  "Other_K13" = "#6a5acd"
)

period_shapes <- c(
  "Pre-Resurgence" = 16,
  "Post-Resurgence" = 17
)

region_shapes <- c(
  "Kigali" = 16,
  "Hotspot" = 17
)

############################################################
# Format metadata
############################################################

samples_info <- samples_info %>%
  mutate(
    period = factor(period),
    region = factor(region),
    K13_status = factor(K13_status)
  )

############################################################
# PCoA function
############################################################

compute_pcoa <- function(gt_matrix, metadata) {

  genind_obj <- df2genind(
    gt_matrix,
    ploidy = 1,
    type = "PA"
  )

  dist_mat <- poppr::prevosti.dist(
    genind_obj
  )

  pcoa <- ape::pcoa(
    dist_mat,
    correction = "cailliez"
  )

  eig <- pcoa$values$Eigenvalues

  var_exp <- eig[eig > 0] /
    sum(eig[eig > 0]) * 100

  df <- as.data.frame(
    pcoa$vectors[,1:10]
  )

  colnames(df) <- paste0(
    "Axis",
    1:10
  )

  df$sample_id <- rownames(gt_matrix)

  df <- left_join(
    df,
    metadata,
    by = c("sample_id" = "SAMPLE")
  )

  list(
    df = df,
    var_exp = var_exp
  )
}

############################################################
# Plot function
############################################################

plot_panels <- function(df,
                        x1, x2, x3, x4,
                        var_exp,
                        color_var,
                        shape_var = NULL,
                        colors,
                        shapes = NULL,
                        filename) {

  make_plot <- function(x, y, xlab, ylab) {

    g <- ggplot(
      df,
      aes(.data[[x]], .data[[y]])
    ) +

      geom_point(
        aes(
          color = .data[[color_var]],
          shape = if (!is.null(shape_var))
            .data[[shape_var]]
        ),
        size = 2,
        alpha = 0.9
      ) +

      stat_ellipse(
        aes(color = .data[[color_var]]),
        type = "norm",
        linewidth = 1
      ) +

      scale_color_manual(
        values = colors
      ) +

      theme_classic(
        base_size = 16
      ) +

      labs(
        x = xlab,
        y = ylab
      )

    if (!is.null(shape_var)) {
      g <- g +
        scale_shape_manual(
          values = shapes
        )
    }

    g
  }

  p1 <- make_plot(
    x1, x2,
    paste0("PCoA1 (", round(var_exp[1],1), "%)"),
    paste0("PCoA2 (", round(var_exp[2],1), "%)")
  )

  p2 <- make_plot(
    x1, x3,
    paste0("PCoA1 (", round(var_exp[1],1), "%)"),
    paste0("PCoA3 (", round(var_exp[3],1), "%)")
  )

  p3 <- make_plot(
    x2, x3,
    paste0("PCoA2 (", round(var_exp[2],1), "%)"),
    paste0("PCoA3 (", round(var_exp[3],1), "%)")
  )

  p4 <- make_plot(
    x1, x4,
    paste0("PCoA1 (", round(var_exp[1],1), "%)"),
    paste0("PCoA4 (", round(var_exp[4],1), "%)")
  )

  final_plot <- (p1 | p2) /
                (p4 | p3)

  ggsave(
    filename,
    final_plot,
    width = 16,
    height = 11,
    dpi = 600
  )
}

############################################################
# Kigali subset
############################################################

meta_kigali <- samples_info %>%
  filter(region == "Kigali")

gt_kigali <- geno_matrix[
  meta_kigali$SAMPLE,
]

pcoa_kigali <- compute_pcoa(
  gt_kigali,
  meta_kigali
)

plot_panels(
  pcoa_kigali$df,
  "Axis1","Axis2","Axis3","Axis4",
  pcoa_kigali$var_exp,
  "period",
  NULL,
  period_colors,
  NULL,
  "PCoA_Kigali_period.svg"
)

plot_panels(
  pcoa_kigali$df,
  "Axis1","Axis2","Axis3","Axis4",
  pcoa_kigali$var_exp,
  "K13_status",
  "period",
  k13_colors,
  period_shapes,
  "PCoA_Kigali_K13.svg"
)

############################################################
# Post-resurgence subset
############################################################

meta_post <- samples_info %>%
  filter(period == "Post-Resurgence")

gt_post <- geno_matrix[
  meta_post$SAMPLE,
]

pcoa_post <- compute_pcoa(
  gt_post,
  meta_post
)

plot_panels(
  pcoa_post$df,
  "Axis1","Axis2","Axis3","Axis4",
  pcoa_post$var_exp,
  "region",
  NULL,
  region_colors,
  NULL,
  "PCoA_Post_region.svg"
)

plot_panels(
  pcoa_post$df,
  "Axis1","Axis2","Axis3","Axis4",
  pcoa_post$var_exp,
  "K13_status",
  "region",
  k13_colors,
  region_shapes,
  "PCoA_Post_K13.svg"
)

cat("PCoA analyses completed.\n")
