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
############################################################
## Load metadata
############################################################

meta <- read_tsv(
    metadata_file,
    show_col_types = FALSE
)

meta <- meta %>%
    rename(sample_id = SAMPLE)

############################################################
## Load VCF
############################################################

vcf <- read.vcfR(
    vcf_file,
    verbose = FALSE
)

############################################################
## Extract genotype matrix
############################################################

gt <- extract.gt(
    vcf,
    element = "GT",
    as.numeric = TRUE
)

gt <- t(gt)

############################################################
## Replace missing calls
############################################################

gt[is.na(gt)] <- 0

############################################################
## Match metadata
############################################################

common_ids <- intersect(
    rownames(gt),
    meta$sample_id
)

gt <- gt[common_ids, ]

meta <- meta %>%
    filter(sample_id %in% common_ids)

meta <- meta[
    match(common_ids, meta$sample_id),
]

############################################################
## Check
############################################################

stopifnot(
    all(
        rownames(gt) ==
        meta$sample_id
    )
)

cat(
    "Samples:",
    nrow(gt),
    "\n"
)

cat(
    "Variants:",
    ncol(gt),
    "\n"
))
############################################################
## Build SNP position table
############################################################

snp_positions <- annot_table %>%

    select(
        SNP_ID,
        CHROM,
        POS,
        Gene
    ) %>%

    distinct()

############################################################
## Gene order
############################################################

gene_order <- c(
    "MSP2",
    "CSP",
    "DHFR",
    "MDR1",
    "CRT",
    "DHPS",
    "K13"
)

############################################################
## Create cumulative positions
############################################################

gap <- 5000

snp_positions <- snp_positions %>%

    mutate(
        Gene = factor(
            Gene,
            levels = gene_order
        )
    ) %>%

    arrange(
        Gene,
        POS
    ) %>%

    group_by(Gene) %>%

    mutate(
        rel_pos = POS - min(POS)
    ) %>%

    ungroup()

gene_offsets <- snp_positions %>%

    group_by(Gene) %>%

    summarise(
        span = max(rel_pos) + gap,
        .groups = "drop"
    ) %>%

    mutate(
        offset = cumsum(
            lag(span, default = 0)
        )
    )

snp_positions <- snp_positions %>%

    left_join(
        gene_offsets,
        by = "Gene"
    ) %>%

    mutate(
        cum_pos = rel_pos + offset
    )

gene_centers <- snp_positions %>%

    group_by(Gene) %>%

    summarise(
        center = mean(cum_pos),
        .groups = "drop"
    )
############################################################
## Per-SNP FST
## Pre vs Post resurgence
############################################################

compute_fst <- function(
    gt_matrix,
    metadata
) {

    idx_pre <- metadata$period ==
        "Pre-Resurgence"

    idx_post <- metadata$period ==
        "Post-Resurgence"

    fst_values <- sapply(
        seq_len(ncol(gt_matrix)),
        function(i) {

            g1 <- gt_matrix[idx_pre, i]

            g2 <- gt_matrix[idx_post, i]

            p1 <- mean(g1)

            p2 <- mean(g2)

            if(
                (p1 == 0 & p2 == 0) |
                (p1 == 1 & p2 == 1)
            ) {
                return(0)
            }

            pbar <- (
                length(g1)*p1 +
                length(g2)*p2
            ) /
            (
                length(g1)+
                length(g2)
            )

            num <- (p1-p2)^2

            den <- pbar*(1-pbar)

            if(
                den == 0
            ) return(NA)

            num/den
        }
    )

    data.frame(
        SNP_ID = colnames(gt_matrix),
        fst = fst_values
    )
}

fst_results <- compute_fst(
    gt,
    meta
)
############################################################
## Allele frequency shift
## Pre-Resurgence vs Post-Resurgence
############################################################

compute_delta_af <- function(
    gt_matrix,
    metadata
){

    idx_pre <- metadata$period == "Pre-Resurgence"

    idx_post <- metadata$period == "Post-Resurgence"

    freq_pre <- apply(
        gt_matrix[idx_pre, , drop = FALSE],
        2,
        mean,
        na.rm = TRUE
    )

    freq_post <- apply(
        gt_matrix[idx_post, , drop = FALSE],
        2,
        mean,
        na.rm = TRUE
    )

    data.frame(

        SNP_ID = colnames(gt_matrix),

        AF_Pre = freq_pre,

        AF_Post = freq_post,

        Delta_AF = freq_post - freq_pre,

        Abs_Delta_AF = abs(freq_post - freq_pre),

        stringsAsFactors = FALSE

    )

}

delta_results <- compute_delta_af(
    gt,
    meta
)

############################################################
## Join annotation
############################################################

delta_results <- delta_results %>%

    left_join(
        annot_table,
        by = "SNP_ID"
    ) %>%

    left_join(
        snp_positions,
        by = c(
            "SNP_ID",
            "CHROM",
            "POS",
            "Gene"
        )
    )

############################################################
## Define empirical threshold
############################################################

delta_threshold <- quantile(

    delta_results$Abs_Delta_AF,

    0.95,

    na.rm = TRUE

)

delta_results <- delta_results %>%

    mutate(

        Outlier =

            Abs_Delta_AF >= delta_threshold

    )

############################################################
## Summary
############################################################

cat(

    "\nAbsolute allele frequency threshold:",

    round(delta_threshold,3),

    "\n"

)

cat(

    "Number of outlier SNPs:",

    sum(delta_results$Outlier),

    "\n"

)
############################################################
## Signed allele frequency shift
############################################################

signed_delta <- delta_results %>%

    mutate(

        Direction = case_when(

            Delta_AF > 0 ~ "Increase",

            Delta_AF < 0 ~ "Decrease",

            TRUE ~ "No change"

        )

    )

############################################################
## Summary
############################################################

cat("\n")

cat(

    "Allele frequency increases:",

    sum(

        signed_delta$Direction == "Increase"

    ),

    "\n"

)

cat(

    "Allele frequency decreases:",

    sum(

        signed_delta$Direction == "Decrease"

    ),

    "\n"

)

############################################################
## Export results
############################################################

write.csv(

    signed_delta,

    file.path(
        output_dir,
        "Allele_frequency_shift.csv"
    ),

    row.names = FALSE

)
############################################################
## Manhattan plot function
############################################################

plot_manhattan <- function(data,
                           value,
                           ylab,
                           outfile,
                           highlight = TRUE){

    ymax <- max(
        data[[value]],
        na.rm = TRUE
    ) * 1.10

    p <- ggplot(
        data,
        aes(
            x = cum_pos,
            y = .data[[value]]
        )
    ) +

        geom_point(
            aes(color = Gene),
            size = 2
        ) +

        scale_color_brewer(
            palette = "Dark2"
        ) +

        scale_x_continuous(
            breaks = gene_centers$center,
            labels = gene_centers$Gene
        ) +

        labs(
            x = "Target gene",
            y = ylab
        ) +

        coord_cartesian(
            ylim = c(0, ymax)
        ) +

        theme_bw(base_size = 14) +

        theme(

            legend.position = "none",

            panel.grid.major.x = element_blank(),

            panel.grid.minor = element_blank(),

            axis.text.x = element_text(
                angle = 45,
                hjust = 1
            )

        )

    if(highlight){

        p <- p +

            geom_text_repel(

                data = subset(
                    data,
                    Outlier
                ),

                aes(
                    label = LABEL
                ),

                size = 3,

                max.overlaps = Inf

            )

    }

    ggsave(

        file.path(
            output_dir,
            outfile
        ),

        p,

        width = 10,

        height = 4,

        dpi = 600

    )

    return(p)

}
############################################################
## FST Manhattan
############################################################

fst_plot_data <- fst_results %>%

    left_join(
        annot_table,
        by = "SNP_ID"
    ) %>%

    left_join(
        snp_positions,
        by = c(
            "SNP_ID",
            "CHROM",
            "POS",
            "Gene"
        )
    )

fst_threshold <- quantile(

    fst_plot_data$fst,

    0.95,

    na.rm = TRUE

)

fst_plot_data <- fst_plot_data %>%

    mutate(

        Outlier =

            fst >= fst_threshold

    )

plot_manhattan(

    fst_plot_data,

    value = "fst",

    ylab = expression(F[ST]),

    outfile = "FST_manhattan.svg"

)
############################################################
## Absolute allele frequency shift
############################################################

plot_manhattan(

    delta_results,

    value = "Abs_Delta_AF",

    ylab = expression("|" * Delta * "AF|"),

    outfile = "Allele_frequency_shift.svg"

)
############################################################
## Signed allele frequency shifts
############################################################

signed_plot <- ggplot(

    signed_delta,

    aes(

        cum_pos,

        Delta_AF,

        color = Gene

    )

) +

    geom_hline(

        yintercept = 0,

        linetype = 2,

        colour = "grey40"

    ) +

    geom_point(

        size = 2

    ) +

    geom_text_repel(

        data = subset(

            signed_delta,

            Outlier

        ),

        aes(

            label = LABEL

        ),

        size = 3,

        max.overlaps = Inf

    ) +

    scale_color_brewer(

        palette = "Dark2"

    ) +

    scale_x_continuous(

        breaks = gene_centers$center,

        labels = gene_centers$Gene

    ) +

    theme_bw(base_size = 14) +

    theme(

        legend.position = "none",

        panel.grid.major.x = element_blank(),

        panel.grid.minor = element_blank(),

        axis.text.x = element_text(

            angle = 45,

            hjust = 1

        )

    ) +

    labs(

        x = "Target gene",

        y = expression(Delta * "AF")

    )

ggsave(

    file.path(

        output_dir,

        "Signed_allele_frequency_shift.svg"

    ),

    signed_plot,

    width = 10,

    height = 4,

    dpi = 600

)

############################################################
## Export summary tables
############################################################

write.csv(

    fst_plot_data,

    file.path(
        output_dir,
        "FST_results.csv"
    ),

    row.names = FALSE

)

write.csv(

    signed_delta,

    file.path(
        output_dir,
        "Allele_frequency_shift_results.csv"
    ),

    row.names = FALSE

)

cat("\nSelection scan completed successfully.\n")
