# All file paths for the R analyses, in one place.
# Run R from the repository root, then: source("config.R")

PROJECT_DIR <- normalizePath(".", mustWork = FALSE)

# --- Inputs tracked in this repository ---
SAMPLE_METADATA <- file.path(PROJECT_DIR, "data", "metadata", "samples_metadata.tsv")
SEQUENCED_DATA  <- file.path(PROJECT_DIR, "data", "metadata", "rwanda_sequenced_data.xlsx")
FLUO_DIR        <- file.path(PROJECT_DIR, "data", "fluorescence_assays")

# --- Processed tables ---
PROCESSED_DIR        <- file.path(PROJECT_DIR, "data", "processed")
POPULATION_VCF       <- file.path(PROCESSED_DIR, "final_population.ann.vcf")
POPULATION_VAR_TABLE <- file.path(PROCESSED_DIR, "population_variant_table.tsv")
DRUG_MARKER_TABLE    <- file.path(PROCESSED_DIR, "drug_marker_table.tsv")

# --- Outputs ---
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIGURE_DIR  <- file.path(RESULTS_DIR, "figures")
TABLE_DIR   <- file.path(RESULTS_DIR, "tables")

for (d in c(PROCESSED_DIR, RESULTS_DIR, FIGURE_DIR, TABLE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FIG_WIDTH <- 10; FIG_HEIGHT <- 7; FIG_DPI <- 300

# --- stage 01: merged per-sample FASTQ files ---
MERGED_DIR <- file.path(PROJECT_DIR, "data", "fastq_merged")

# Selection-scan inputs (Figure 5)
SELECTION_VCF <- "/Volumes/MTD_drive/Manuscript_review_analysis_04-08-2026/final.vcf.gz"
SELECTION_ANNOTATION <- "/Volumes/MTD_drive/Manuscript_review_analysis_04-08-2026/snpeff_annotation_extraction/final_annotation/annotated_variants_snpeff_final.tsv"

# Figure 1: GADM shapefiles (Rwanda plus DR Congo, Uganda, Tanzania, Burundi)
GADM_DIR <- "~/Downloads"

# Supplementary Figure S1b: monthly malaria case counts
MALARIA_CASES_FILE <- "/Volumes/MTD_drive/Manuscript_review_analysis_04-08-2026/Malaria cases Kigali 2023vs2024.xlsx"
