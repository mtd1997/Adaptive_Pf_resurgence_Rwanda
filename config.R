# All file paths for the R analyses, in one place.
#
# Paths under data/ and results/ are relative to the repository. The large
# inputs are not in the repository; set them in config.local.R before running
# anything. See config.local.R.example and the README.

PROJECT_DIR <- normalizePath(".", mustWork = FALSE)

# ---- inputs tracked in this repository ----
SAMPLE_METADATA <- file.path(PROJECT_DIR, "data", "metadata", "samples_metadata.tsv")
SEQUENCED_DATA  <- file.path(PROJECT_DIR, "data", "metadata", "rwanda_sequenced_data.xlsx")
FLUO_DIR        <- file.path(PROJECT_DIR, "data", "fluorescence_assays")

# ---- processed tables ----
PROCESSED_DIR        <- file.path(PROJECT_DIR, "data", "processed")
POPULATION_VCF       <- file.path(PROCESSED_DIR, "final_population.ann.vcf")
POPULATION_VAR_TABLE <- file.path(PROCESSED_DIR, "population_variant_table.tsv")
DRUG_MARKER_TABLE    <- file.path(PROCESSED_DIR, "drug_marker_table.tsv")
COVERAGE_TABLE       <- file.path(PROCESSED_DIR, "coverage_all_samples.tsv")

# One row per sample and mutation, de-identified. Figures 2b, 2c, 3a, 3b, 3d
# and S3 are generated from this file, which is tracked in the repository.
MANUSCRIPT_DATA <- Sys.getenv(
  "PF_MANUSCRIPT_DATA",
  unset = file.path(PROCESSED_DIR, "Manuscript_data_29-07-2026.tsv")
)

# ---- outputs ----
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIGURE_DIR  <- file.path(RESULTS_DIR, "figures")
TABLE_DIR   <- file.path(RESULTS_DIR, "tables")

for (d in c(PROCESSED_DIR, RESULTS_DIR, FIGURE_DIR, TABLE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FIG_WIDTH <- 10; FIG_HEIGHT <- 7; FIG_DPI <- 300

# ---- stage 01: merged per-sample FASTQ files ----
MERGED_DIR <- file.path(PROJECT_DIR, "data", "fastq_merged")

# ---- large inputs, not in the repository ----
# Each may be set by an environment variable or, more simply, in
# config.local.R, which is read below and is not tracked by git.

UNSET <- "<SET_THIS>"
from_env <- function(name) Sys.getenv(name, unset = UNSET)

# Selection scan, Figure 5.
SELECTION_VCF        <- from_env("PF_SELECTION_VCF")
SELECTION_ANNOTATION <- from_env("PF_SELECTION_ANNOTATION")

# Population VCF for the IBS distance matrix, Figure 4.
POPGEN_VCF <- from_env("PF_POPGEN_VCF")

# Monthly malaria case counts, Supplementary Figure S1b.
MALARIA_CASES_FILE <- from_env("PF_MALARIA_CASES")

# Saved Kigali tree the Figure 4d script checks itself against.
FIGURE04D_REFERENCE_TREE <- from_env("PF_REFERENCE_TREE")

# GADM shapefiles for Figures 1, 3d and S1a. Download the Rwanda level-2 file
# and the level-0 files for COD, UGA, TZA and BDI from gadm.org and unzip them
# into one directory. They are not redistributed with this repository.
GADM_DIR <- Sys.getenv("GADM_DIR", unset = file.path(PROJECT_DIR, "data", "gadm"))

EXTERNAL_PATHS <- c(
  "SELECTION_VCF", "SELECTION_ANNOTATION",
  "POPGEN_VCF", "MALARIA_CASES_FILE"
)

local_config <- file.path(PROJECT_DIR, "config.local.R")
if (file.exists(local_config)) source(local_config)

unset_paths <- Filter(function(n) identical(get(n), UNSET), EXTERNAL_PATHS)
if (length(unset_paths) > 0) {
  message(
    "config.R: no path set for ", paste(unset_paths, collapse = ", "),
    "\n  Copy config.local.R.example to config.local.R and fill it in.",
    "\n  Scripts that do not need these paths will still run."
  )
}
