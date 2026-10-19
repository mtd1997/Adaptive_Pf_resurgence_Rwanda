# rename_fastq.R — rename merged FASTQ files from barcode names to sample IDs
#
# Inputs:  BARCODE_MAP   CSV with two columns:
#          SampleID  study sample identifier, e.g. A25DBL004
#          barcode   merged file name, e.g. barcode01_run1.fastq.gz
#          One row per sample, covering every sequencing run.
# Outputs: The same folder, with files renamed to <SampleID>.fastq.gz
#
# Usage:   From the repository root:
#          Rscript scripts/01_preprocessing/rename_fastq.R

if (!file.exists("config.R")) {
  stop("config.R not found. Run this script from the repository root:\n",
       "  Rscript scripts/01_preprocessing/rename_fastq.R", call. = FALSE)
}
source("config.R")

# --- dry-run flag ---------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

if (dry_run) {
  message("DRY RUN: showing what would happen, renaming nothing.\n")
}

# --- read the mapping table -------------------------------------------------------
if (!file.exists(BARCODE_MAP)) {
  stop("Mapping table not found: ", BARCODE_MAP, call. = FALSE)
}

# read.csv with sep = "" lets R detect comma or tab automatically
mapping <- read.csv(BARCODE_MAP, stringsAsFactors = FALSE,
                    sep = if (grepl("\t", readLines(BARCODE_MAP, n = 1))) "\t" else ",")

required_cols <- c("SampleID", "barcode")
missing_cols  <- setdiff(required_cols, names(mapping))
if (length(missing_cols) > 0) {
  stop("Mapping table is missing column(s): ", paste(missing_cols, collapse = ", "),
       "\nFound: ", paste(names(mapping), collapse = ", "), call. = FALSE)
}

mapping$SampleID <- trimws(mapping$SampleID)
mapping$barcode  <- trimws(mapping$barcode)

message("Mapping table : ", BARCODE_MAP)
message("Rows          : ", nrow(mapping))
message("FASTQ folder  : ", MERGED_DIR)

# --- validate the table before renaming anything ---------------------------------
dup_samples  <- mapping$SampleID[duplicated(mapping$SampleID)]
dup_barcodes <- mapping$barcode[duplicated(mapping$barcode)]

if (length(dup_samples) > 0) {
  stop("Duplicate SampleID values in the mapping table: ",
       paste(unique(dup_samples), collapse = ", "),
       "\nEach sample must appear once, or renaming would overwrite data.",
       call. = FALSE)
}
if (length(dup_barcodes) > 0) {
  stop("Duplicate barcode values in the mapping table: ",
       paste(unique(dup_barcodes), collapse = ", "),
       "\nThis usually means a run identifier is missing from the barcode names.",
       call. = FALSE)
}

# --- compare the table against what is actually on disk ----------------------------
fastq_present <- list.files(MERGED_DIR, pattern = "\\.fastq\\.gz$")

expected <- mapping$barcode
found    <- intersect(expected, fastq_present)
absent   <- setdiff(expected, fastq_present)
unlisted <- setdiff(fastq_present, c(expected, paste0(mapping$SampleID, ".fastq.gz")))

message("\n--- before renaming ---")
message("Files in folder        : ", length(fastq_present))
message("Listed and present     : ", length(found))
message("Listed but not on disk : ", length(absent))
message("On disk but not listed : ", length(unlisted))

if (length(absent) > 0) {
  message("\nListed in the table but missing from the folder:")
  message(paste("   ", head(absent, 20), collapse = "\n"))
  if (length(absent) > 20) message("    ... and ", length(absent) - 20, " more")
}
if (length(unlisted) > 0) {
  message("\nIn the folder but absent from the table (will be left alone):")
  message(paste("   ", head(unlisted, 20), collapse = "\n"))
  if (length(unlisted) > 20) message("    ... and ", length(unlisted) - 20, " more")
}

if (length(found) == 0) {
  stop("\nNo files to rename. Check that MERGED_DIR in config.R points at the ",
       "output of concatenate_fastq.sh, and that the barcode column matches ",
       "the merged file names.", call. = FALSE)
}

# --- refuse to overwrite an existing destination -------------------------------------
to_rename   <- mapping[mapping$barcode %in% found, ]
destinations <- file.path(MERGED_DIR, paste0(to_rename$SampleID, ".fastq.gz"))
clashes      <- destinations[file.exists(destinations)]

if (length(clashes) > 0) {
  stop("\nThese destination files already exist:\n",
       paste("   ", basename(clashes), collapse = "\n"),
       "\nRenaming would overwrite them. Resolve this before continuing.",
       call. = FALSE)
}

# --- rename -----------------------------------------------------------------------------
message("\n--- renaming ", nrow(to_rename), " files ---")

renamed <- 0
failed  <- character(0)

for (i in seq_len(nrow(to_rename))) {

  old_path <- file.path(MERGED_DIR, to_rename$barcode[i])
  new_path <- file.path(MERGED_DIR, paste0(to_rename$SampleID[i], ".fastq.gz"))

  if (dry_run) {
    message("  would rename: ", to_rename$barcode[i], " -> ",
            to_rename$SampleID[i], ".fastq.gz")
    renamed <- renamed + 1
    next
  }

  if (file.rename(old_path, new_path)) {
    message("  ", to_rename$barcode[i], " -> ", to_rename$SampleID[i], ".fastq.gz")
    renamed <- renamed + 1
  } else {
    failed <- c(failed, to_rename$barcode[i])
    message("  FAILED: ", to_rename$barcode[i])
  }
}

# --- summary -------------------------------------------------------------------------------
message("\n============================================================")
if (dry_run) {
  message("DRY RUN complete. ", renamed, " files would be renamed.")
  message("Re-run without --dry-run to apply.")
} else {
  message("Renamed : ", renamed)
  if (length(failed) > 0) {
    message("Failed  : ", length(failed), " -> ", paste(failed, collapse = ", "))
  }
  message("Next    : bash scripts/01_preprocessing/nanoplot_qc.sh")
}
message("============================================================")
