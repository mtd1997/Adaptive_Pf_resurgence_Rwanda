# Figure 4d — neighbour-joining tree of the Kigali samples
#
# Inputs:  POPGEN_VCF (config.R); the metadata roster in data/metadata/
#          --reference  the saved Newick the published panel was drawn from
#
# Outputs: results/figures/figure04/figure04d_kigali_tree.nwk
#          results/tables/figure04/figure04d_tip_metadata.tsv
#          results/tables/figure04/figure04d_itol_{period,k13,district}.txt
#
# Usage:   Rscript scripts/06_figures/figure04d_distance_tree.R \
#            --reference /path/to/..._Kigali_tree.nwk
#
# The published panel was rendered in iTOL from the Newick file and the
# annotation files below; no script draws it. Two scripts in the original
# analysis built a Kigali tree with the same genotype matrix but different
# distances, so both are built here and the saved tree decides which was
# used.


# ---- setup ------------------------------------------------------------------

options(error = NULL)

required_packages <- c(
  "dplyr", "readr", "stringr", "tibble", "vcfR", "adegenet", "ape"
)

absent_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(absent_packages) > 0) {
  stop(
    "These packages are needed but not installed:\n  ",
    paste(absent_packages, collapse = ", "),
    "\n\nInstall them with:\n  install.packages(c(",
    paste0("\"", absent_packages, "\"", collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(vcfR)
  library(adegenet)
  library(ape)
})


find_project_root <- function() {

  candidate <- normalizePath(getwd(), mustWork = TRUE)

  for (level in seq_len(5)) {
    if (file.exists(file.path(candidate, "config.R"))) return(candidate)
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }

  stop("Could not find config.R. Run this from the top of the repository.",
       call. = FALSE)
}

PROJECT_ROOT <- find_project_root()

source(file.path(PROJECT_ROOT, "config.R"))


arguments <- commandArgs(trailingOnly = TRUE)

get_argument <- function(flag, default = NULL) {
  position <- match(flag, arguments)
  if (is.na(position)) return(default)
  if (position == length(arguments)) {
    stop("The ", flag, " flag needs a value after it.", call. = FALSE)
  }
  arguments[[position + 1]]
}

VCF_FILE <- get_argument(
  "--vcf",
  if (exists("POPGEN_VCF", inherits = TRUE)) {
    get("POPGEN_VCF", inherits = TRUE)
  } else NULL
)

if (is.null(VCF_FILE) || !file.exists(VCF_FILE)) {
  stop("VCF not found. Set POPGEN_VCF in config.R or pass --vcf PATH.\n",
       "  tried: ", if (is.null(VCF_FILE)) "(unset)" else VCF_FILE,
       call. = FALSE)
}

METADATA_DIR <- file.path(PROJECT_ROOT, "data", "metadata")

METADATA_FILE <- get_argument(
  "--metadata",
  {
    preferred <- file.path(METADATA_DIR, "samples_metadata_figure04.tsv")
    fallback  <- file.path(METADATA_DIR, "samples_metadata.tsv")
    if (file.exists(preferred)) preferred else fallback
  }
)

if (!file.exists(METADATA_FILE)) {
  stop("Metadata not found:\n  ", METADATA_FILE, call. = FALSE)
}

REFERENCE_TREE <- get_argument("--reference", NULL)
FORCED_METHOD  <- get_argument("--method", NULL)

if (!is.null(FORCED_METHOD) &&
    !FORCED_METHOD %in% c("euclidean", "manhattan")) {
  stop("--method must be euclidean or manhattan.", call. = FALSE)
}

FIGURE_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figure04")
TABLE_DIR  <- file.path(PROJECT_ROOT, "results", "tables",  "figure04")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)


KIGALI_DISTRICTS <- c("Gasabo", "Kicukiro", "Nyarugenge")

REFERENCE_SAMPLE_NAMES <- c("3D7", "Pf3D7", "PF3D7")

# Relative, because write.tree() rounds branch lengths on save.
TOLERANCE <- 1e-6

# Labels and colours as they appear on the published panel.
PERIOD_DISPLAY <- tibble::tribble(
  ~code,              ~label,             ~colour,
  "Pre-Resurgence",   "Baseline-Kigali",  "#52C1CB",
  "Baseline-Kigali",  "Baseline-Kigali",  "#52C1CB",
  "Post-Resurgence",  "Resurgence",       "#FF757D",
  "Resurgence",       "Resurgence",       "#FF757D"
)

K13_DISPLAY <- tibble::tribble(
  ~code,          ~label,           ~colour,
  "WT",           "Wild Type",      "#51E19D",
  "R561H_only",   "R561H only",     "#FFC033",
  "K189T_only",   "K189T only",     "#E35199",
  "R561H+K189T",  "R561H + K189T",  "#FE8133",
  "Other_K13",    "Other K13",      "#937EF3"
)

DISTRICT_DISPLAY <- tibble::tribble(
  ~code,         ~label,        ~colour,
  "Gasabo",      "Gasabo",      "#1B9E77",
  "Kicukiro",    "Kicukiro",    "#D95F02",
  "Nyarugenge",  "Nyarugenge",  "#7570B3"
)

UNASSIGNED_COLOUR <- "#BDBDBD"


report_lines <- character(0)

add_to_report <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  message(line)
}


# ---- read -------------------------------------------------------------------

add_to_report(strrep("=", 72))
add_to_report("FIGURE 4D — KIGALI NEIGHBOUR-JOINING TREE")
add_to_report(strrep("=", 72))
add_to_report("VCF:       ", VCF_FILE)
add_to_report("Metadata:  ", METADATA_FILE)
add_to_report("Reference: ",
              if (is.null(REFERENCE_TREE)) "(none given)" else REFERENCE_TREE)
add_to_report("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
add_to_report("")

metadata <- readr::read_tsv(
  METADATA_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)

# The roster uses lower-case names; the original scripts used SAMPLE /
# District / K13_status. Accept either.
normalise_column <- function(data, candidates, target) {

  present <- intersect(candidates, names(data))

  if (length(present) == 0) {
    stop("The metadata has none of these columns: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }

  data[[target]] <- stringr::str_trim(as.character(data[[present[[1]]]]))
  data
}

metadata <- metadata %>%
  normalise_column(c("SAMPLE", "sample"), "sample_id") %>%
  normalise_column(c("District", "district"), "district_name") %>%
  normalise_column(c("PERIOD", "period"), "period_code") %>%
  normalise_column(c("K13_status", "k13_status"), "k13_code")

if (any(c("REGION", "region") %in% names(metadata))) {
  metadata <- normalise_column(metadata, c("REGION", "region"), "region_name")
} else {
  metadata$region_name <- NA_character_
}

add_to_report(sprintf("Metadata rows: %d", nrow(metadata)))

vcf <- vcfR::read.vcfR(VCF_FILE, verbose = FALSE)

genlight_all <- vcfR::vcfR2genlight(vcf)

genlight_all <- genlight_all[
  !adegenet::indNames(genlight_all) %in% REFERENCE_SAMPLE_NAMES,
]

vcf_samples <- adegenet::indNames(genlight_all)

add_to_report(sprintf(
  "VCF samples after dropping the reference: %d", length(vcf_samples)
))


# ---- the Kigali subset ------------------------------------------------------

kigali_metadata <- metadata %>%
  dplyr::filter(
    sample_id %in% vcf_samples,
    district_name %in% KIGALI_DISTRICTS
  ) %>%
  dplyr::arrange(match(sample_id, vcf_samples))

if (nrow(kigali_metadata) < 3) {
  stop(
    "The Kigali subset has fewer than three samples.\n",
    "Districts present in the metadata: ",
    paste(sort(unique(metadata$district_name)), collapse = ", "),
    call. = FALSE
  )
}

add_to_report("")
add_to_report(sprintf("Kigali by district: %d", nrow(kigali_metadata)))

if (!all(is.na(metadata$region_name))) {
  by_region <- metadata %>%
    dplyr::filter(sample_id %in% vcf_samples, region_name == "Kigali")
  add_to_report(sprintf("Kigali by region:   %d  (differs by %d sample(s))",
                        nrow(by_region),
                        length(union(
                          setdiff(kigali_metadata$sample_id, by_region$sample_id),
                          setdiff(by_region$sample_id, kigali_metadata$sample_id)
                        ))))
}


# ---- build both candidates --------------------------------------------------
#
# NA.method = "mean" fills a missing genotype with its column mean. The
# Figure 4 PCoA panels instead set missing genotypes to zero, so the two
# matrices are not interchangeable.

build_candidates <- function(sample_ids) {

  subset_genlight <- genlight_all[
    adegenet::indNames(genlight_all) %in% sample_ids,
  ]

  subset_genlight <- subset_genlight[
    match(sample_ids, adegenet::indNames(subset_genlight)),
  ]

  stopifnot(identical(adegenet::indNames(subset_genlight), sample_ids))

  genotype_matrix <- adegenet::tab(subset_genlight, NA.method = "mean")

  list(
    matrix = genotype_matrix,
    trees = list(
      euclidean = ape::nj(stats::dist(genotype_matrix)),
      manhattan = ape::nj(stats::dist(genotype_matrix, method = "manhattan"))
    )
  )
}

built <- build_candidates(kigali_metadata$sample_id)

add_to_report(sprintf("Genotype matrix: %d samples x %d sites",
                      nrow(built$matrix), ncol(built$matrix)))

candidates <- built$trees


# ---- which candidate is the published tree ----------------------------------

compare_to_reference <- function(candidate, reference) {

  shared <- intersect(candidate$tip.label, reference$tip.label)

  if (length(shared) != length(candidate$tip.label) ||
      length(shared) != length(reference$tip.label)) {
    return(list(comparable = FALSE))
  }

  candidate_distances <- stats::cophenetic(candidate)[shared, shared]
  reference_distances <- stats::cophenetic(reference)[shared, shared]

  largest <- max(abs(candidate_distances - reference_distances))

  list(
    comparable = TRUE,
    largest_difference = largest,
    relative_difference = largest / max(reference_distances),
    correlation = stats::cor(as.vector(candidate_distances),
                             as.vector(reference_distances))
  )
}

chosen_method <- NULL

if (!is.null(REFERENCE_TREE)) {

  if (!file.exists(REFERENCE_TREE)) {
    stop("Reference tree not found:\n  ", REFERENCE_TREE, call. = FALSE)
  }

  reference <- ape::read.tree(REFERENCE_TREE)

  add_to_report("")
  add_to_report(sprintf("Saved tree: %d tips", length(reference$tip.label)))

  only_in_saved <- setdiff(reference$tip.label, kigali_metadata$sample_id)
  only_in_built <- setdiff(kigali_metadata$sample_id, reference$tip.label)

  if (length(only_in_saved) > 0) {
    stop(
      sprintf("The saved tree contains %d sample(s) absent from this subset:\n  %s\n\n",
              length(only_in_saved),
              paste(utils::head(only_in_saved, 20), collapse = ", ")),
      "The roster or the VCF is not the one the saved tree was built from.",
      call. = FALSE
    )
  }

  # A neighbour-joining topology depends on every tip, so the candidates are
  # rebuilt on exactly the saved tree's samples rather than pruned.
  if (length(only_in_built) > 0) {

    add_to_report(sprintf(
      "This subset has %d sample(s) the saved tree does not: %s",
      length(only_in_built), paste(only_in_built, collapse = ", ")
    ))
    add_to_report("Rebuilding on the saved tree's samples for comparison.")

    readr::write_tsv(
      tibble::tibble(sample = only_in_built),
      file.path(TABLE_DIR, "figure04d_samples_not_in_saved_tree.tsv")
    )

    kigali_metadata <- kigali_metadata %>%
      dplyr::filter(sample_id %in% reference$tip.label) %>%
      dplyr::arrange(match(sample_id, vcf_samples))

    built <- build_candidates(kigali_metadata$sample_id)
    candidates <- built$trees

    add_to_report(sprintf("Rebuilt on %d samples (%d sites).",
                          nrow(built$matrix), ncol(built$matrix)))
  }

  add_to_report("")

  for (method in names(candidates)) {

    result <- compare_to_reference(candidates[[method]], reference)

    if (!isTRUE(result$comparable)) {
      add_to_report(sprintf("  %-10s tip sets differ", method))
      next
    }

    matched <- result$relative_difference < TOLERANCE

    add_to_report(sprintf(
      "  %-10s largest difference %.4g (relative %.3e)  correlation %.6f  %s",
      method, result$largest_difference, result$relative_difference,
      result$correlation, if (matched) "MATCH" else "no"
    ))

    if (matched) chosen_method <- method
  }

  if (is.null(chosen_method)) {
    stop(
      "Neither candidate reproduces the saved tree. Nothing has been written.",
      call. = FALSE
    )
  }

  add_to_report("")
  add_to_report(sprintf("Using the %s tree.", chosen_method))

} else if (!is.null(FORCED_METHOD)) {

  chosen_method <- FORCED_METHOD
  add_to_report("")
  add_to_report(sprintf(
    "No reference given; building the %s tree because --method said so.",
    chosen_method
  ))

} else {

  stop(
    "Two distances could have built the published tree and without the\n",
    "saved tree there is no way to tell which.\n\n",
    "  --reference /path/to/..._Kigali_tree.nwk\n",
    "or force one with --method euclidean | manhattan.",
    call. = FALSE
  )
}

tree <- candidates[[chosen_method]]


# ---- tip metadata and iTOL annotations --------------------------------------

lookup <- function(codes, table) {
  index <- match(codes, table$code)
  list(
    label  = ifelse(is.na(index), codes, table$label[index]),
    colour = ifelse(is.na(index), UNASSIGNED_COLOUR, table$colour[index])
  )
}

tip_metadata <- tibble::tibble(sample_id = tree$tip.label) %>%
  dplyr::left_join(
    kigali_metadata %>%
      dplyr::select(sample_id, period_code, district_name, k13_code),
    by = "sample_id"
  )

if (any(is.na(tip_metadata$period_code))) {
  stop(sprintf("%d tip(s) have no metadata.",
               sum(is.na(tip_metadata$period_code))), call. = FALSE)
}

period   <- lookup(tip_metadata$period_code,   PERIOD_DISPLAY)
k13      <- lookup(tip_metadata$k13_code,      K13_DISPLAY)
district <- lookup(tip_metadata$district_name, DISTRICT_DISPLAY)

tip_metadata$period_label   <- period$label
tip_metadata$k13_label      <- k13$label
tip_metadata$district_label <- district$label

unmapped <- unique(c(
  tip_metadata$period_code[!tip_metadata$period_code %in% PERIOD_DISPLAY$code],
  tip_metadata$k13_code[!tip_metadata$k13_code %in% K13_DISPLAY$code]
))

if (length(unmapped) > 0) {
  add_to_report("")
  add_to_report("Values with no colour assigned, drawn grey: ",
                paste(unmapped, collapse = ", "))
}

write_itol_colorstrip <- function(colours, display_table, dataset_label,
                                  filename) {

  # PERIOD_DISPLAY maps two codes onto each label, so collapse first.
  used <- dplyr::distinct(display_table, label, colour)

  writeLines(
    c(
      "DATASET_COLORSTRIP",
      "SEPARATOR TAB",
      paste0("DATASET_LABEL\t", dataset_label),
      "COLOR\t#000000",
      "STRIP_WIDTH\t25",
      "MARGIN\t5",
      "SHOW_INTERNAL\t0",
      paste0("LEGEND_TITLE\t", dataset_label),
      paste0("LEGEND_SHAPES\t", paste(rep(1, nrow(used)), collapse = "\t")),
      paste0("LEGEND_COLORS\t", paste(used$colour, collapse = "\t")),
      paste0("LEGEND_LABELS\t", paste(used$label, collapse = "\t")),
      "DATA",
      paste(tip_metadata$sample_id, colours, sep = "\t")
    ),
    con = filename
  )
}

write_itol_colorstrip(
  period$colour, PERIOD_DISPLAY, "Period",
  file.path(TABLE_DIR, "figure04d_itol_period.txt")
)

write_itol_colorstrip(
  k13$colour, K13_DISPLAY, "K13 genotype",
  file.path(TABLE_DIR, "figure04d_itol_k13.txt")
)

write_itol_colorstrip(
  district$colour, DISTRICT_DISPLAY, "District",
  file.path(TABLE_DIR, "figure04d_itol_district.txt")
)


# ---- write ------------------------------------------------------------------

ape::write.tree(tree, file = file.path(FIGURE_DIR, "figure04d_kigali_tree.nwk"))

readr::write_tsv(
  tip_metadata %>%
    dplyr::select(
      sample   = sample_id,
      period   = period_label,
      district = district_label,
      k13      = k13_label
    ),
  file.path(TABLE_DIR, "figure04d_tip_metadata.tsv")
)

add_to_report("")
add_to_report(strrep("=", 72))
add_to_report("WRITTEN")
add_to_report(strrep("=", 72))
add_to_report("  figures/figure04/figure04d_kigali_tree.nwk")
add_to_report("  tables/figure04/figure04d_tip_metadata.tsv")
add_to_report("  tables/figure04/figure04d_itol_period.txt")
add_to_report("  tables/figure04/figure04d_itol_k13.txt")
add_to_report("  tables/figure04/figure04d_itol_district.txt")
add_to_report("")
add_to_report("Upload the .nwk to iTOL, then drag on the annotation files.")
add_to_report("The published panel has Period as the inner ring and")
add_to_report("K13 genotype as the outer one.")

writeLines(report_lines, file.path(TABLE_DIR, "figure04d_report.txt"))

capture.output(
  sessionInfo(),
  file = file.path(TABLE_DIR, "figure04d_sessionInfo.txt")
)
