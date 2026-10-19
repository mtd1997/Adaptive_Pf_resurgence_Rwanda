# Supplementary Figure S1 — assembly
#
# Inputs:  results/figures/panels/figureS1a_baseline_map.rds
#          results/figures/panels/figureS1b_malaria_cases.rds
# Outputs: results/figures/figureS1/figureS1_horizontal.svg / .pdf / .png
#          results/figures/figureS1/figureS1_vertical.svg / .pdf / .png
#
# Usage:   From the top of the repository, after the two upstream scripts:

options(error = NULL)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})


find_project_root <- function() {

  candidate <- normalizePath(getwd(), mustWork = TRUE)

  for (level in seq_len(5)) {

    if (file.exists(file.path(candidate, "config.R"))) {
      return(candidate)
    }

    parent <- dirname(candidate)

    if (identical(parent, candidate)) {
      break
    }

    candidate <- parent
  }

  stop(
    "Could not find config.R. Run this from the top of the repository.",
    call. = FALSE
  )
}

PROJECT_ROOT <- find_project_root()

source(file.path(PROJECT_ROOT, "config.R"))

PANEL_DIR    <- file.path(PROJECT_ROOT, "results", "figures", "panels")
FIGURES1_DIR <- file.path(PROJECT_ROOT, "results", "figures", "figureS1")

dir.create(FIGURES1_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- read the panels ---------------------------------------------------

read_panel <- function(filename, produced_by) {

  path <- file.path(PANEL_DIR, filename)

  if (!file.exists(path)) {
    stop(
      "Panel file not found:\n  ", path,
      "\n\nRun this first:\n  Rscript scripts/06_figures/", produced_by,
      call. = FALSE
    )
  }

  readRDS(path)
}

panel_a_data <- read_panel(
  "figureS1a_baseline_map.rds",
  "figure01_sampling_maps.R"
)

panel_b_data <- read_panel(
  "figureS1b_malaria_cases.rds",
  "figureS1b_malaria_cases.R"
)

message("Panels read:")
message(sprintf(
  "  A  baseline sampling map   (%d samples)", panel_a_data$n
))
message("  B  monthly malaria cases")


# ---- panel a -----------------------------------------------------------

panel_a <- panel_a_data$panel +

  ggplot2::scale_size_continuous(
    name   = "Number of samples",
    range  = c(3, 9),
    breaks = panel_a_data$breaks
  ) +

  ggplot2::labs(tag = "A") +

  ggplot2::guides(
    size = ggplot2::guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  ) +

  ggplot2::theme(
    plot.tag          = ggplot2::element_text(
      size = 24, face = "bold", colour = "black"
    ),
    plot.tag.position = c(0.01, 0.99),

    legend.title      = ggplot2::element_text(size = 12, face = "bold"),
    legend.text       = ggplot2::element_text(size = 11),
    legend.position   = c(0.72, 0.10),
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.background = ggplot2::element_rect(fill = "white", colour = NA),
    legend.margin     = ggplot2::margin(0, 0, 0, 0),

    plot.margin       = ggplot2::margin(5, 5, 5, 5)
  )


# ---- panel b -----------------------------------------------------------

panel_b <- panel_b_data$panel +

  ggplot2::labs(tag = "B") +

  ggplot2::theme(
    plot.tag          = ggplot2::element_text(
      size = 24, face = "bold", colour = "black"
    ),
    plot.tag.position = c(0.01, 0.99),

    axis.title        = ggplot2::element_text(size = 14, face = "bold"),
    axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1, size = 9.5),
    axis.text.y       = ggplot2::element_text(size = 10.5),

    strip.text        = ggplot2::element_text(size = 13, face = "bold"),

    legend.title      = ggplot2::element_text(size = 12, face = "bold"),
    legend.text       = ggplot2::element_text(size = 11),
    legend.position   = "bottom",

    plot.margin       = ggplot2::margin(10, 10, 10, 10)
  )


# ---- the two layouts ---------------------------------------------------

figure_horizontal <- (panel_a | panel_b) +
  patchwork::plot_layout(widths = c(1.15, 1.55)) &
  ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10))

figure_vertical <- (panel_a / panel_b) +
  patchwork::plot_layout(heights = c(1.0, 1.30)) &
  ggplot2::theme(plot.margin = ggplot2::margin(8, 10, 8, 10))


# ---- save --------------------------------------------------------------

save_layout <- function(plot_object, basename, width, height) {

  for (extension in c("svg", "pdf", "png")) {

    ggplot2::ggsave(
      filename  = file.path(
        FIGURES1_DIR, paste0(basename, ".", extension)
      ),
      plot      = plot_object,
      width     = width,
      height    = height,
      units     = "in",
      dpi       = 600,
      bg        = "white",
      limitsize = FALSE
    )
  }
}

save_layout(figure_horizontal, "figureS1_horizontal", 18, 8)
save_layout(figure_vertical,   "figureS1_vertical",   12, 15)

capture.output(
  sessionInfo(),
  file = file.path(FIGURES1_DIR, "figureS1_sessionInfo.txt")
)


# ---- verify and summarise ----------------------------------------------

expected_files <- c(
  outer(
    c("figureS1_horizontal", "figureS1_vertical"),
    c("svg", "pdf", "png"),
    function(base, extension) {
      file.path(FIGURES1_DIR, paste0(base, ".", extension))
    }
  )
)

file_status <- data.frame(
  file   = basename(expected_files),
  exists = file.exists(expected_files)
)

message("\n", strrep("=", 72))
message("SUPPLEMENTARY FIGURE S1 ASSEMBLED")
message(strrep("=", 72))

print(file_status, row.names = FALSE)

if (!all(file_status$exists)) {
  warning(
    "Some output files were not written. Check the messages above.",
    call. = FALSE
  )
}

message(sprintf("\nOutput directory:\n  %s", FIGURES1_DIR))
message("\nTwo layouts were produced. Pick whichever suits the column width:")
message("  figureS1_horizontal   A beside B,  18 x 8 in")
message("  figureS1_vertical     A above B,   12 x 15 in")
message(strrep("=", 72))
