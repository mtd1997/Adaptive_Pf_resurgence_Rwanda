################################################################################
# SUPPLEMENTARY FIGURE S1 — ASSEMBLY
################################################################################
#
# WHAT IT DOES
#   Combines two panels that were drawn by other scripts:
#
#     A   Baseline-Kigali sampling map, 2021-2023
#         from figure01_sampling_maps.R
#
#     B   Monthly reported malaria cases in four Kigali facilities,
#         2023 versus 2024
#         from figureS1b_malaria_cases.R
#
#   It produces both a horizontal layout (A beside B) and a vertical one
#   (A above B), because which fits better depends on the journal's column
#   width.
#
# WHY IT READS FILES INSTEAD OF USING THE SESSION
#   The original version expected three objects -- `map_baseline`,
#   `baseline_breaks` and `p` -- to already exist in the R session, and
#   stopped with an error if they did not. That required three scripts to be
#   run in one session in the right order, and meant the figure could not be
#   rebuilt on its own later.
#
#   The two upstream scripts now save their panels to
#   results/figures/panels/, and this one reads them. Same figure, but each
#   script can be run independently and in any order.
#
# INPUTS
#   results/figures/panels/figureS1a_baseline_map.rds
#   results/figures/panels/figureS1b_malaria_cases.rds
#
# OUTPUTS
#   results/figures/figureS1/figureS1_horizontal.svg / .pdf / .png
#   results/figures/figureS1/figureS1_vertical.svg / .pdf / .png
#
#   The individual panels are written by their own scripts, so they are not
#   duplicated here.
#
# HOW TO RUN
#   From the top of the repository, after the two upstream scripts:
#
#     Rscript scripts/06_figures/figure01_sampling_maps.R
#     Rscript scripts/06_figures/figureS1b_malaria_cases.R
#     Rscript scripts/06_figures/figureS1_assemble.R
#
################################################################################


# ==============================================================================
# 1. SETUP
# ==============================================================================

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


# ==============================================================================
# 2. READ THE PANELS
# ==============================================================================

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


# ==============================================================================
# 3. PANEL A
# ==============================================================================
#
# The size scale is restated so the legend can sit inside the map, in the
# lower right, where there is empty space. The breaks travel with the panel
# from the upstream script rather than being recomputed, so they cannot
# differ from the ones used when the map was drawn.
# ==============================================================================

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


# ==============================================================================
# 4. PANEL B
# ==============================================================================

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


# ==============================================================================
# 5. THE TWO LAYOUTS
# ==============================================================================
#
# Panel B is given the greater share in both layouts because it carries four
# facets while panel A is a single map.
# ==============================================================================

figure_horizontal <- (panel_a | panel_b) +
  patchwork::plot_layout(widths = c(1.15, 1.55)) &
  ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10))

figure_vertical <- (panel_a / panel_b) +
  patchwork::plot_layout(heights = c(1.0, 1.30)) &
  ggplot2::theme(plot.margin = ggplot2::margin(8, 10, 8, 10))


# ==============================================================================
# 6. SAVE
# ==============================================================================
#
# Nothing is printed to the graphics device: a large sf panel can fail when
# the plot pane is small, and printing under Rscript leaves a stray
# Rplots.pdf behind.
# ==============================================================================

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


# ==============================================================================
# 7. VERIFY AND SUMMARISE
# ==============================================================================

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
