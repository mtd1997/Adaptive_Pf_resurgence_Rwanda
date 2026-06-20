
# Statistical analyses for K13 mutations in Rwanda

# This script summarizes the statistical analysis used to compare K13 mutation profiles.
# The main focus is the comparison between wild type parasites and the R561H + K189T double mutant.
# The script calculates proportions, runs exact binomial tests, compares groups with Fisher tests,
# adds 95% confidence intervals, and uses a simple logistic regression to confirm the same trends.
# The plots are generated to make the differences between periods and locations easier to visualize.


# 1. Packages and paths


packages <- c(
  "tidyverse",
  "readxl",
  "data.table",
  "stringr",
  "purrr",
  "dplyr",
  "ggtext",
  "tidyr"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

input_file <- "/Users/diallo/Desktop/Rwanda_sequenced_data_24-12-2025.xlsx"
output_dir <- "/Users/diallo/Desktop"

df <- read_excel(input_file)


# 2. Small helper functions


# Convert a p-value into the usual star notation used on the plots.
p_to_stars <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}

# Basic check to make sure the counts match the expected total.
check_counts <- function(counts_df, total_n, label) {
  cat("\n============================\n")
  cat("CHECK:", label, "\n")
  cat("============================\n")
  
  print(counts_df %>% arrange(desc(n)))
  
  cat("\nExpected total =", total_n, "\n")
  cat("Observed sum =", sum(counts_df$n), "\n")
  
  stopifnot(sum(counts_df$n) == total_n)
  
  wt <- counts_df %>% filter(category == "Wild type") %>% pull(n)
  mut <- counts_df %>% filter(category == "R561H + K189T") %>% pull(n)
  
  if (length(wt) != 1) stop("Problem with Wild type count")
  if (length(mut) != 1) stop("Problem with R561H + K189T count")
  stopifnot(wt >= 0)
  
  cat("Wild type =", wt, "\n")
  cat("R561H + K189T =", mut, "\n")
}

# This classification keeps the main K13 groups used in the figures.
classify_k13 <- function(muts) {
  muts <- unique(muts)
  
  has_R561H <- "R561H" %in% muts
  has_K189T <- "K189T" %in% muts
  other_muts <- setdiff(muts, c("R561H", "K189T"))
  
  if (has_R561H && has_K189T) return("R561H + K189T")
  if (has_R561H && length(other_muts) == 0) return("R561H only")
  if (has_K189T && length(other_muts) == 0) return("K189T only")
  if (has_R561H && length(other_muts) > 0) return("R561H + other")
  if (has_K189T && length(other_muts) > 0) return("K189T + other")
  
  if ("H384R" %in% muts) return("H384R")
  if ("A675V" %in% muts) return("A675V")
  if ("P574L" %in% muts) return("P574L")
  if ("C469Y" %in% muts) return("C469Y")
  if ("P441L" %in% muts) return("P441L")
  if ("G449A" %in% muts) return("G449A")
  if ("C496F" %in% muts) return("C496F")
  if ("N537D" %in% muts) return("N537D")
  if ("R515K" %in% muts) return("R515K")
  
  return("Other (non-WHO)")
}

# Build one mutation category per sample.
make_k13_categories <- function(data) {
  data %>%
    filter(GENE_NAME == "K13") %>%
    group_by(SAMPLE) %>%
    summarise(
      K13_mutations = paste(unique(HGVS_PROTEIN), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(
      category = map_chr(
        K13_mutations,
        ~ classify_k13(str_split(.x, ",\\s*")[[1]])
      )
    )
}

# Add samples with no K13 mutation as wild type.
add_wild_type <- function(k13_cat_df, total_samples) {
  counts <- k13_cat_df %>%
    count(category, name = "n")
  
  wt_n <- total_samples - sum(counts$n)
  
  bind_rows(
    counts,
    tibble(category = "Wild type", n = wt_n)
  )
}

# Safe extraction for one category count.
get_count <- function(counts_df, category_name) {
  value <- counts_df %>%
    filter(category == category_name) %>%
    pull(n)
  
  if (length(value) == 0) 0 else value
}

# Run the exact binomial test between wild type and double mutant.
run_wt_double_binom <- function(counts_df) {
  wild_n <- get_count(counts_df, "Wild type")
  double_n <- get_count(counts_df, "R561H + K189T")
  
  test <- binom.test(
    x = wild_n,
    n = wild_n + double_n,
    p = 0.5,
    alternative = "two.sided"
  )
  
  list(
    wild_n = wild_n,
    double_n = double_n,
    test = test,
    wild_pct = round(100 * wild_n / (wild_n + double_n), 1),
    double_pct = round(100 * double_n / (wild_n + double_n), 1),
    p_value = signif(test$p.value, 3)
  )
}

# Print the exact CI for wild type and the complementary CI for the double mutant.
print_binom_ci <- function(test_obj, wild_count, mutant_count, label) {
  total_tested <- wild_count + mutant_count
  wild_prop <- wild_count / total_tested
  wild_ci_low <- test_obj$conf.int[1]
  wild_ci_high <- test_obj$conf.int[2]
  
  mutant_prop <- mutant_count / total_tested
  mutant_ci_low <- 1 - wild_ci_high
  mutant_ci_high <- 1 - wild_ci_low
  
  cat("\n====================================\n")
  cat(label, "\n")
  cat("====================================\n")
  cat("p-value =", test_obj$p.value, "\n")
  
  cat("\nWild type proportion =", round(wild_prop * 100, 2), "%\n")
  cat(
    "95% CI Wild type =",
    round(wild_ci_low * 100, 2), "% -",
    round(wild_ci_high * 100, 2), "%\n"
  )
  
  cat("\nR561H + K189T proportion =", round(mutant_prop * 100, 2), "%\n")
  cat(
    "95% CI R561H + K189T =",
    round(mutant_ci_low * 100, 2), "% -",
    round(mutant_ci_high * 100, 2), "%\n"
  )
}

# Common order and colors for the K13 plots.
category_order <- c(
  "Wild type",
  "H384R",
  "A675V",
  "P441L",
  "P574L",
  "R561H + K189T",
  "R561H only",
  "R561H + other",
  "K189T only",
  "K189T + other",
  "C469Y",
  "G449A",
  "C496F",
  "N537D",
  "R515K",
  "Other (non-WHO)"
)

k13_colors <- c(
  "Wild type" = "#2ecc71",
  "H384R" = "#4C6EF5",
  "A675V" = "#D9D9D9",
  "P441L" = "#D8A7FF",
  "P574L" = "#1F5AA6",
  "K189T + other" = "#F4A6C1",
  "K189T only" = "#d62f7d",
  "R561H + K189T" = "#ff7f0e",
  "R561H only" = "#f4b400",
  "R561H + other" = "#C49A00",
  "C469Y" = "#8dd3c7",
  "G449A" = "#fb8072",
  "C496F" = "#80b1d3",
  "N537D" = "#bebada",
  "R515K" = "#b3de69",
  "Other (non-WHO)" = "#5B4CC4"
)


# 3. Total samples and gene-level prevalence by year


total_par_year <- df %>%
  group_by(YEAR) %>%
  summarise(Total_Samples = n_distinct(SAMPLE), .groups = "drop")

write.table(
  total_par_year,
  file = file.path(output_dir, "sample_collected_accross_years.tsv"),
  row.names = FALSE,
  sep = "\t"
)

mut_sum <- df %>%
  group_by(YEAR, GENE_NAME) %>%
  summarise(Mutated_Sample_Count = n_distinct(SAMPLE), .groups = "drop")

mut_sum_percent <- mut_sum %>%
  left_join(total_par_year, by = "YEAR") %>%
  mutate(Percentage = round((Mutated_Sample_Count / Total_Samples) * 100, 2)) %>%
  arrange(YEAR, desc(Mutated_Sample_Count))

mut_sum_percent2 <- mut_sum_percent %>%
  dplyr::select(-Total_Samples) %>%
  left_join(total_par_year, by = "YEAR") %>%
  mutate(FacetLabel = paste0(YEAR, "\nN=", Total_Samples))

write.table(
  mut_sum_percent2,
  file = file.path(output_dir, "figure2.tsv"),
  row.names = FALSE,
  sep = "\t"
)

mut_sum_percent2$FacetLabel <- factor(
  mut_sum_percent2$FacetLabel,
  levels = mut_sum_percent2 %>%
    distinct(YEAR, FacetLabel) %>%
    arrange(YEAR) %>%
    pull(FacetLabel)
)

custom_gene_colors <- c("#B903AA", "#D16100", "#6B7900", "orangered", "#00C2A0")
gene_levels <- unique(mut_sum_percent2$GENE_NAME)
gene_color_mapping <- setNames(custom_gene_colors[seq_along(gene_levels)], gene_levels)

mut_sum_percent2_plot <- ggplot(
  mut_sum_percent2,
  aes(x = GENE_NAME, y = Percentage, fill = GENE_NAME)
) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.3, size = 3, colour = "black") +
  facet_wrap(~ FacetLabel, scales = "free", nrow = 1) +
  scale_fill_manual(values = gene_color_mapping) +
  labs(
    title = "Prevalence of mutated samples per drug resistance gene across years",
    x = "Drug Resistance Genes",
    y = "Percentage of mutated samples"
  ) +
  theme(
    axis.title = element_text(size = 14, colour = "black", face = "bold"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 14, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", size = 20, color = "black", hjust = 0.5),
    strip.text.x = element_text(size = 14, color = "black", face = "bold"),
    strip.background = element_rect(fill = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

print(mut_sum_percent2_plot)


# 4. KFH comparison: 2021-2023 vs 2024-2025 Gasabo


# First comparison kept as in the original analysis:
# 2021-2023 all available samples vs 2024-2025 Gasabo/KFH samples.

df_2123 <- df %>%
  filter(YEAR %in% c(2021, 2022, 2023))

df_2425_kfh <- df %>%
  filter(YEAR %in% c(2024, 2025), District %in% c("Gasabo"))

k13_cat_2123 <- make_k13_categories(df_2123)
k13_cat_2425_kfh <- make_k13_categories(df_2425_kfh)

total_samples_2123 <- df_2123 %>% distinct(SAMPLE) %>% nrow()
total_samples_2425 <- df_2425_kfh %>% distinct(SAMPLE) %>% nrow()

counts_2123 <- add_wild_type(k13_cat_2123, total_samples_2123)
counts_2425 <- add_wild_type(k13_cat_2425_kfh, total_samples_2425)

kfh_2123_test <- run_wt_double_binom(counts_2123)
kfh_2425_test <- run_wt_double_binom(counts_2425)

wild_n <- kfh_2123_test$wild_n
r561h_k189t_n <- kfh_2123_test$double_n
k13_prop_test <- kfh_2123_test$test

wild_n_2425 <- kfh_2425_test$wild_n
r561h_k189t_n_2425 <- kfh_2425_test$double_n
k13_prop_test_2425 <- kfh_2425_test$test

cat(
  "KFH 2021-2023:",
  "Wild type =", kfh_2123_test$wild_pct, "% vs mutant =", kfh_2123_test$double_pct,
  "%, p =", kfh_2123_test$p_value, "\n"
)

cat(
  "KFH 2024-2025:",
  "Wild type =", kfh_2425_test$wild_pct, "% vs mutant =", kfh_2425_test$double_pct,
  "%, p =", kfh_2425_test$p_value, "\n"
)

plot_df_kfh <- bind_rows(
  counts_2123 %>% mutate(period = "2021-2023"),
  counts_2425 %>% mutate(period = "2024-2025")
) %>%
  group_by(period) %>%
  mutate(prevalence = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = category_order)) %>%
  left_join(
    tibble(
      period = c("2021-2023", "2024-2025"),
      total_n = c(total_samples_2123, total_samples_2425)
    ),
    by = "period"
  ) %>%
  mutate(period_label = paste0(period, " KFH (N = ", total_n, ")"))

write.table(
  plot_df_kfh,
  file = file.path(output_dir, "figure4.tsv"),
  row.names = FALSE,
  sep = "\t"
)

sig_kfh <- plot_df_kfh %>%
  group_by(period, period_label) %>%
  summarise(y_max = max(prevalence, na.rm = TRUE), .groups = "drop") %>%
  arrange(period) %>%
  mutate(
    y = y_max * 1.12,
    y_low = y_max * 1.05,
    y_text = y_max * 1.16,
    xstart = 1,
    xend = 6,
    p_value = case_when(
      grepl("^2021", period) ~ k13_prop_test$p.value,
      grepl("^2024", period) ~ k13_prop_test_2425$p.value
    ),
    stars = p_to_stars(p_value),
    label = stars,
    legend_line = paste0(stars, " : p = ", signif(p_value, 3))
  )

K13_evolution_KFH <- ggplot(plot_df_kfh, aes(x = category, y = prevalence, fill = category)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_segment(
    data = sig_kfh,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_segment(
    data = sig_kfh,
    aes(x = xstart, xend = xstart, y = y_low, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_segment(
    data = sig_kfh,
    aes(x = xend, xend = xend, y = y_low, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_text(
    data = sig_kfh,
    aes(x = (xstart + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  geom_text(
    data = sig_kfh,
    aes(x = xend + 0.25, y = y_text + 4, label = legend_line),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(~ period_label, scales = "free", nrow = 1) +
  scale_fill_manual(values = k13_colors) +
  expand_limits(y = max(sig_kfh$y_text, na.rm = TRUE) * 1.15) +
  coord_cartesian(clip = "off") +
  labs(
    x = "K13 mutation category",
    y = "Prevalence (%)"
  ) +
  theme(
    axis.title = element_text(size = 10, colour = "black", face = "bold"),
    axis.text = element_text(size = 7, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

print(K13_evolution_KFH)


# 5. KFH period comparison for wild type and double mutant


period_comparison_df <- bind_rows(
  counts_2123 %>% mutate(period = "2021-2023", total_n = total_samples_2123),
  counts_2425 %>% mutate(period = "2024-2025", total_n = total_samples_2425)
) %>%
  filter(category %in% c("Wild type", "R561H + K189T")) %>%
  mutate(
    prevalence = n / total_n * 100,
    category = factor(category, levels = c("Wild type", "R561H + K189T")),
    period = factor(period, levels = c("2021-2023", "2024-2025"))
  )

wt_2123_n <- get_count(counts_2123, "Wild type")
wt_2425_n <- get_count(counts_2425, "Wild type")
double_2123_n <- get_count(counts_2123, "R561H + K189T")
double_2425_n <- get_count(counts_2425, "R561H + K189T")

test_wt_period <- fisher.test(matrix(
  c(
    wt_2123_n, total_samples_2123 - wt_2123_n,
    wt_2425_n, total_samples_2425 - wt_2425_n
  ),
  nrow = 2,
  byrow = TRUE
))

test_double_period <- fisher.test(matrix(
  c(
    double_2123_n, total_samples_2123 - double_2123_n,
    double_2425_n, total_samples_2425 - double_2425_n
  ),
  nrow = 2,
  byrow = TRUE
))

sig_period <- period_comparison_df %>%
  group_by(category) %>%
  summarise(y_max = max(prevalence, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    x = as.numeric(category),
    xstart = x - 0.22,
    xend = x + 0.22,
    y = y_max * 1.15,
    y_low = y_max * 1.07,
    y_text = y_max * 1.22,
    p_value = case_when(
      category == "Wild type" ~ test_wt_period$p.value,
      category == "R561H + K189T" ~ test_double_period$p.value
    ),
    stars = p_to_stars(p_value),
    legend_line = paste0(category, " : ", stars, " p = ", signif(p_value, 3))
  )

legend_period <- paste(sig_period$legend_line, collapse = "\n")

K13_period_comparison_KFH <- ggplot(
  period_comparison_df,
  aes(x = category, y = prevalence, fill = category, group = period)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "black") +
  scale_fill_manual(values = c("Wild type" = "#2ecc71", "R561H + K189T" = "#ff7f0e")) +
  geom_text(
    aes(label = period),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 3,
    fontface = "bold"
  ) +
  geom_segment(
    data = sig_period,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_segment(
    data = sig_period,
    aes(x = xstart, xend = xstart, y = y_low, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_segment(
    data = sig_period,
    aes(x = xend, xend = xend, y = y_low, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_text(
    data = sig_period,
    aes(x = x, y = y_text, label = stars),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = legend_period,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  expand_limits(y = max(sig_period$y_text, na.rm = TRUE) * 1.25) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Wild type and double mutant prevalence\nbetween 2021-2023 and 2024-2025",
    x = NULL,
    y = "Prevalence (%)"
  ) +
  theme(
    axis.title = element_text(size = 12, colour = "black", face = "bold"),
    axis.text = element_text(size = 11, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    plot.title = ggtext::element_textbox_simple(
      size = 14,
      face = "bold",
      halign = 0.5,
      fill = "grey75",
      box.color = "black",
      padding = margin(6, 6, 6, 6),
      margin = margin(0, 0, 12, 0)
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )

print(K13_period_comparison_KFH)


# 6. Kigali and hotspot analysis for 2024-2025


df_2425_kigali <- df %>%
  filter(YEAR %in% c(2024, 2025), District %in% c("Gasabo", "Kicukiro", "Nyarugenge"))

df_2425_hotspot <- df %>%
  filter(
    YEAR %in% c(2024, 2025),
    District %in% c("Bugesera", "Gisagara", "Kirehe", "Karongi", "Nyagatare", "Musanze", "Rutsiro")
  )

k13_cat_2425_kigali <- make_k13_categories(df_2425_kigali)
k13_cat_2425_hotspot <- make_k13_categories(df_2425_hotspot)

total_kigali <- df_2425_kigali %>% distinct(SAMPLE) %>% nrow()
total_hotspot <- df_2425_hotspot %>% distinct(SAMPLE) %>% nrow()

counts_2425_kigali <- add_wild_type(k13_cat_2425_kigali, total_kigali)
counts_2425_hotspot <- add_wild_type(k13_cat_2425_hotspot, total_hotspot)

kigali_results <- run_wt_double_binom(counts_2425_kigali)
hotspot_results <- run_wt_double_binom(counts_2425_hotspot)

wild_kigali <- kigali_results$wild_n
mut_kigali <- kigali_results$double_n
kigali_test <- kigali_results$test

wild_hotspot <- hotspot_results$wild_n
mut_hotspot <- hotspot_results$double_n
hotspot_test <- hotspot_results$test

cat(
  "Kigali 2024-2025:",
  "Wild type =", kigali_results$wild_pct, "% vs mutant =", kigali_results$double_pct,
  "%, p =", kigali_results$p_value, "\n"
)

cat(
  "Hotspots 2024-2025:",
  "Wild type =", hotspot_results$wild_pct, "% vs mutant =", hotspot_results$double_pct,
  "%, p =", hotspot_results$p_value, "\n"
)

plot_df_all <- bind_rows(
  counts_2123 %>% mutate(period = "KFH 2021-2023"),
  counts_2425_kigali %>% mutate(period = "Kigali 2024-2025"),
  counts_2425_hotspot %>% mutate(period = "Hotspots 2024-2025")
) %>%
  group_by(period) %>%
  mutate(prevalence = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = category_order)) %>%
  left_join(
    tibble(
      period = c("KFH 2021-2023", "Kigali 2024-2025", "Hotspots 2024-2025"),
      total_n = c(total_samples_2123, total_kigali, total_hotspot)
    ),
    by = "period"
  ) %>%
  mutate(
    period_label = paste0(period, "(N = ", total_n, ")"),
    period_label = factor(
      period_label,
      levels = c(
        paste0("KFH 2021-2023", "(N = ", total_samples_2123, ")"),
        paste0("Kigali 2024-2025", "(N = ", total_kigali, ")"),
        paste0("Hotspots 2024-2025", "(N = ", total_hotspot, ")")
      )
    )
  )

write.table(
  plot_df_all,
  file = file.path(output_dir, "figure5.tsv"),
  row.names = FALSE,
  sep = "\t"
)

sig_all <- plot_df_all %>%
  group_by(period, period_label) %>%
  summarise(y_max = max(prevalence, na.rm = TRUE), .groups = "drop") %>%
  arrange(period_label) %>%
  mutate(
    y = y_max * 1.12,
    y_low = y_max * 1.05,
    y_text = y_max * 1.16,
    xstart = 1,
    xend = case_when(
      grepl("Hotspots", period_label) ~ 5,
      TRUE ~ 6
    ),
    p_value = case_when(
      grepl("^KFH 2021", period_label) ~ k13_prop_test$p.value,
      grepl("^Kigali 2024", period_label) ~ kigali_test$p.value,
      grepl("^Hotspots 2024", period_label) ~ hotspot_test$p.value
    ),
    stars = p_to_stars(p_value),
    label = stars,
    legend_line = paste0(stars, " : p = ", signif(p_value, 3))
  )

K13_evolution_all <- ggplot(plot_df_all, aes(x = category, y = prevalence, fill = category)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_segment(
    data = sig_all,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_segment(
    data = sig_all,
    aes(x = xstart, xend = xstart, y = y_low, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_segment(
    data = sig_all,
    aes(x = xend, xend = xend, y = y_low, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_text(
    data = sig_all,
    aes(x = (xstart + xend) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  geom_text(
    data = sig_all,
    aes(x = xend + 0.25, y = y_text + 4, label = legend_line),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(~ period_label, scales = "free", nrow = 1) +
  scale_fill_manual(values = k13_colors) +
  expand_limits(y = max(sig_all$y_text, na.rm = TRUE) * 1.15) +
  coord_cartesian(clip = "off") +
  labs(
    x = "K13 mutation category",
    y = "Prevalence (%)"
  ) +
  theme(
    axis.title = element_text(size = 10, colour = "black", face = "bold"),
    axis.text = element_text(size = 7, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

print(K13_evolution_all)


# 7. Kigali vs hotspot comparison in 2024-2025


location_comparison_df <- bind_rows(
  counts_2425_kigali %>%
    mutate(location_group = "Kigali 2024-2025", total_n = total_kigali),
  counts_2425_hotspot %>%
    mutate(location_group = "Hotspots 2024-2025", total_n = total_hotspot)
) %>%
  filter(category %in% c("Wild type", "R561H + K189T")) %>%
  mutate(
    prevalence = n / total_n * 100,
    category = factor(category, levels = c("Wild type", "R561H + K189T")),
    location_group = factor(
      location_group,
      levels = c("Kigali 2024-2025", "Hotspots 2024-2025")
    )
  )

wt_kigali_n <- location_comparison_df %>%
  filter(category == "Wild type", location_group == "Kigali 2024-2025") %>%
  pull(n)

wt_hotspots_n <- location_comparison_df %>%
  filter(category == "Wild type", location_group == "Hotspots 2024-2025") %>%
  pull(n)

double_kigali_n <- location_comparison_df %>%
  filter(category == "R561H + K189T", location_group == "Kigali 2024-2025") %>%
  pull(n)

double_hotspots_n <- location_comparison_df %>%
  filter(category == "R561H + K189T", location_group == "Hotspots 2024-2025") %>%
  pull(n)

total_kigali_2425 <- location_comparison_df %>%
  filter(location_group == "Kigali 2024-2025") %>%
  distinct(total_n) %>%
  pull(total_n)

total_hotspots_2425 <- location_comparison_df %>%
  filter(location_group == "Hotspots 2024-2025") %>%
  distinct(total_n) %>%
  pull(total_n)

test_wt_location <- fisher.test(matrix(
  c(
    wt_kigali_n, total_kigali_2425 - wt_kigali_n,
    wt_hotspots_n, total_hotspots_2425 - wt_hotspots_n
  ),
  nrow = 2,
  byrow = TRUE
))

test_double_location <- fisher.test(matrix(
  c(
    double_kigali_n, total_kigali_2425 - double_kigali_n,
    double_hotspots_n, total_hotspots_2425 - double_hotspots_n
  ),
  nrow = 2,
  byrow = TRUE
))

sig_location <- location_comparison_df %>%
  group_by(category) %>%
  summarise(y_max = max(prevalence, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    x = as.numeric(category),
    xstart = x - 0.22,
    xend = x + 0.22,
    y = y_max * 1.15,
    y_low = y_max * 1.07,
    y_text = y_max * 1.22,
    p_value = case_when(
      category == "Wild type" ~ test_wt_location$p.value,
      category == "R561H + K189T" ~ test_double_location$p.value
    ),
    stars = p_to_stars(p_value),
    legend_line = paste0(category, " : ", stars, " p = ", signif(p_value, 3))
  )

legend_location <- paste(sig_location$legend_line, collapse = "\n")

K13_location_comparison_2425 <- ggplot(
  location_comparison_df,
  aes(x = category, y = prevalence, fill = category, group = location_group)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "black") +
  scale_fill_manual(values = c("Wild type" = "#2ecc71", "R561H + K189T" = "#ff7f0e")) +
  geom_text(
    aes(label = location_group),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 3,
    fontface = "bold"
  ) +
  geom_segment(
    data = sig_location,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_segment(
    data = sig_location,
    aes(x = xstart, xend = xstart, y = y_low, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_segment(
    data = sig_location,
    aes(x = xend, xend = xend, y = y_low, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_text(
    data = sig_location,
    aes(x = x, y = y_text, label = stars),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = legend_location,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  expand_limits(y = max(sig_location$y_text, na.rm = TRUE) * 1.25) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Wild type and double mutant prevalence<br>between Kigali and Hotspots in 2024-2025",
    x = NULL,
    y = "Prevalence (%)"
  ) +
  theme(
    axis.title = element_text(size = 12, colour = "black", face = "bold"),
    axis.text = element_text(size = 11, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = ggtext::element_textbox_simple(
      size = 14,
      face = "bold",
      halign = 0.5,
      fill = "grey75",
      box.color = "black",
      padding = margin(6, 6, 6, 6),
      margin = margin(0, 0, 12, 0)
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )

print(K13_location_comparison_2425)


# 8. Export figures


ggsave(
  "Evolution on k13 mutations in Rwanda_all.svg",
  plot = K13_evolution_all,
  path = output_dir,
  width = 20,
  height = 7,
  units = "in",
  dpi = 600
)

ggsave(
  "Evolution on k13 mutations in kigali.svg",
  plot = K13_evolution_KFH,
  path = output_dir,
  width = 15,
  height = 7,
  units = "in",
  dpi = 600
)

ggsave(
  "Comparison_wildtype_double_mutant_periods.svg",
  plot = K13_period_comparison_KFH,
  path = output_dir,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)

ggsave(
  "Comparison_wildtype_double_mutant_Kigali_vs_Hotspots_2024_2025.svg",
  plot = K13_location_comparison_2425,
  path = output_dir,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)


# 9. Confidence intervals for the binomial tests

# The binomial test gives the 95% CI for the wild-type proportion.
# For the double mutant, the CI is obtained as the complementary interval.

print_binom_ci(
  test_obj = k13_prop_test,
  wild_count = wild_n,
  mutant_count = r561h_k189t_n,
  label = "KFH 2021-2023"
)

print_binom_ci(
  test_obj = k13_prop_test_2425,
  wild_count = wild_n_2425,
  mutant_count = r561h_k189t_n_2425,
  label = "KFH 2024-2025"
)

print_binom_ci(
  test_obj = kigali_test,
  wild_count = wild_kigali,
  mutant_count = mut_kigali,
  label = "Kigali 2024-2025"
)

print_binom_ci(
  test_obj = hotspot_test,
  wild_count = wild_hotspot,
  mutant_count = mut_hotspot,
  label = "Hotspots 2024-2025"
)


# 10. Logistic regression: Wild type vs R561H + K189T

# Here the model is used as another way to compare the double-mutant
# proportion with the wild-type proportion in each group.

logistic_input <- tibble(
  comparison = c(
    "KFH 2021-2023",
    "KFH 2024-2025",
    "Kigali 2024-2025",
    "Hotspots 2024-2025"
  ),
  wild_type = c(
    wild_n,
    wild_n_2425,
    wild_kigali,
    wild_hotspot
  ),
  double_mutant = c(
    r561h_k189t_n,
    r561h_k189t_n_2425,
    mut_kigali,
    mut_hotspot
  )
)

run_logistic_regression <- function(wild_type, double_mutant) {
  model_data <- tibble(
    double_mutant = double_mutant,
    wild_type = wild_type
  )
  
  model <- glm(
    cbind(double_mutant, wild_type) ~ 1,
    data = model_data,
    family = binomial(link = "logit")
  )
  
  model_summary <- summary(model)$coefficients
  
  log_odds <- model_summary[1, "Estimate"]
  se <- model_summary[1, "Std. Error"]
  z_value <- model_summary[1, "z value"]
  p_value <- model_summary[1, "Pr(>|z|)"]
  
  odds_ratio <- exp(log_odds)
  ci_low <- exp(log_odds - 1.96 * se)
  ci_high <- exp(log_odds + 1.96 * se)
  
  tibble(
    log_odds = log_odds,
    odds_ratio = odds_ratio,
    ci_low = ci_low,
    ci_high = ci_high,
    z_value = z_value,
    p_value = p_value
  )
}

logistic_results <- logistic_input %>%
  rowwise() %>%
  mutate(
    total_compared = wild_type + double_mutant,
    wild_type_percent = 100 * wild_type / total_compared,
    double_mutant_percent = 100 * double_mutant / total_compared,
    model_result = list(run_logistic_regression(wild_type, double_mutant))
  ) %>%
  unnest(model_result) %>%
  ungroup() %>%
  mutate(
    odds_ratio = round(odds_ratio, 3),
    ci_low = round(ci_low, 3),
    ci_high = round(ci_high, 3),
    wild_type_percent = round(wild_type_percent, 2),
    double_mutant_percent = round(double_mutant_percent, 2),
    p_value = signif(p_value, 3),
    interpretation = case_when(
      odds_ratio > 1 & p_value < 0.05 ~ "R561H + K189T higher than wild type",
      odds_ratio < 1 & p_value < 0.05 ~ "Wild type higher than R561H + K189T",
      TRUE ~ "No significant difference"
    )
  ) %>%
  select(
    comparison,
    wild_type,
    double_mutant,
    total_compared,
    wild_type_percent,
    double_mutant_percent,
    odds_ratio,
    ci_low,
    ci_high,
    p_value,
    interpretation
  )

print(logistic_results)

print_logistic_pvalue <- function(label, wild_type, double_mutant) {
  model <- glm(
    cbind(double_mutant, wild_type) ~ 1,
    family = binomial(link = "logit")
  )
  
  p_value <- summary(model)$coefficients[1, "Pr(>|z|)"]
  
  cat("\n")
  cat(label, "\n")
  cat("p-value =", format(p_value, scientific = TRUE, digits = 16), "\n")
}

cat("\n============================================\n")
cat("Logistic regression: WT vs R561H + K189T\n")
cat("============================================\n")

print_logistic_pvalue(
  label = "KFH 2021-2023",
  wild_type = wild_n,
  double_mutant = r561h_k189t_n
)

print_logistic_pvalue(
  label = "KFH 2024-2025",
  wild_type = wild_n_2425,
  double_mutant = r561h_k189t_n_2425
)

print_logistic_pvalue(
  label = "Kigali 2024-2025",
  wild_type = wild_kigali,
  double_mutant = mut_kigali
)

print_logistic_pvalue(
  label = "Hotspots 2024-2025",
  wild_type = wild_hotspot,
  double_mutant = mut_hotspot
)

write.table(
  logistic_results,
  file = file.path(output_dir, "logistic_regression_WT_vs_R561H_K189T.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
