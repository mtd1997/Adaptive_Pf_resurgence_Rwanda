
# Comparison of survival (%) at 700 nM DHA between KF strains
# The script follows the same calculation strategy as the previous plot:
# first the fluorescence values are corrected with the blank,
# then they are normalized using the No drug condition as the 100% reference.
#
# The main goal is to compare the survival of the sensitive strains
# (KFHS49 and KFHS50) with the other KF strains after treatment with 700 nM DHA.



# Short description of the analysis
 
# This script analyses fluorescence data obtained for different KF strains.
# The fluorescence signal is first corrected by subtracting the mean blank value.
# Then, survival is expressed as a percentage of the No drug condition, which is
# used as the untreated control.
#
# The final objective is to visualize and compare the survival percentages
# between: 1) No drug, 2) 700 nM DHA for the sensitive strains KFHS49/KFHS50,
# and 3) 700 nM DHA for the other strains considered as partial resistant strains.
# The script also performs t-tests and adds the corresponding significance
# stars directly on the figure.

# List of packages required for reading Excel files, manipulating data,
# performing reshaping steps, and creating the final plot.
packages <- c("readxl", "dplyr", "ggplot2", "purrr", "stringr", "tidyr")

# This loop loads each package. If a package is missing, it is installed first.
# This makes the script easier to run on another computer.
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Folder containing the Excel files exported from the fluorescence experiment.
# This path must be changed if the files are stored in another folder.
input_dir <- "/Users/diallo/Desktop/fluodata"

# All Excel files in the folder are collected automatically.
# full.names = TRUE keeps the complete file path, which is useful when reading the files.
excel_files <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)

# The script stops here if no Excel file is found.
# This avoids running the full analysis on an empty input folder.
if (length(excel_files) == 0) {
  stop("No .xlsx file was found in the fluodata folder.")
}


# Function used to read and clean one strain file


read_one_file <- function(file) {
  
  # The Excel sheet is read without column names because the useful table
  # starts lower in the file and does not contain standard headers.
  raw <- readxl::read_excel(file, sheet = "Results", col_names = FALSE)
  
  # The strain name is taken from the file name, without the .xlsx extension.
  strain_name <- tools::file_path_sans_ext(basename(file))
  
  # Only the rows containing the values of interest are selected.
  # Generic column names are added to make the next steps easier.
  dat <- raw[11:17, ]
  colnames(dat) <- paste0("V", 1:ncol(dat))
  
  # The useful columns are renamed with clear names.
  # The fluorescence values are converted to numeric format to allow calculations.
  df_wide <- dat %>%
    dplyr::transmute(
      replicate = as.character(V5),
      sample = as.character(V6),
      d1 = as.numeric(V7),
      d2 = as.numeric(V8),
      d3 = as.numeric(V9),
      d4 = as.numeric(V10),
      d5 = as.numeric(V11),
      d6 = as.numeric(V12),
      d7 = as.numeric(V13),
      d8 = as.numeric(V14),
      d9 = as.numeric(V15),
      blank = as.numeric(V16),
      control_value = as.numeric(V17),
      control_label = as.character(V18)
    )
  
  # I use the mean blank value from the extract wells as the background signal.
  df_extracts <- df_wide %>%
    dplyr::filter(sample %in% c("Infusion A", "Lyoph.A Afra"))
  
  global_blank_mean <- mean(df_extracts$blank, na.rm = TRUE)
  
  # The corrected No drug condition is used as the reference value, corresponding to 100% survival.
  mean_no_drug <- df_wide %>%
    dplyr::filter(!is.na(control_label), control_label == "No drug") %>%
    dplyr::mutate(control_net = control_value - global_blank_mean) %>%
    dplyr::summarise(mean_no_drug = mean(control_net, na.rm = TRUE)) %>%
    dplyr::pull(mean_no_drug)
  
  # Only the two conditions needed for this comparison are kept: No drug and 700 nM DHA.
  df_conditions <- df_wide %>%
    dplyr::filter(control_label %in% c("No drug", "700 nM DHA")) %>%
    dplyr::transmute(
      file_name = basename(file),
      strain = strain_name,
      replicate = replicate,
      condition = control_label,
      fluorescence = control_value,
      global_blank_mean = global_blank_mean,
      mean_no_drug = mean_no_drug,
      fluorescence_net = fluorescence - global_blank_mean,
      survival_pct = ((fluorescence - global_blank_mean) / mean_no_drug) * 100
    )
  
  return(df_conditions)
}


# Read all strain files and prepare the final dataset


# The same cleaning function is applied to all Excel files.
# map_df() combines the results into one single table.
df_survival_all <- purrr::map_df(excel_files, read_one_file) %>%
  dplyr::filter(
    !is.na(strain),
    strain != "",
    stringr::str_detect(strain, "^KF"),
    !is.na(survival_pct),
    condition %in% c("No drug", "700 nM DHA")
  ) %>%
  dplyr::mutate(
    # The samples are separated into three groups for the plot:
    # No drug, sensitive strains under DHA, and partial resistant strains under DHA.
    dha_group = dplyr::case_when(
      condition == "No drug" ~ "No drug",
      condition == "700 nM DHA" & strain %in% c("KFHS49", "KFHS50") ~ "700 nM DHA\nSensitives",
      condition == "700 nM DHA" & !strain %in% c("KFHS49", "KFHS50") ~ "700 nM DHA\nPartial resistants"
    ),
    dha_group = factor(
      dha_group,
      levels = c(
        "No drug",
        "700 nM DHA\nSensitives",
        "700 nM DHA\nPartial resistants"
      )
    ),
    condition = factor(condition, levels = c("No drug", "700 nM DHA"))
  )

print(df_survival_all)

df_survival_all <- df_survival_all %>%
  dplyr::mutate(
    x_pos = dplyr::case_when(
      dha_group == "No drug" ~ 1,
      dha_group == "700 nM DHA\nSensitives" ~ 1.65,
      dha_group == "700 nM DHA\nPartial resistants" ~ 2.30
    )
  )



# Define one color for each strain

strain_levels <- sort(unique(as.character(df_survival_all$strain)))

vivid_palette <- c(
  "#E41A1C",  # red
  "#377EB8",  # blue
  "#4DAF4A",  # green
  "#984EA3",  # purple
  "#FF7F00",  # orange
  "#FFFF33",  # yellow
  "#A65628",  # brown
  "#F781BF",  # pink
  "#00BFC4",  # cyan
  "#7CAE00",  # bright green
  "#C77CFF",  # violet
  "#F8766D"   # coral
)

# If there are more strains than colors in the base palette,
# colorRampPalette() creates enough intermediate colors.
custom_strain_colors <- setNames(
  colorRampPalette(vivid_palette)(length(strain_levels)),
  strain_levels
)


# Prepare strain-level means: one point corresponds to the mean of the triplicates for one strain


# For the statistical analysis, the triplicates are summarized per strain.
# This avoids treating technical replicates as fully independent biological observations.
df_survival_strain_mean <- df_survival_all %>%
  dplyr::group_by(strain, condition, dha_group) %>%
  dplyr::summarise(
    mean_survival = mean(survival_pct, na.rm = TRUE),
    sd_survival = sd(survival_pct, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

print(df_survival_strain_mean)

df_survival_strain_mean <- df_survival_strain_mean %>%
  dplyr::mutate(
    x_pos = dplyr::case_when(
      dha_group == "No drug" ~ 1,
      dha_group == "700 nM DHA\nSensitives" ~ 1.65,
      dha_group == "700 nM DHA\nPartial resistants" ~ 2.30
    )
  )


# Helper function used to convert p-values into significance stars


# This small function gives the usual star notation used on the figure.
p_to_star <- function(p) {
  if (is.na(p)) {
    "ns"
  } else if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else {
    "ns"
  }
}


# Statistical tests used for the three main comparisons


df_test_split <- df_survival_strain_mean %>%
  dplyr::mutate(
    # The samples are separated into three groups for the plot:
    # No drug, sensitive strains under DHA, and partial resistant strains under DHA.
    dha_group = dplyr::case_when(
      condition == "No drug" ~ "No drug",
      condition == "700 nM DHA" & strain %in% c("KFHS49", "KFHS50") ~ "700 nM DHA\nKFHS49 + KFHS50",
      condition == "700 nM DHA" & !strain %in% c("KFHS49", "KFHS50") ~ "700 nM DHA\nOther strains"
    )
  )

# The table is converted to wide format so that each strain has
# one column for No drug and one column for 700 nM DHA.
df_test_wide <- df_test_split %>%
  dplyr::select(strain, condition, mean_survival) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = mean_survival
  ) %>%
  dplyr::mutate(
    dha_category = dplyr::if_else(
      strain %in% c("KFHS49", "KFHS50"),
      "700 nM DHA\nKFHS49 + KFHS50",
      "700 nM DHA\nOther strains"
    )
  )

# Test 1: paired comparison between No drug and 700 nM DHA for KFHS49 and KFHS50
test_no_drug_vs_700_ref <- t.test(
  df_test_wide$`No drug`[df_test_wide$dha_category == "700 nM DHA\nKFHS49 + KFHS50"],
  df_test_wide$`700 nM DHA`[df_test_wide$dha_category == "700 nM DHA\nKFHS49 + KFHS50"],
  paired = TRUE
)

# Test 2: paired comparison between No drug and 700 nM DHA for the other strains
test_no_drug_vs_700_other <- t.test(
  df_test_wide$`No drug`[df_test_wide$dha_category == "700 nM DHA\nOther strains"],
  df_test_wide$`700 nM DHA`[df_test_wide$dha_category == "700 nM DHA\nOther strains"],
  paired = TRUE
)

# Test 3: unpaired comparison between sensitive strains and other strains under 700 nM DHA
test_700_ref_vs_700_other <- t.test(
  df_test_wide$`700 nM DHA`[df_test_wide$dha_category == "700 nM DHA\nKFHS49 + KFHS50"],
  df_test_wide$`700 nM DHA`[df_test_wide$dha_category == "700 nM DHA\nOther strains"],
  paired = FALSE
)

# The p-values are extracted from the test objects and then converted into stars.
p1 <- test_no_drug_vs_700_ref$p.value
p2 <- test_no_drug_vs_700_other$p.value
p3 <- test_700_ref_vs_700_other$p.value

lab1 <- p_to_star(p1)
lab2 <- p_to_star(p2)
lab3 <- p_to_star(p3)

print(test_no_drug_vs_700_ref)
print(test_no_drug_vs_700_other)
print(test_700_ref_vs_700_other)

cat("No drug vs 700 nM DHA KFHS49+KFHS50 p =", p1, "\n")
cat("No drug vs 700 nM DHA Other strains p =", p2, "\n")
cat("700 nM DHA KFHS49+KFHS50 vs Other strains p =", p3, "\n")

# Automatic y-axis limits and positions for the statistical annotations


# The y-axis range is calculated from the data instead of being fixed manually.
# This makes the plot more robust if new values are added later.
y_min_plot_A <- min(df_survival_all$survival_pct, na.rm = TRUE)
y_max_plot_A <- max(df_survival_all$survival_pct, na.rm = TRUE)
y_range_plot_A <- y_max_plot_A - y_min_plot_A

y_upper_limit_A <- y_max_plot_A + 0.90 * y_range_plot_A
y_lower_limit_A <- y_min_plot_A - 0.08 * y_range_plot_A

p_bar_y_A  <- y_upper_limit_A - 0.22 * y_range_plot_A
p_tick_y_A <- y_upper_limit_A - 0.26 * y_range_plot_A
p_text_y_A <- y_upper_limit_A - 0.18 * y_range_plot_A

# Each point represents one individual value from the Excel files

# Calculate the mean and standard deviation shown above each group
labels_replicates <- df_survival_all %>%
  dplyr::group_by(dha_group, x_pos) %>%
  dplyr::summarise(
    mean_survival = mean(survival_pct, na.rm = TRUE),
    sd_survival = sd(survival_pct, na.rm = TRUE),
    y_label = max(survival_pct, na.rm = TRUE) + 8,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = paste0(
      "Mean = ", round(mean_survival, 1),
      "%\nSD = ", round(sd_survival, 1)
    )
  )

# Final plot showing the individual replicate values, the boxplots,
# the mean/SD labels, and the statistical comparisons.
plot_points_replicates <- ggplot(
  df_survival_all,
  aes(x = x_pos, y = survival_pct, color = strain))+
  geom_boxplot(
    aes(group = x_pos),
    width = 0.28,
    outlier.shape = NA,
    fill = "grey90",
    color = "black",
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.2,
    alpha = 0.85,
    shape = 16
  ) +
  # Comparison 1: No drug vs 700 nM DHA KFHS49 + KFHS50
  # Comparison 1: No drug vs 700 nM DHA Sensitives
  annotate(
    "segment",
    x = 1,
    xend = 1.65,
    y = p_bar_y_A,
    yend = p_bar_y_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = p_tick_y_A,
    yend = p_bar_y_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "segment",
    x = 1.65,
    xend = 1.65,
    y = p_tick_y_A,
    yend = p_bar_y_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "text",
    x = 1.325,
    y = p_text_y_A,
    label = lab1,
    size = 8.5,
    fontface = "bold",
    color = "black"
  ) +
  
  # Comparison 2: No drug vs 700 nM DHA Partial resistants
  annotate(
    "segment",
    x = 1,
    xend = 2.30,
    y = p_bar_y_A + 0.13 * y_range_plot_A,
    yend = p_bar_y_A + 0.13 * y_range_plot_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = p_tick_y_A + 0.13 * y_range_plot_A,
    yend = p_bar_y_A + 0.13 * y_range_plot_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "segment",
    x = 2.30,
    xend = 2.30,
    y = p_tick_y_A + 0.13 * y_range_plot_A,
    yend = p_bar_y_A + 0.13 * y_range_plot_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "text",
    x = 1.65,
    y = p_text_y_A + 0.13 * y_range_plot_A,
    label = lab2,
    size = 8.5,
    fontface = "bold",
    color = "black"
  ) +
  
  # Comparison 3: 700 nM DHA Sensitives vs 700 nM DHA Partial resistants
  annotate(
    "segment",
    x = 1.65,
    xend = 2.30,
    y = p_bar_y_A - 0.13 * y_range_plot_A,
    yend = p_bar_y_A - 0.13 * y_range_plot_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "segment",
    x = 1.65,
    xend = 1.65,
    y = p_tick_y_A - 0.13 * y_range_plot_A,
    yend = p_bar_y_A - 0.13 * y_range_plot_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "segment",
    x = 2.30,
    xend = 2.30,
    y = p_tick_y_A - 0.13 * y_range_plot_A,
    yend = p_bar_y_A - 0.13 * y_range_plot_A,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "text",
    x = 1.975,
    y = p_text_y_A - 0.13 * y_range_plot_A,
    label = lab3,
    size = 6.5,
    fontface = "bold",
    color = "black"
  ) +
  geom_text(
    data = labels_replicates,
    aes(
      x = x_pos,
      y = y_label,
      label = label
    ),
    inherit.aes = FALSE,
    size = 6.5,
    fontface = "bold",
    color = "black",
    vjust = 0,
    lineheight = 1.05
  ) +
  scale_color_manual(
    values = custom_strain_colors,
    name = "Strain"
  ) +
  scale_x_continuous(
    breaks = c(1, 1.65, 2.30),
    labels = c(
      "No drug",
      "700 nM DHA\nSensitives",
      "700 nM DHA\nPartial resistants"
    ),
    limits = c(0.7, 2.6)
  ) +
  scale_y_continuous(
    limits = c(y_lower_limit_A, y_upper_limit_A),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Distribution of survival (%) at No drug and 700 nM DHA",
    x = "",
    y = "Survival (%)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 20),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(size = 14)
  )

print(plot_points_replicates)

# The final figure is saved in the same folder as the input files.
ggsave(
  filename = "plot_survival_points_replicates_NoDrug_vs_700nMDHA.png",
  plot = plot_points_replicates,
  path = input_dir,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)