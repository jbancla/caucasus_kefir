# StrainPhlAn t__SGB7985 Lactococcus lactis Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/regional_phylogenetic_structure
# 17/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkarian' = '#4FAFA8', 'Karachay-Cherkess' = '#E8403D')


# Load the data -----------------------------------------------------------

llactis_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7985_Lactococcus_lactis_phylogenetic_distance.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
llactis_mat_df <- llactis_data |>
  bind_rows(
    llactis_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
llactis_all_ids <- sort(unique(c(llactis_mat_df$strain_id_1, llactis_mat_df$strain_id_2)))

llactis_diag_df <- tibble(
  strain_id_1 = llactis_all_ids,
  strain_id_2 = llactis_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
llactis_mat_df <- bind_rows(llactis_diag_df, llactis_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)


# II. Within-Region vs Between-Region Structure ---------------------------

# Long-format distances for REGION comparisons (unique pairs only; remove self)
llactis_dist_long <- llactis_mat_df |>
  filter(strain_id_1 < strain_id_2) |>
  mutate(
    region_1 = if_else(
      strain_id_1 %in% kab_samples,
      'Kabardino-Balkarian',
      'Karachay-Cherkess'
    ),
    region_2 = if_else(
      strain_id_2 %in% kab_samples,
      'Kabardino-Balkarian',
      'Karachay-Cherkess'
    ),
    comparison_type = case_when(
      region_1 == 'Kabardino-Balkarian' & region_2 == 'Kabardino-Balkarian' ~ 'Within KB',
      region_1 == 'Karachay-Cherkess'   & region_2 == 'Karachay-Cherkess'   ~ 'Within KC',
      TRUE ~ 'Between Regions'
    )
  ) |>
  filter(phylogenetic_distance > 0)

# Order comparison classes
llactis_dist_long$comparison_type <- factor(
  llactis_dist_long$comparison_type,
  levels = c('Within KB', 'Within KC', 'Between Regions')
)

# Summary table (prints in console)
llactis_region_summary <- llactis_dist_long |>
  group_by(comparison_type) |>
  summarise(
    n_comparisons = n(),
    mean_distance = mean(phylogenetic_distance),
    sd_distance   = sd(phylogenetic_distance),
    .groups = 'drop'
  )

print(llactis_region_summary)

# Optional: save the summary table
write_tsv(llactis_region_summary, 'data/results/06_strain_profiling/llactis_region_structure_summary.tsv')

# Overall test (prints in console)
llactis_region_kw <- kruskal.test(
  phylogenetic_distance ~ comparison_type,
  data = llactis_dist_long
)

print(llactis_region_kw)

# Preparation for saving
llactis_region_kw_df <- tibble(
  test = 'Kruskal-Wallis',
  chi_squared = llactis_region_kw$statistic,
  df = llactis_region_kw$parameter,
  p_value = llactis_region_kw$p.value
)

write_tsv(llactis_region_kw_df, 'data/results/06_strain_profiling/llactis_region_kw_test.tsv')

# Format Kruskal–Wallis p-value
llactis_kw_p <- signif(llactis_region_kw$p.value, 3)

# Create plotmath label (italic p)
llactis_kw_label <- paste0(
  "'Kruskal-Wallis: ' ~ italic(p) == ",
  llactis_kw_p
)


# Boxplot of Distances by Comparison Class --------------------------------

llactis_region_plot <- ggplot(
  llactis_dist_long,
  aes(x = comparison_type, y = phylogenetic_distance, fill = comparison_type)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # Add Kruskal–Wallis annotation
  annotate(
    'text',
    x = 2,
    y = max(llactis_dist_long$phylogenetic_distance) * 1.05,
    label = llactis_kw_label,
    parse = TRUE,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      'Within KB' = '#4FAFA8',
      'Within KC' = '#E8403D',
      'Between Regions' = 'grey60'
    )
  ) +
  
  labs(
    title = '(B)',
    subtitle = expression(bolditalic('Lactococcus lactis')),
    x = NULL,
    y = 'Phylogenetic Distance'
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(face = 'bold', size = 12, margin = margin(r = 10)),
    legend.position = 'none'
  )

# Check the plot
llactis_region_plot

# Save the plot
ggsave(
  llactis_region_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_region/llactis_region_structure_boxplot.tiff',
  height = 15,
  width = 20,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)