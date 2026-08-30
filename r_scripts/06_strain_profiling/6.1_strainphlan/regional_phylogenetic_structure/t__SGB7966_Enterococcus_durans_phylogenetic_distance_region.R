# StrainPhlAn t__SGB7966 Enterococcus durans Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/regional_phylogenetic_structure
# 17/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkarian' = '#4FAFA8', 'Karachay-Cherkess' = '#E8403D')


# Load the data -----------------------------------------------------------

edurans_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7966_Enterococcus_durans_phylogenetic_distance_new.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
edurans_mat_df <- edurans_data |>
  bind_rows(
    edurans_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
edurans_all_ids <- sort(unique(c(edurans_mat_df$strain_id_1, edurans_mat_df$strain_id_2)))

edurans_diag_df <- tibble(
  strain_id_1 = edurans_all_ids,
  strain_id_2 = edurans_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
edurans_mat_df <- bind_rows(edurans_diag_df, edurans_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)


# II. Within-Region vs Between-Region Structure ---------------------------

# Long-format distances for REGION comparisons (unique pairs only; remove self)
edurans_dist_long <- edurans_mat_df |>
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
edurans_dist_long$comparison_type <- factor(
  edurans_dist_long$comparison_type,
  levels = c('Within KB', 'Within KC', 'Between Regions')
)

# Summary table (prints in console)
edurans_region_summary <- edurans_dist_long |>
  group_by(comparison_type) |>
  summarise(
    n_comparisons = n(),
    mean_distance = mean(phylogenetic_distance),
    sd_distance   = sd(phylogenetic_distance),
    .groups = 'drop'
  )

print(edurans_region_summary)

# Optional: save the summary table
write_tsv(edurans_region_summary, 'data/results/06_strain_profiling/edurans_region_structure_summary.tsv')

# Overall test (prints in console)
edurans_region_kw <- kruskal.test(
  phylogenetic_distance ~ comparison_type,
  data = edurans_dist_long
)

print(edurans_region_kw)

# Preparation for saving
edurans_region_kw_df <- tibble(
  test = 'Kruskal-Wallis',
  chi_squared = edurans_region_kw$statistic,
  df = edurans_region_kw$parameter,
  p_value = edurans_region_kw$p.value
)

write_tsv(edurans_region_kw_df, 'data/results/06_strain_profiling/edurans_region_kw_test.tsv')

# Format Kruskal–Wallis p-value
edurans_kw_p <- signif(edurans_region_kw$p.value, 3)

# Create plotmath label (italic p)
edurans_kw_label <- paste0(
  "'Kruskal-Wallis: ' ~ italic(p) == ",
  edurans_kw_p
)


# Boxplot of Distances by Comparison Class --------------------------------

edurans_region_plot <- ggplot(
  edurans_dist_long,
  aes(x = comparison_type, y = phylogenetic_distance, fill = comparison_type)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # Add Kruskal–Wallis annotation
  annotate(
    'text',
    x = 2,
    y = max(edurans_dist_long$phylogenetic_distance) * 1.05,
    label = edurans_kw_label,
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
    # title = 'Regional Phylogenetic Structure',
    title = '(A)',
    subtitle = expression(bolditalic('Enterococcus durans')),
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
edurans_region_plot

# Save the plot
ggsave(
  edurans_region_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_region/edurans_region_structure_boxplot.tiff',
  height = 15,
  width = 20,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)