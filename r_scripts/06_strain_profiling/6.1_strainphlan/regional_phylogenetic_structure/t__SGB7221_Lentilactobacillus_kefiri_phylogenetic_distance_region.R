# StrainPhlAn t__SGB7221 Lentilactobacillus kefiri Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/regional_phylogenetic_structure
# 17/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c(
  'k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
  'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g'
)

site_cols <- c('Kabardino-Balkarian' = '#4FAFA8', 'Karachay-Cherkess' = '#E8403D')


# Load the data -----------------------------------------------------------

lkefiri_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7221_Lentilactobacillus_kefiri_phylogenetic_distance.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
lkefiri_mat_df <- lkefiri_data |>
  bind_rows(
    lkefiri_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
lkefiri_all_ids <- sort(unique(c(lkefiri_mat_df$strain_id_1, lkefiri_mat_df$strain_id_2)))

lkefiri_diag_df <- tibble(
  strain_id_1 = lkefiri_all_ids,
  strain_id_2 = lkefiri_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
lkefiri_mat_df <- bind_rows(lkefiri_diag_df, lkefiri_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)


# II. Within-Region vs Between-Region Structure ---------------------------

# Long-format distances for REGION comparisons (unique pairs only; remove self)
lkefiri_dist_long <- lkefiri_mat_df |>
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
lkefiri_dist_long$comparison_type <- factor(
  lkefiri_dist_long$comparison_type,
  levels = c('Within KB', 'Within KC', 'Between Regions')
)

# Summary table (prints in console)
lkefiri_region_summary <- lkefiri_dist_long |>
  group_by(comparison_type) |>
  summarise(
    n_comparisons = n(),
    mean_distance = mean(phylogenetic_distance),
    sd_distance   = sd(phylogenetic_distance),
    .groups = 'drop'
  )

print(lkefiri_region_summary)

# Optional: save the summary table
write_tsv(
  lkefiri_region_summary,
  'data/results/06_strain_profiling/lkefiri_region_structure_summary.tsv'
)

# Overall test (prints in console) - only if >1 group exists
if (n_distinct(lkefiri_dist_long$comparison_type) > 1) {
  
  lkefiri_region_kw <- kruskal.test(
    phylogenetic_distance ~ comparison_type,
    data = lkefiri_dist_long
  )
  
  print(lkefiri_region_kw)
  
  # Preparation for saving
  lkefiri_region_kw_df <- tibble(
    test = 'Kruskal-Wallis',
    chi_squared = lkefiri_region_kw$statistic,
    df = lkefiri_region_kw$parameter,
    p_value = lkefiri_region_kw$p.value
  )
  
  write_tsv(
    lkefiri_region_kw_df,
    'data/results/06_strain_profiling/lkefiri_region_kw_test.tsv'
  )
  
  # Format Kruskal–Wallis p-value
  lkefiri_kw_p <- signif(lkefiri_region_kw$p.value, 3)
  
  # Create plotmath label (italic p)
  lkefiri_kw_label <- paste0(
    "'Kruskal-Wallis: ' ~ italic(p) == ",
    lkefiri_kw_p
  )
  
} else {
  
  message('Only one region present. Kruskal-Wallis test not applicable.')
  lkefiri_kw_label <- NULL
}


# Boxplot of Distances by Comparison Class --------------------------------

lkefiri_region_plot <- ggplot(
  lkefiri_dist_long,
  aes(x = comparison_type, y = phylogenetic_distance, fill = comparison_type)
) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # Add Kruskal–Wallis annotation (only if applicable)
  {
    if (!is.null(lkefiri_kw_label)) {
      annotate(
        'text',
        x = 2,
        y = max(lkefiri_dist_long$phylogenetic_distance) * 1.05,
        label = lkefiri_kw_label,
        parse = TRUE,
        size = 4
      )
    }
  } +
  
  scale_fill_manual(
    values = c(
      'Within KB' = '#4FAFA8',
      'Within KC' = '#E8403D',
      'Between Regions' = 'grey60'
    )
  ) +
  
  labs(
    title = '(D)',
    subtitle = expression(bolditalic('Lentilactobacillus kefiri')),
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
lkefiri_region_plot

# Save the plot
ggsave(
  lkefiri_region_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_region/lkefiri_region_structure_boxplot.png',
  height = 15,
  width = 20,
  units = c('cm'),
  dpi = 300
)

# Combine the Plots of Four Dominant Bacteria -----------------------------

library(patchwork)

combined_phylogenetic_distance_plot <- (edurans_hm | llactis_hm) / (lkefiranofaciens_hm | lkefiri_hm)

# Check the plot
combined_phylogenetic_distance_plot

# Save the plot
ggsave(combined_phylogenetic_distance_plot,
       filename = 'plots/06_strain_profiling/6.1_strainphlan_annotation/combined_phylogenetic_distance_plot.png',
       height = 30,
       width = 40,
       units = c('cm'),
       dpi = 300
)