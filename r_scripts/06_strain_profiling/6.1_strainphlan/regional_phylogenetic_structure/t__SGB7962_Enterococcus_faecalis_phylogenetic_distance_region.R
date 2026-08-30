# StrainPhlAn t__SGB7962 Enterococcus faecalis Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/regional_phylogenetic_structure
# 17/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------
# NOTE: E. faecalis has ONLY KC data here (so everything is Karachay-Cherkess)

site_cols <- c('Kabardino-Balkarian' = '#4FAFA8', 'Karachay-Cherkess' = '#E8403D')


# Load the data -----------------------------------------------------------

efaecalis_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7962_Enterococcus_faecalis_phylogenetic_distance.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
efaecalis_mat_df <- efaecalis_data |>
  bind_rows(
    efaecalis_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
efaecalis_all_ids <- sort(unique(c(efaecalis_mat_df$strain_id_1, efaecalis_mat_df$strain_id_2)))

efaecalis_diag_df <- tibble(
  strain_id_1 = efaecalis_all_ids,
  strain_id_2 = efaecalis_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
efaecalis_mat_df <- bind_rows(efaecalis_diag_df, efaecalis_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)


# II. Within-Region vs Between-Region Structure ---------------------------

# Long-format distances for REGION comparisons (unique pairs only; remove self)
efaecalis_dist_long <- efaecalis_mat_df |>
  filter(strain_id_1 < strain_id_2) |>
  mutate(
    # Only KC present here
    region_1 = 'Karachay-Cherkess',
    region_2 = 'Karachay-Cherkess',
    comparison_type = 'Within KC'
  ) |>
  filter(phylogenetic_distance > 0)

# Order comparison classes (keep same global order for consistency)
efaecalis_dist_long$comparison_type <- factor(
  efaecalis_dist_long$comparison_type,
  levels = c('Within KB', 'Within KC', 'Between Regions')
)

# Summary table (prints in console)
efaecalis_region_summary <- efaecalis_dist_long |>
  group_by(comparison_type) |>
  summarise(
    n_comparisons = n(),
    mean_distance = mean(phylogenetic_distance),
    sd_distance   = sd(phylogenetic_distance),
    .groups = 'drop'
  )

print(efaecalis_region_summary)

# Optional: save the summary table
write_tsv(
  efaecalis_region_summary,
  'data/results/06_strain_profiling/efaecalis_region_structure_summary.tsv'
)

# Kruskal–Wallis only if >=2 groups are present
efaecalis_kw_label <- NULL

efaecalis_n_groups <- efaecalis_dist_long |>
  distinct(comparison_type) |>
  nrow()

if (efaecalis_n_groups >= 2) {
  
  efaecalis_region_kw <- kruskal.test(
    phylogenetic_distance ~ comparison_type,
    data = efaecalis_dist_long
  )
  
  print(efaecalis_region_kw)
  
  # Preparation for saving
  efaecalis_region_kw_df <- tibble(
    test = 'Kruskal-Wallis',
    chi_squared = efaecalis_region_kw$statistic,
    df = efaecalis_region_kw$parameter,
    p_value = efaecalis_region_kw$p.value
  )
  
  write_tsv(
    efaecalis_region_kw_df,
    'data/results/06_strain_profiling/efaecalis_region_kw_test.tsv'
  )
  
  # Format Kruskal–Wallis p-value
  efaecalis_kw_p <- signif(efaecalis_region_kw$p.value, 3)
  
  # Create plotmath label (italic p)
  efaecalis_kw_label <- paste0(
    "'Kruskal-Wallis: ' ~ italic(p) == ",
    efaecalis_kw_p
  )
  
} else {
  
  message('Only one region present. Kruskal-Wallis test not applicable.')
  efaecalis_kw_label <- "'Only one region present; Kruskal-Wallis not applicable'"
}


# Boxplot of Distances by Comparison Class --------------------------------

efaecalis_region_plot <- ggplot(
  efaecalis_dist_long,
  aes(x = comparison_type, y = phylogenetic_distance, fill = comparison_type)
) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # Add annotation (KW if possible, otherwise message)
  annotate(
    'text',
    x = 1,
    y = max(efaecalis_dist_long$phylogenetic_distance) * 1.08,
    label = efaecalis_kw_label,
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
    title = '(E)',
    subtitle = expression(bolditalic('Enterococcus faecalis')),
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
efaecalis_region_plot

# Save the plot
ggsave(
  efaecalis_region_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_region/efaecalis_region_structure_boxplot.png',
  height = 15,
  width = 20,
  units = c('cm'),
  dpi = 300
)