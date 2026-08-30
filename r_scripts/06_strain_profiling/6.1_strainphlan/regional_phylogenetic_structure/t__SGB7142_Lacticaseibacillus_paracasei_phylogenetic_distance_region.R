# StrainPhlAn t__SGB7142 Lacticaseibacillus paracasei Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/regional_phylogenetic_structure
# 17/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkarian' = '#4FAFA8', 'Karachay-Cherkess' = '#E8403D')


# Load the data -----------------------------------------------------------

lparacasei_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7142_Lacticaseibacillus_paracasei_phylogenetic_distance_new.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
lparacasei_mat_df <- lparacasei_data |>
  bind_rows(
    lparacasei_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
lparacasei_all_ids <- sort(unique(c(lparacasei_mat_df$strain_id_1, lparacasei_mat_df$strain_id_2)))

lparacasei_diag_df <- tibble(
  strain_id_1 = lparacasei_all_ids,
  strain_id_2 = lparacasei_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
lparacasei_mat_df <- bind_rows(lparacasei_diag_df, lparacasei_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)


# II. Within-Region vs Between-Region Structure ---------------------------

# Long-format distances for REGION comparisons (unique pairs only; remove self)
lparacasei_dist_long <- lparacasei_mat_df |>
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
lparacasei_dist_long$comparison_type <- factor(
  lparacasei_dist_long$comparison_type,
  levels = c('Within KB', 'Within KC', 'Between Regions')
)

# Summary table (prints in console)
lparacasei_region_summary <- lparacasei_dist_long |>
  group_by(comparison_type) |>
  summarise(
    n_comparisons = n(),
    mean_distance = mean(phylogenetic_distance),
    sd_distance   = sd(phylogenetic_distance),
    .groups = 'drop'
  )

print(lparacasei_region_summary)

# Optional: save the summary table
write_tsv(
  lparacasei_region_summary,
  'data/results/06_strain_profiling/lparacasei_region_structure_summary.tsv'
)

# Kruskal–Wallis only if >= 2 comparison groups exist
lparacasei_n_groups <- n_distinct(lparacasei_dist_long$comparison_type)

if (lparacasei_n_groups >= 2) {
  
  # Overall test (prints in console)
  lparacasei_region_kw <- kruskal.test(
    phylogenetic_distance ~ comparison_type,
    data = lparacasei_dist_long
  )
  
  print(lparacasei_region_kw)
  
  # Preparation for saving
  lparacasei_region_kw_df <- tibble(
    test = 'Kruskal-Wallis',
    chi_squared = as.numeric(lparacasei_region_kw$statistic),
    df = as.numeric(lparacasei_region_kw$parameter),
    p_value = as.numeric(lparacasei_region_kw$p.value)
  )
  
  write_tsv(
    lparacasei_region_kw_df,
    'data/results/06_strain_profiling/lparacasei_region_kw_test.tsv'
  )
  
  # Format Kruskal–Wallis p-value
  lparacasei_kw_p <- signif(lparacasei_region_kw$p.value, 3)
  
  # Create plotmath label (italic p)
  lparacasei_kw_label <- paste0(
    "'Kruskal-Wallis: ' ~ italic(p) == ",
    lparacasei_kw_p
  )
  
} else {
  
  message('Only one region present. Kruskal-Wallis test not applicable.')
  lparacasei_kw_label <- NULL
}


# Pairwise Comparison -----------------------------------------------------

lparacasei_pairwise <- pairwise.wilcox.test(
  lparacasei_dist_long$phylogenetic_distance,
  lparacasei_dist_long$comparison_type,
  p.adjust.method = 'BH'
)

# Save pairwise p-values
lparacasei_pairwise_df <- as.data.frame(as.table(lparacasei_pairwise$p.value)) |>
  drop_na() |>
  rename(
    group_1 = Var1,
    group_2 = Var2,
    p_value = Freq
  )

write_tsv(
  lparacasei_pairwise_df,
  'data/results/06_strain_profiling/lparacasei_pairwise_wilcox_BH.tsv'
)

# Convert pairwise results into "bracket" rows (only significant)
lparacasei_sig_df <- lparacasei_pairwise_df |>
  mutate(
    group_1 = as.character(group_1),
    group_2 = as.character(group_2),
    p_value = as.numeric(p_value),
    sig = case_when(
      p_value < 0.001 ~ '***',
      p_value < 0.01  ~ '**',
      p_value < 0.05  ~ '*',
      TRUE ~ 'ns'
    )
  ) |>
  filter(sig != 'ns') |>
  mutate(
    x1 = as.numeric(factor(group_2, levels = levels(lparacasei_dist_long$comparison_type))),
    x2 = as.numeric(factor(group_1, levels = levels(lparacasei_dist_long$comparison_type)))
  )


# Boxplot of Distances by Comparison Class --------------------------------

y_max <- max(lparacasei_dist_long$phylogenetic_distance)
base_y <- y_max * 1.05
step_y <- y_max * 0.08

lparacasei_region_plot <- ggplot(
  lparacasei_dist_long,
  aes(x = comparison_type, y = phylogenetic_distance, fill = comparison_type)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  # Pairwise significance brackets (only significant ones)
  {
    if (nrow(lparacasei_sig_df) > 0) {
      lapply(seq_len(nrow(lparacasei_sig_df)), function(i) {
        yi <- base_y + (i - 1) * step_y
        
        list(
          annotate(
            'segment',
            x = lparacasei_sig_df$x1[i], xend = lparacasei_sig_df$x2[i],
            y = yi, yend = yi
          ),
          annotate(
            'segment',
            x = lparacasei_sig_df$x1[i], xend = lparacasei_sig_df$x1[i],
            y = yi - step_y * 0.3, yend = yi
          ),
          annotate(
            'segment',
            x = lparacasei_sig_df$x2[i], xend = lparacasei_sig_df$x2[i],
            y = yi - step_y * 0.3, yend = yi
          ),
          annotate(
            'text',
            x = (lparacasei_sig_df$x1[i] + lparacasei_sig_df$x2[i]) / 2,
            y = yi + step_y * 0.15,
            label = lparacasei_sig_df$sig[i],
            size = 6
          )
        )
      })
    }
  } +
  
  # Add Kruskal–Wallis annotation (only if available)
  {
    if (!is.null(lparacasei_kw_label)) {
      annotate(
        'text',
        x = 2,
        y = y_max * 1.20,
        label = lparacasei_kw_label,
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
    title = '(C)',
    subtitle = expression(bolditalic('Lacticaseibacillus paracasei')),
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
lparacasei_region_plot

# Save the plot
ggsave(
  lparacasei_region_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_region/lparacasei_region_structure_boxplot.tiff',
  height = 15,
  width = 20,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine Plots of Bacteria for Regional Structuring ----------------------

library(gridExtra)

layout <- rbind(
  c(1, 2, 3)
)

combined_region_plot <- grid.arrange(
  edurans_region_plot, 
  llactis_region_plot,
  lparacasei_region_plot,
  layout_matrix = layout,
  widths = c(1, 1, 1)
)

# Save the plot
ggsave(
  combined_region_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_region/combined_region_plot_gridextra.tiff',
  height = 15,
  width = 30,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)