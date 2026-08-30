# HUMAnN3 Diversity Analysis without Experimental Kefir Samples
# ~/Documents/RStudio/mk_caucasus/scripts/04_functional_profiling/5.1_humann3
# 14/07/2025

# Load the necessary packages
library(tidyverse)


# Load and Check Data -----------------------------------------------------

hm3_species_data <- read_tsv(
  'data/05_functional_profiling/humann3_latest/caucasus_stratified_pathabundance_cpm/hm3_caucasus_joined_pathabundance_cpm_unstratified.tsv', 
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_Abundance-CPM$')) |> # remove suffix '_Abundance-CPM$' in all column names
  filter(!pathway %in% c('UNMAPPED', 'UNINTEGRATED')) # removing UNMAPPED and UNINTEGRATED rows in pathway column

# Set pathway column as rownames and transpose
hm3_species_data_tp <- hm3_species_data |>
  column_to_rownames(var = "pathway") |>  # Make pathway names the rownames
  t() |>  # Transpose so rows = samples, cols = pathways
  as.data.frame()


# Alpha Diversity ---------------------------------------------------------

library(vegan)

# Shannon diversity
hm3_shannon <- diversity(hm3_species_data_tp, index = 'shannon')

# Simpson diversity
hm3_simpson <- diversity(hm3_species_data_tp, index = 'simpson')

# Combine results
hm3_alpha_df <- data.frame(
  sample = rownames(hm3_species_data_tp),
  shannon = hm3_shannon,
  simpson = hm3_simpson
)

# Create a column for sample type and collection site
hm3_alpha_df <- hm3_alpha_df |> 
  mutate(
    sample_type = ifelse(grepl('_g$', hm3_alpha_df$sample), 'Kefir Grains', 'Milk Kefir'),
    collection_site = case_when(
      sample %in% c('k1_m', 'k2_m', 'k3_m', 'k5_m', 'k6_m', 'k7_m', 'k8_m', 'k9_m',
                    'k1_g', 'k3_g', 'k5_g', 'k6_g', 'k7_g', 'k9_g') ~ 'KB',
      TRUE ~ 'KC'
    )
  )

# Use only Shannon and Simpson Metrics for Alpha Diversity
hm3_alpha_df_long <- hm3_alpha_df |> 
  select(sample, sample_type, collection_site, shannon, simpson) |> 
  pivot_longer(cols = c(shannon, simpson),
               names_to = 'metric',
               values_to = 'value') |> 
  mutate(metric = factor(metric, levels = c('shannon', 'simpson')))

# Filter the dataframe to include only the grains and their corresponding milks 
hm3_alpha_kgmk <- hm3_alpha_df_long |> 
  filter(sample %in% c('k1_g', 'k3_g', 'k5_g', 'k6_g', 'k7_g', 'k9_g',
                       'k1_m', 'k3_m', 'k5_m', 'k6_m', 'k7_m', 'k9_m'))

# Filter the dataframe to include only milk kefir samples
hm3_alpha_milk_kefir <- hm3_alpha_df_long |>
  filter(sample_type == 'Milk Kefir')


## Sample Type ------------------------------------------------------------

library(ggpubr)

# Comparison between kefir grains and their corresponding milks
comparison_a <- list(c('Kefir Grains', 'Milk Kefir')) 

# Plot the data
hm3_alpha_plot_st <- ggplot(hm3_alpha_kgmk, aes(x = sample_type, y = value, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
  stat_compare_means(comparisons = comparison_a, method = 'wilcox.test', label = 'p.signif', size = 5) +
  facet_wrap(~ metric, scales = 'free_y',
             labeller = as_labeller(c(shannon = 'Shannon', simpson = 'Simpson'))) +
  scale_fill_manual(values = c('Kefir Grains' = '#624E88', 'Milk Kefir' = '#FAB12F')) +
  labs(
    # title = 'Alpha Diversity of Caucasus Kefir Grains and their Corresponding Milk (HUMAnN3)',
    # title = 'Functional',
    subtitle = '(B)',
    x = NULL,
    y = 'Diversity Index',
    fill = 'Sample Type') +
  theme_bw() +
  theme(
    # plot.title = element_text(face = 'bold', size = 15, hjust = 0.5),
    plot.subtitle = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    legend.position = 'bottom',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Check the plot
hm3_alpha_plot_st

# Save the plot
ggsave(
  hm3_alpha_plot_st,
  filename = 'plots/05_functional_profiling/humann3_latest/hm3_alpha_kgmk_for_combining_no_richness.tiff',
  height = 15,
  width = 20,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


## Collection Site --------------------------------------------------------

comparison_b <- list(c('KB', 'KC')) 


# Milk Kefir by Collection Site -------------------------------------------

hm3_alpha_plot_mkcs <- ggplot(hm3_alpha_milk_kefir, aes(x = collection_site, y = value, fill = collection_site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
  stat_compare_means(comparisons = comparison_b, method = 'wilcox.test', label = 'p.signif', size = 5) +
  facet_wrap(~ metric, scales = 'free_y',
             labeller = as_labeller(c(shannon = 'Shannon', simpson = 'Simpson'))) +
  scale_fill_manual(values = c('KB' = '#4FACAC', 'KC' = '#E8403D')) +
  labs(
    # title = 'Alpha Diversity of Caucasus Milk Kefirs by Collection Site (HUMAnN3)',
    title = '(C)',
    x = NULL,
    y = 'Diversity Index',
    fill = 'Caucasus Region (milk kefir only)') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    legend.position = 'bottom',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Check the plot
hm3_alpha_plot_mkcs

# Save the plot
ggsave(
  hm3_alpha_plot_mkcs,
  filename = 'plots/05_functional_profiling/humann3_latest/hm3_alpha_mkcs_for_combining_no_richness.tiff',
  height = 15,
  width = 20,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


## Statistical Table - Pairwise Wilcoxon Test -----------------------------

library(rstatix)

# Kefir grain vs Milk kefir (kgmk)
hm3_alpha_pairwise_kgmk <- hm3_alpha_kgmk |> 
  group_by(metric) |> 
  pairwise_wilcox_test(value ~ sample_type, p.adjust.method = 'none')

# View the table
print(hm3_alpha_pairwise_kgmk)

# Save the table
write_csv(hm3_alpha_pairwise_kgmk, 'data/results/05_functional_profiling/humann3_latest/diversity_analysis/hm3_alpha_pairwise_kgmk.csv')

# Milk kefir by collection site

hm3_alpha_pairwise_mkcs <- hm3_alpha_milk_kefir |> 
  group_by(metric) |> 
  pairwise_wilcox_test(value ~ collection_site, p.adjust.method = 'none')

# View the table
print(hm3_alpha_pairwise_mkcs)

# Save the table
write_csv(hm3_alpha_pairwise_mkcs, 'data/results/05_functional_profiling/humann3_latest/diversity_analysis/hm3_alpha_pairwise_mkcs.csv')


# Beta Diversity and Ordination -------------------------------------------

# Bray-Curtis dissimilarity
hm3_bray_dist <- vegdist(hm3_species_data_tp, method = 'bray')


# Interaction between Sample Type and Collection Site ---------------------

# Create metadata with both sample_type and collection_site
hm3_beta_df <- data.frame(sample = rownames(hm3_species_data_tp)) |>
  mutate(
    sample_type = ifelse(grepl('_g$', sample), 'kefir grain', 'milk kefir'),
    collection_site = case_when(
      sample %in% c('k1_m', 'k2_m', 'k3_m', 'k5_m', 'k6_m', 'k7_m', 'k8_m', 'k9_m',
                    'k1_g', 'k3_g', 'k5_g', 'k6_g', 'k7_g', 'k9_g') ~ 'KB',
      TRUE ~ 'KC'
    ),
    group = interaction(sample_type, collection_site)
  )


## PCoA Ordination --------------------------------------------------------

library(ape)

# Principal Coordinates Analysis (PCoA)
hm3_pcoa_group <- ape::pcoa(hm3_bray_dist)

  # Extract coordinates and convert to dataframe
  hm3_pcoa_df_coord <- as.data.frame(hm3_pcoa_group$vectors) # converting to data frame
  hm3_pcoa_df_coord$sample <- rownames(hm3_pcoa_df_coord) # adding sample column
  
  # Merge hm3_pcoa_df_group with metadata from hm3_beta_df
  hm3_pcoa_df_group <- left_join(hm3_pcoa_df_coord, hm3_beta_df, by = 'sample')
  hm3_pcoa_df_group$shape <- ifelse(grepl('_g$', hm3_pcoa_df_group$sample), 'kefir grain', 'milk kefir')

  # How to check variance explained by each axis
  round(hm3_pcoa_group$values$Relative_eig[1:5] * 100, 2)
  # Returns 56.08 25.46  9.17  2.77  2.18
  # Axis 1 explains 56.08% of variation, Axis 2 explains 25.46%
  # Plotting just axis 1 & 2 captures ~82% - often enough to reveal strong structure

# Create the label text manually based on PERMANOVA result
group_pcoa_perm_label <- data.frame(
  Axis.1 = min(hm3_pcoa_df_group$Axis.1) - 0.25,
  Axis.2 = max(hm3_pcoa_df_group$Axis.2) - 0.60,
  label = "bold('PERMANOVA')"
)

group_pcoa_r2p_label <- data.frame(
  Axis.1 = min(hm3_pcoa_df_group$Axis.1) - 0.25,
  Axis.2 = max(hm3_pcoa_df_group$Axis.2) - 0.63,
  label = "italic(R^2)*' = 0.58, '*italic(p)*' = 0.001'"
)

# Plot the data
hm3_group_pcoa_plot <- ggplot(hm3_pcoa_df_group, aes(x = Axis.1, y = Axis.2)) +
  stat_ellipse(geom = 'polygon', aes(group = group, fill = group, color = group), 
               type = 'norm', level = 0.95, linetype = 'dashed', alpha = 0.5) +
  geom_point(aes(color = group, shape = shape), size = 3) +
  scale_shape_manual(values = c('kefir grain' = 17, 'milk kefir' = 15),
                     labels = c('kefir grain' = 'Kefir Grain', 'milk kefir' = 'Milk Kefir')) +
  scale_fill_manual(values = c('kefir grain.KB' = '#F5CBCB',
                               'milk kefir.KB' = '#BBDCE5',
                               'kefir grain.KC' = '#D4B996',
                               'milk kefir.KC' = '#CADCAE'),
                    labels = c('kefir grain.KB' = 'Kefir Grain - KB',
                               'milk kefir.KB' = 'Milk Kefir - KB',
                               'kefir grain.KC' = 'Kefir Grain - KC',
                               'milk kefir.KC' = 'Milk Kefir - KC')) +
  scale_color_manual(values = c('kefir grain.KB' = '#A4193D',
                                'milk kefir.KB' = '#0063B2',
                                'kefir grain.KC' = '#A07855',
                                'milk kefir.KC' = '#2C5F2D'),
                     labels = c('kefir grain.KB' = 'Kefir Grain - KB',
                                'milk kefir.KB' = 'Milk Kefir - KB',
                                'kefir grain.KC' = 'Kefir Grain - KC',
                                'milk kefir.KC' = 'Milk Kefir - KC')) +
  labs(
    # title = 'Bray-Curtis - PCoA Ordination by Sample Type and Collection Site (HUMAnN3)', 
    title = '(D)',
    x = 'PCoA 1 (56.08%)', 
    y = 'PCoA 2 (25.46%)',
    color = 'Group',
    shape = 'Sample Type',
    fill = 'Group') +
  guides(
    shape = guide_legend(order = 1),  # Sample Type on top
    color = guide_legend(order = 2),  # Group below
    fill  = guide_legend(order = 2)
  ) +
  geom_text(data = group_pcoa_perm_label, aes(x = Axis.1, y = Axis.2, label = label),
            parse = TRUE, hjust = 0, vjust = 1, size = 5) +
  geom_text(data = group_pcoa_r2p_label, aes(x = Axis.1, y = Axis.2, label = label),
            parse = TRUE, hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    legend.box = 'vertical',
    legend.box.spacing = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.position = 'bottom', # move the legends at the bottom
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Check the plot
hm3_group_pcoa_plot

# Save the plot
ggsave(
  hm3_group_pcoa_plot,
  filename = 'plots/05_functional_profiling/humann3_latest/hm3_pcoa_group_for_combining.tiff',
  height = 20,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


## PERMANOVA --------------------------------------------------------------

# Combined model
hm3_permanova_group <-  adonis2(hm3_bray_dist ~ sample_type * collection_site, data = hm3_beta_df)

# Write to CSV
write.csv(as.data.frame(hm3_permanova_group),'data/results/05_functional_profiling/humann3_latest/diversity_analysis/hm3_permanova_group.csv', row.names = TRUE)

# Marginal PERMANOVA
hm3_permanova_margin <- adonis2(hm3_species_data_tp ~ sample_type + collection_site, 
                                data = hm3_beta_df,
                                method = 'bray',
                                by = 'margin')

# The same as above but distance is pre-calculated in hm3_bray_dist
hm3_permanova_margin_v2 <- adonis2(hm3_bray_dist ~ sample_type + collection_site, 
                                   data = hm3_beta_df,
                                   by = 'margin')

# Write to CSV
write.csv(hm3_permanova_margin_v2, 'data/results/05_functional_profiling/humann3_latest/diversity_analysis/hm3_permanova_margin.csv', row.names = TRUE)


# Combine All Functional Information --------------------------------------

library(gridExtra)

# Option 1
layout <- rbind(
  c(1, 1, 1), # top plot spans 3 columns
  c(2, 3, 4) # bottom plot spans 3 columns
)

functional_combined <- grid.arrange(
  hm3_strat_contributor_plot,
  hm3_alpha_plot_st,
  hm3_alpha_plot_mkcs,
  hm3_group_pcoa_plot,
  layout_matrix = layout,
  heights = c(1.2, 1),
  widths = c(1.3, 1.3, 2)
)

# Save the plot
ggsave(
  functional_combined,
  filename = 'plots/05_functional_profiling/humann3_latest/functional_combined_plot.tiff',
  height = 45,
  width = 50,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)