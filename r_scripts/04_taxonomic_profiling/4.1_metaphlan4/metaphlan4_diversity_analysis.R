# Metaphlan4 Diversity Analysis
# ~/Documents/RStudio/mk_caucasus/scripts/04_taxonomic_profiling/4.1_metaphlan4
# 06/04/2025

# Load the necessary packages.
library(tidyverse)


# Load and Check Data -----------------------------------------------------

# Load the data
mp4_diversity_species_data <- read.table('data/04_taxonomic_profiling/metaphlan4/mp4_caucasus_merged_species.tsv',
                               header = TRUE, row.names = 1, sep = '\t')

# Transpose so rows = samples, cols = taxa
mp4_species_data_tp <- t(mp4_diversity_species_data) # tp: transpose


# Alpha Diversity ---------------------------------------------------------

library(vegan)

# Shannon diversity
mp4_shannon <- diversity(mp4_species_data_tp, index = 'shannon')

# Simpson diversity
mp4_simpson <- diversity(mp4_species_data_tp, index = 'simpson')

# Combine results
mp4_alpha_df <- data.frame(
  sample = rownames(mp4_species_data_tp),
  shannon = mp4_shannon,
  simpson = mp4_simpson
)

# Create a column for sample type and collection site
mp4_alpha_df <- mp4_alpha_df |> 
  mutate(
    sample_type = ifelse(grepl('_g$', mp4_alpha_df$sample), 'Kefir Grains', 'Milk Kefir'),
    collection_site = case_when(
      sample %in% c('k1_m', 'k2_m', 'k3_m', 'k5_m', 'k6_m', 'k7_m', 'k8_m', 'k9_m',
                    'k1_g', 'k3_g', 'k5_g', 'k6_g', 'k7_g', 'k9_g') ~ 'KB',
      TRUE ~ 'KC'
    )
  )

# Use only Shannon and Simpson Metrics for Alpha Diversity
mp4_alpha_df_long <- mp4_alpha_df |> 
  select(sample, sample_type, collection_site, shannon, simpson) |> 
  pivot_longer(cols = c(shannon, simpson),
               names_to = 'metric',
               values_to = 'value') |> 
  mutate(metric = factor(metric, levels = c('shannon', 'simpson')))

# Filter the dataframe to include only the grains and their corresponding milks 
mp4_alpha_kgmk <- mp4_alpha_df_long |> 
  filter(sample %in% c('k1_g', 'k3_g', 'k5_g', 'k6_g', 'k7_g', 'k9_g',
                       'k1_m', 'k3_m', 'k5_m', 'k6_m', 'k7_m', 'k9_m'))

# Filter the dataframe to include only milk kefir samples (for diversity comparison later)
mp4_alpha_milk_kefir <- mp4_alpha_df_long %>%
  filter(sample_type == 'Milk Kefir')


## Sample Type ------------------------------------------------------------

library(ggpubr)

# Comparison between kefir grains and their corresponding milks
comparison_a <- list(c('Kefir Grains', 'Milk Kefir'))  

# Plot the data
mp4_alpha_plot_st <- ggplot(mp4_alpha_kgmk, aes(x = sample_type, y = value, fill = sample_type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
  stat_compare_means(comparisons = comparison_a, method = 'wilcox.test', label = 'p.signif', size = 5) + # disables the "exact test" attempt
  facet_wrap(~ metric, scales = 'free_y',
             labeller = as_labeller(c(shannon = 'Shannon', simpson = 'Simpson'))) + # strip text name
  # scale_fill_manual(values = c('kefir grain' = '#7b9acc', 'milk kefir' = '#FCF6F5')) + # specify box colors
  scale_fill_manual(values = c('Kefir Grains' = '#624E88', 'Milk Kefir' = '#FAB12F')) + # specify box colors
  labs(
    # title = 'Alpha Diversity of Caucasus Kefir Grains and their Corresponding Milk (MetaPhlAn4)',
    # title = 'Taxonomic',
    subtitle = '(B)',
    x = NULL,
    y = 'Diversity Index',
    fill = 'Sample Type') +
  theme_bw() +
  theme(
    # plot.title = element_text(face = 'bold', size = 15, hjust = 0.5), # center the plot title
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
    legend.key.size = unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Check the plot
mp4_alpha_plot_st

# Save the plot
ggsave(
  mp4_alpha_plot_st,
  filename = 'plots/04_taxonomic_profiling/metaphlan4/mp4_alpha_kgmk_for_combining_no_richness.tiff',
  height = 15,
  width = 20,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


## Collection Site --------------------------------------------------------

comparison_b <- list(c('KB', 'KC'))  


## Milk Kefir by Collection Site ------------------------------------------

# Plot the data
mp4_alpha_plot_mkcs <- ggplot(mp4_alpha_milk_kefir, aes(x = collection_site, y = value, fill = collection_site)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
  stat_compare_means(comparisons = comparison_b, method = 'wilcox.test', label = 'p.signif', size = 5) +
  facet_wrap(~ metric, scales = 'free_y',
             labeller = as_labeller(c(shannon = 'Shannon', simpson = 'Simpson'))) + # strip text name
  # scale_fill_manual(values = c('Kabardino-Balkarian' = '#7b9acc', 'Karachay-Cherkess' = '#FCF6F5')) + # specify box colors
  scale_fill_manual(values = c('KB' = '#4FACAC', 'KC' = '#E8403D')) + # specify box colors
  labs(
    # title = 'Alpha Diversity of Caucasus Milk Kefirs by Collection Site (MetaPhlAn4)',
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
    legend.key.size = unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Check the plot
mp4_alpha_plot_mkcs

# Save the plot
ggsave(
  mp4_alpha_plot_mkcs,
  filename = 'plots/04_taxonomic_profiling/metaphlan4/mp4_alpha_mkcs_for_combining_no_richness.tiff',
  height = 15,
  width = 20,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


## Statistical Table - Pairwise Wilcoxon Test -----------------------------

library(rstatix)

# Kefir grain vs Milk kefir (kgmk)
mp4_alpha_pairwise_kgmk <- mp4_alpha_kgmk |> 
  group_by(metric) |> 
  pairwise_wilcox_test(value ~ sample_type, p.adjust.method = 'none')

# View the table
print(mp4_alpha_pairwise_kgmk)

# Save the table
write_csv(mp4_alpha_pairwise_kgmk, 'data/results/04_taxonomic_profiling/metaphlan4/diversity_analysis/mp4_alpha_pairwise_kgmk.csv')

# Milk kefir by collection site

mp4_alpha_pairwise_mkcs <- mp4_alpha_milk_kefir |> 
  group_by(metric) |> 
  pairwise_wilcox_test(value ~ collection_site, p.adjust.method = 'none')

# View the table
print(mp4_alpha_pairwise_mkcs)

# Save the table
write_csv(mp4_alpha_pairwise_mkcs, 'data/results/04_taxonomic_profiling/metaphlan4/diversity_analysis/mp4_alpha_pairwise_mkcs.csv')


# Beta Diversity and Ordination -------------------------------------------

# Bray-Curtis dissimilarity
mp4_bray_dist <- vegdist(mp4_species_data_tp, method = 'bray')


# Interaction between Sample Type and Collection Site ---------------------

# Create metadata with both sample_type and collection_site
mp4_beta_df <- data.frame(sample = rownames(mp4_species_data_tp)) |>
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
mp4_pcoa_group <- ape::pcoa(mp4_bray_dist)

  # Extract coordinates and convert to dataframe
  mp4_pcoa_df_coord <- as.data.frame(mp4_pcoa_group$vectors) # converting to data frame
  mp4_pcoa_df_coord$sample <- rownames(mp4_pcoa_df_coord) # adding sample column
  
  # Merge mp4_pcoa_df_group with metadata from mp4_beta_df
  mp4_pcoa_df_group <- left_join(mp4_pcoa_df_coord, mp4_beta_df, by = 'sample')
  mp4_pcoa_df_group$shape <- ifelse(grepl('_g$', mp4_pcoa_df_group$sample), 'kefir grain', 'milk kefir')

  # How to check variance explained by each axis
  round(mp4_pcoa_group$values$Relative_eig[1:5] * 100, 2)
  # Returns 54.32 38.75  3.55  2.82  2.35
  # Axis 1 explains 54.95% of variation, Axis 2 explains 38.53%
  # Plotting just axis 1 & 2 captures ~94% - often enough to reveal strong structure

# Create the label text manually based on PERMANOVA result
group_pcoa_perm_label <- data.frame(
  Axis.1 = min(mp4_pcoa_df_group$Axis.1) - 0.5,
  Axis.2 = max(mp4_pcoa_df_group$Axis.2) + 0.55,
  label = "bold('PERMANOVA')"
)

group_pcoa_r2p_label <- data.frame(
  Axis.1 = min(mp4_pcoa_df_group$Axis.1) - 0.5,
  Axis.2 = max(mp4_pcoa_df_group$Axis.2) + 0.49,
  label = "italic(R^2)*' = 0.70, '*italic(p)*' = 0.001'"
)
  
# Plot the data
mp4_group_pcoa_plot <- ggplot(mp4_pcoa_df_group, aes(x = Axis.1, y = Axis.2)) +
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
    # title = 'Bray-Curtis - PCoA Ordination by Sample Type and Collection Site (MetaPhlAn4)', 
    title = '(D)',
    x = 'PCoA 1 (54.32%)',
    y = 'PCoA 2 (38.75%)',
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
    legend.key.size = unit(0.8, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Check the plot
mp4_group_pcoa_plot

# Save the plot
ggsave(
  mp4_group_pcoa_plot,
  filename = 'plots/04_taxonomic_profiling/metaphlan4/mp4_pcoa_group_for_combining.tiff',
  height = 20,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


## PERMANOVA --------------------------------------------------------------

# Combined model
mp4_permanova_group <-  adonis2(mp4_bray_dist ~ sample_type * collection_site, data = mp4_beta_df)

# Write to CSV
write.csv(as.data.frame(mp4_permanova_group),
          'data/results/04_taxonomic_profiling/metaphlan4/diversity_analysis/mp4_permanova_group.csv', row.names = TRUE)

# Marginal PERMANOVA
mp4_permanova_margin <- adonis2(mp4_species_data_tp ~ sample_type + collection_site, 
                                data = mp4_beta_df,
                                method = 'bray',
                                by = 'margin')

# The same as above but distance is pre-calculated in mp4_bray_dist
mp4_permanova_margin_v2 <- adonis2(mp4_bray_dist ~ sample_type + collection_site, 
                                   data = mp4_beta_df,
                                   by = 'margin')

# Write to CSV
write.csv(mp4_permanova_margin_v2,
          'data/results/04_taxonomic_profiling/metaphlan4/diversity_analysis/mp4_permanova_margin.csv', row.names = TRUE)


# Combine All Taxonomic Information ---------------------------------------

library(gridExtra)

# Layout
layout <- rbind(
  c(1, 1, 1), # top plot spans 3 columns
  c(2, 3, 4) # bottom plot spans 3 columns
)

# Plot the data
taxonomic_combined <- grid.arrange(
  mp4_species_plot,
  mp4_alpha_plot_st,
  mp4_alpha_plot_mkcs,
  mp4_group_pcoa_plot,
  layout_matrix = layout,
  heights = c(1.2, 1),
  widths = c(1.3, 1.3, 2)
)

# Save the plot
ggsave(
  taxonomic_combined,
  filename = 'plots/04_taxonomic_profiling/metaphlan4/taxonomic_combined_plot.tiff',
  height = 45,
  width = 50,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)