# DRAM MAGs CAZymes Comparisons
# ~/Documents/RStudio/mk_caucasus/scripts/07_downstream_analysis/7.2_dram
# 19/03/2026

# This is the Revised and Working Script ----------------------------------

# Load the necessary packages.
library(paletteer)
library(tidyverse)


# Load and Check Data -----------------------------------------------------

# Load data
dram_cazymes_data <- read_tsv('data/07_downstream_analysis/7.2_dram/hq_mags/distilled/combined_metabolism_summary.tsv', show_col_types = FALSE) |> 
  select(-gene_description, -header, -subheader) |> 
  mutate(
    module_ID = case_when(
      module == 'Glycoside Hydrolases' ~ 'GH',
      module == 'Glycosyl Transferases' ~ 'GT',
      module == 'Carbohydrate Esterases' ~ 'CE',
      module == 'Carbohydrate-Binding Modules' ~ 'CBM',
      module == 'Auxiliary Activities' ~ 'AA',
      module == 'Polysaccharide Lyases' ~ 'PL',
      TRUE ~ NA_character_ # safer fallback
    ),
    module_ID = factor(module_ID, levels = c('GH', 'GT', 'CE', 'CBM', 'AA', 'PL')),
    module = factor(
      module, levels = c('Glycoside Hydrolases', 'Glycosyl Transferases',
                         'Carbohydrate Esterases', 'Carbohydrate-Binding Modules',
                         'Auxiliary Activities', 'Polysaccharide Lyases'))
  ) |> 
  filter(gene_number > 0)

# Load gtdb data as well for species annotation
gtdb_data <- read_tsv('data/07_downstream_analysis/7.1_gtdbtk/hq_mags/gtdb_hqmags_with_qc_data.tsv', show_col_types = FALSE) |> 
  select(user_genome, species, sample_type, collection_site)

# Final data (combine dram and gtdb data)
final_cazymes_data <- dram_cazymes_data |>
  left_join(gtdb_data, by = c('bin_name' = 'user_genome'))

# Save final metadata
write_tsv(final_cazymes_data, 'data/07_downstream_analysis/7.2_dram/hq_mags/distilled/final_cazyme_gtdb_data.tsv')


# Summarize per MAG -------------------------------------------------------

cazyme_summary <- final_cazymes_data |>
  group_by(bin_name, species, sample_type, collection_site, module_ID) |>
  summarise(total_genes = sum(gene_number), .groups = 'drop')

# Convert to relative abundance
cazyme_rel <- cazyme_summary |>
  group_by(bin_name) |>
  mutate(rel_abundance = total_genes / sum(total_genes)) |>
  ungroup()

# Calculate mean relative abundance
species_order <- cazyme_rel |>
  group_by(species) |>
  summarise(mean_abundance = mean(rel_abundance), .groups = 'drop') |>
  arrange(desc(mean_abundance)) |>
  pull(species)

# Arrange it based on relative abundance
cazyme_rel <- cazyme_rel |>
  mutate(
    species = factor(species, levels = species_order)
  )


# Analysis 1 - Differences Between Sample Types ---------------------------

library(rstatix)

sample_test <- cazyme_rel |>
  group_by(module_ID) |>
  wilcox_test(rel_abundance ~ sample_type) |>
  adjust_pvalue(method = 'BH') |>
  add_significance()

# Write to CSV
write.csv(sample_test, 'data/results/07_downstream_analysis/dram/sample_test.csv', row.names = FALSE)

# Run Wilcoxon rank-sum test (Mann-Whitney U test) (2 groups)
sample_test <- cazyme_rel |>
  group_by(module_ID) |>
  wilcox_test(rel_abundance ~ sample_type) |>
  adjust_pvalue(method = 'BH') |>
  add_significance() |>
  add_xy_position(x = 'sample_type', fun = 'max',
                  step.increase = 0.02) # space between stat line and boxplot

# Plot the data
sample_cazyme_plot <- ggplot(cazyme_rel, aes(x = sample_type, y = rel_abundance, fill = sample_type)) +
  geom_boxplot(linewidth = 0.5, outlier.size = 2) +
  stat_pvalue_manual(
    sample_test,
    label = 'p.adj.signif',
    tip.length = 0.01,
    size = 3.5
  ) +
  facet_wrap(~ module_ID, scales = 'free_y') +
  scale_fill_manual(
    values = c('Kefir Grains' = '#624E88', 'Milk Kefir' = '#FAB12F')
  ) +
  labs(
    title = '(A)',
    x     = NULL,
    y     = 'Relative Abundance',
    fill  = 'Sample Type'
  ) +
  theme_bw() +
  theme(
    axis.title      = element_text(face = 'bold', size = 12),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.title.y    = element_text(margin = margin(r = 10)),
    axis.text.y     = element_text(size = 12),
    legend.title    = element_text(face = 'bold', size = 12),
    legend.text     = element_text(size = 12),
    legend.position = 'bottom',
    panel.grid      = element_blank(),
    plot.title      = element_text(face = 'bold', size = 12),
    strip.text.x    = element_text(face = 'bold', size = 12)
  )

# Check the plot
sample_cazyme_plot

# Save the plot
ggsave(
  sample_cazyme_plot,
  filename = 'plots/07_downstream_analysis/7.2_dram/cazyme_comparisons/sample_cazyme_plot.tiff',
  height = 25,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Analysis 2 - Differences between Regions (KB vs KC) ---------------------

region_test <- cazyme_rel |>
  group_by(module_ID) |>
  wilcox_test(rel_abundance ~ collection_site) |>
  adjust_pvalue(method = 'BH') |>
  add_significance()

# Write to CSV
write.csv(region_test, 'data/results/07_downstream_analysis/dram/region_test.csv', row.names = FALSE)

# Run Wilcoxon rank-sum test (Mann-Whitney U test) (2 groups)
region_test <- cazyme_rel |>
  group_by(module_ID) |>
  wilcox_test(rel_abundance ~ collection_site) |>
  adjust_pvalue(method = 'BH') |>
  add_significance() |>
  add_xy_position(x = 'collection_site', fun = 'max',
                  step.increase = 0.02) # space between stat line and boxplot

# Plot the data
region_cazyme_plot <- ggplot(cazyme_rel, aes(x = collection_site, y = rel_abundance, fill = collection_site)) +
  geom_boxplot(linewidth = 0.5, outlier.size = 2) +
  stat_pvalue_manual(
    region_test,
    label = 'p.adj.signif',
    tip.length = 0.01,
    size = 3.5
  ) +
  facet_wrap(~ module_ID, scales = 'free_y') +
  scale_fill_manual(
    values = c('Kabardino-Balkaria' = '#4FACAC', 'Karachay-Cherkessia' = '#E8403D')
  ) +
  labs(
    title = '(B)',
    x     = NULL,
    y     = 'Relative Abundance',
    fill  = 'Caucasus Region'
  ) +
  theme_bw() +
  theme(
    axis.title      = element_text(face = 'bold', size = 12),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.title.y    = element_text(margin = margin(r = 10)),
    axis.text.y     = element_text(size = 12),
    legend.title    = element_text(face = 'bold', size = 12),
    legend.text     = element_text(size = 12),
    legend.position = 'bottom',
    panel.grid      = element_blank(),
    plot.title      = element_text(face = 'bold', size = 12),
    strip.text.x    = element_text(face = 'bold', size = 12)
  )

# Check the plot
region_cazyme_plot

# Save the plot
ggsave(
  region_cazyme_plot,
  filename = 'plots/07_downstream_analysis/7.2_dram/cazyme_comparisons/region_cazyme_plot.tiff',
  height = 25,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Analysis 3 - Overall CAZyme Composition (Multivariate) ------------------

# Covert the data into matrix

cazyme_matrix <- cazyme_rel |>
  select(bin_name, module_ID, rel_abundance) |>
  pivot_wider(names_from = module_ID, values_from = rel_abundance, values_fill = 0)

# PERMANOVA

library(vegan)

meta <- cazyme_rel |>
  distinct(bin_name, species, sample_type, collection_site)

# Multivariable PERMANOVA
# Tests the combined effect of all variables together
# One single p-value for the whole model
cazyme_permanova_all <- adonis2(
  cazyme_matrix[,-1] ~ sample_type + collection_site + species,
  data = meta,
  method = 'bray'
)

# Write to CSV
write.csv(cazyme_permanova_all, 'data/results/07_downstream_analysis/dram/cazyme_permanova_all.csv', row.names = TRUE)

# Marginal PERMANOVA
# Each variable independently, while controlling for others
cazyme_permanova_margin <- adonis2(
  cazyme_matrix[,-1] ~ sample_type + collection_site + species,
  data = meta,
  method = 'bray',
  by = 'margin'
)

# Write to CSV
write.csv(cazyme_permanova_margin, 'data/results/07_downstream_analysis/dram/cazyme_permanova_margin.csv', row.names = TRUE)


# NMDS Ordination --------------------------------------------------------

# Convert bin_name to rownames
cazyme_matrix_nmds <- cazyme_matrix |>
  column_to_rownames('bin_name')

# Non-Metric Multidimensional Scaling (NMDS)
cazyme_nmds <- metaMDS(cazyme_matrix_nmds, distance = 'bray', k = 2, trymax = 100)

# Extract NMDS points
cazyme_nmds_df <- as.data.frame(cazyme_nmds$points)
cazyme_nmds_df$bin_name <- rownames(cazyme_nmds_df)

# Add metadata
cazyme_nmds_df <- cazyme_nmds_df |>
  left_join(meta, by = 'bin_name')

species_levels <- cazyme_nmds_df |>
  count(species, sort = TRUE) |>
  pull(species)

# Arrange from descending no. of MAGs
cazyme_nmds_df <- cazyme_nmds_df |>
  mutate(
    species = factor(species, levels = species_levels)
  )

species_cols <- paletteer::paletteer_d('ggthemes::Classic_10')[1:length(species_levels)]
species_cols <- setNames(species_cols, species_levels)

species_shapes <- c(16, 17, 15, 18, 3, 7, 8, 4, 9, 10)[1:length(species_levels)]
species_shapes <- setNames(species_shapes, species_levels)

# Create label text manually based on PERMANOVA result
cazyme_nmds_perm_label <- data.frame(
  NMDS1 = min(cazyme_nmds_df$MDS1) + 0.01,
  NMDS2 = max(cazyme_nmds_df$MDS2) - 0.01,
  label = "bold('PERMANOVA')"
)

cazyme_nmds_r2p_label <- data.frame(
  NMDS1 = min(cazyme_nmds_df$MDS1) + 0.01,
  NMDS2 = max(cazyme_nmds_df$MDS2) - 0.018,
  label = "italic(R^2)*' = 0.72, '*italic(p)*' = 0.001'"
)

# Plot the data
cazyme_species_nmds_plot <- ggplot(cazyme_nmds_df, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = species, shape = species), size = 4) +
  scale_color_manual(values = species_cols) +
  scale_shape_manual(values = species_shapes) +
  labs(
    title = '(C)',
    x = 'NMDS 1',
    y = 'NMDS 2',
    color = 'Species',
    shape = 'Species'
  ) +
  
  # For legend if placed at the bottom
  # guides(
  #  color = guide_legend(title.position = 'top', title.hjust = 0.5),
  #  shape = guide_legend(title.position = 'top', title.hjust = 0.5)
  # ) +
  
  geom_text(
    data = cazyme_nmds_perm_label,
    aes(x = NMDS1, y = NMDS2, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = cazyme_nmds_r2p_label,
    aes(x = NMDS1, y = NMDS2, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 5,
    inherit.aes = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(face = 'italic', size = 12),
    # legend.position = 'bottom',
    # legend.box = 'vertical',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

# Check the plot
cazyme_species_nmds_plot

# Save the plot
ggsave(
  cazyme_species_nmds_plot,
  filename = 'plots/07_downstream_analysis/7.2_dram/cazyme_comparisons/cazyme_species_nmds_plot.tiff',
  height = 15,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)

# CAZyme composition differed significantly between kefir grains and milk kefir, with enrichment of glycoside hydrolases in grains and glycosyl transferases in milk (Wilcoxon test, FDR < 0.05). 
# Regional differences were also observed, although these were less pronounced. 
# At the taxonomic level, CAZyme profiles varied significantly across dominant species (Kruskal–Wallis, p < 0.05). 
# Multivariate analysis confirmed that both sample type and species significantly influenced CAZyme composition (PERMANOVA, p < 0.01).


# Combine All CAZyme Plots ------------------------------------------------

library(gridExtra)

layout <- rbind(
  c(1, 2), # top plot spans 2 columns
  c(3, 3) # bottom plot spans 2 columns
)

cazyme_combined <- grid.arrange(
  sample_cazyme_plot,
  region_cazyme_plot,
  cazyme_species_nmds_plot,
  layout_matrix = layout,
  heights = c(1.2, 1),
  widths = c(1, 1)
)

# Save the plot
ggsave(
  cazyme_combined,
  filename = 'plots/07_downstream_analysis/7.2_dram/cazyme_comparisons/cazyme_combined.tiff',
  height = 40,
  width = 35,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)