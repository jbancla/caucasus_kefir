# HUMAnN3 pathway abundance and stratified contributors
# ~/Documents/RStudio/mk_caucasus/scripts/05_functional_profiling/5.1_humann3
# 08/04/26

# Load the necessary packages --------------------------------------------
library(tidyverse)
library(paletteer)
library(ggtext)


# Define metadata ---------------------------------------------------------

# Define sample order
sample_order <- c('k1_g','k3_g','k5_g','k6_g','k7_g','k9_g',
                  'k1_m','k2_m','k3_m','k5_m','k6_m','k7_m',
                  'k8_m','k9_m','k4_m','k10_m','k11_m','k12_m',
                  'k13_m','k14_m','k15_m','k16_m')

# Define Kabardino-Balkaria samples
kab_samples <- c('k1_g','k3_g','k5_g','k6_g','k7_g','k9_g',
                 'k1_m','k2_m','k3_m','k5_m','k6_m','k7_m',
                 'k8_m','k9_m')


# Unstratified pathway abundance ------------------------------------------

# Load unstratified pathway abundance data
hm3_unstrat_data <- read_tsv(
  'data/05_functional_profiling/humann3/caucasus_stratified_pathabundance_cpm/hm3_caucasus_joined_pathabundance_cpm_unstratified.tsv',
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_Abundance-CPM$')) |>
  filter(!pathway %in% c('UNMAPPED', 'UNINTEGRATED'))

# Convert to long format
hm3_unstrat_long <- hm3_unstrat_data |>
  pivot_longer(cols = -pathway, names_to = 'sample', values_to = 'abundance') |>
  mutate(
    sample = factor(sample, levels = sample_order),
    sample_type = if_else(str_detect(sample, '_g$'), 'Kefir Grains', 'Milk Kefir'),
    collection_site = if_else(sample %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia')
  )

# Compare kefir grains vs milk kefir using Wilcoxon rank-sum test
da_results <- hm3_unstrat_long |>
  group_by(pathway) |>
  summarise(
    p_value = wilcox.test(abundance ~ sample_type)$p.value,
    mean_grain = mean(abundance[sample_type == 'Kefir Grains'], na.rm = TRUE),
    mean_milk  = mean(abundance[sample_type == 'Milk Kefir'], na.rm = TRUE),
    .groups = 'drop'
  ) |>
  mutate(
    p_adj = p.adjust(p_value, method = 'BH'),
    enriched_in = if_else(mean_milk > mean_grain, 'Milk Kefir', 'Kefir Grains'),
    log2FC = log2((mean_milk + 1) / (mean_grain + 1))
  ) |>
  arrange(p_adj)

# Keep significant pathways
sig_pathways <- da_results |>
  filter(p_adj < 0.05)

# Save significant pathways
write_tsv(
  sig_pathways,
  'data/results/05_functional_profiling/humann3/diff_abundance_significant_pathways.tsv'
)

# Identify top 10 milk- and grain-enriched pathways
top10_milk <- sig_pathways |>
  filter(enriched_in == 'Milk Kefir') |>
  arrange(p_adj) |>
  slice_head(n = 10) |>
  pull(pathway)

top10_grain <- sig_pathways |>
  filter(enriched_in == 'Kefir Grains') |>
  arrange(p_adj) |>
  slice_head(n = 10) |>
  pull(pathway)

top20_pathways <- c(top10_grain, top10_milk)


# Stratified pathway contributors -----------------------------------------

# Selecting target species
keep_species <- c(
  'Lactococcus lactis',
  'Enterococcus durans',
  'Lactobacillus kefiranofaciens',
  'Lactobacillus kefiri',
  'Lactobacillus casei group',
  'Enterococcus faecalis',
  'Leuconostoc mesenteroides',
  'Pseudomonas helleri',
  'Lactobacillus helveticus',
  'Lactobacillus parakefiri',
  'Hafnia paralvei',
  'Enterococcus gilvus',
  'Bifidobacterium mongoliense',
  'Hafnia alvei',
  'Saccharomyces cerevisiae',
  'Serratia liquefaciens',
  'Enterococcus malodoratus'
)

# Load stratified pathway abundance data
hm3_strat_data <- read_tsv(
  'data/05_functional_profiling/humann3/caucasus_stratified_pathabundance_cpm/hm3_caucasus_joined_pathabundance_cpm_stratified.tsv',
  show_col_types = FALSE
) |>
  rename_with(~ str_remove(., '_Abundance-CPM$')) |>
  filter(!str_detect(pathway, 'UNINTEGRATED|unclassified')) |>
  separate(
    pathway,
    into = c('pathway', 'taxonomy'),
    sep = '\\|',
    remove = TRUE
  ) |>
  mutate(
    genus = str_extract(taxonomy, '(?<=g__)[^\\.]+'),
    species = str_extract(taxonomy, '(?<=s__).*'),
    species = str_replace_all(species, '_', ' '),
    species = if_else(species %in% keep_species, species, 'Others')
  ) |>
  select(-taxonomy)

# Keep only top 10 milk- and grain-enriched pathways
hm3_strat_top <- hm3_strat_data |>
  filter(pathway %in% top20_pathways)

# Convert to long format and add metadata
hm3_strat_long <- hm3_strat_top |>
  pivot_longer(
    cols = -c(pathway, genus, species),
    names_to = 'sample',
    values_to = 'abundance'
  ) |>
  mutate(
    sample = factor(sample, levels = sample_order),
    sample_type = if_else(str_detect(sample, '_g$'), 'Kefir Grains', 'Milk Kefir'),
    collection_site = if_else(sample %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia')
  )

# Summarise mean contribution per pathway, sample type, and species
hm3_strat_summary <- hm3_strat_long |>
  group_by(pathway, sample_type, species) |>
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = 'drop'
  )

# Keep top 10 contributing taxa per pathway and sample type
hm3_strat_plot <- hm3_strat_summary |>
  group_by(pathway, sample_type) |>
  arrange(desc(mean_abundance), .by_group = TRUE) |>
  mutate(rank = row_number()) |>
  mutate(
    species_plot = if_else(rank <= 10, species, 'Others')
  ) |>
  group_by(pathway, sample_type, species_plot) |>
  summarise(
    mean_abundance = sum(mean_abundance),
    .groups = 'drop'
  )

# Order pathways
hm3_strat_plot <- hm3_strat_plot |>
  mutate(
    pathway = factor(pathway, levels = top20_pathways),
    sample_type = factor(sample_type, levels = c('Kefir Grains', 'Milk Kefir'))
  )

# Order species in legend by total abundance across plotted data
species_order <- hm3_strat_plot |>
  group_by(species_plot) |>
  summarise(total_abundance = sum(mean_abundance, na.rm = TRUE), .groups = 'drop') |>
  arrange(desc(total_abundance)) |>
  pull(species_plot)

# Move 'Others' to the end
species_order <- c(setdiff(species_order, 'Others'), 'Others')

hm3_strat_plot <- hm3_strat_plot |>
  mutate(
    species_plot = factor(species_plot, levels = species_order)
  )

# Shortened Pathway Name
hm3_strat_plot <- hm3_strat_plot |>
  mutate(
    pathway_short = factor(
      str_extract(pathway, '^[^:]+'),
      levels = str_extract(top20_pathways, '^[^:]+')
    )
  )

# Create palette from Classic_20
n_species <- length(species_order)
species_cols <- paletteer_d('ggthemes::Classic_20', n = n_species)
names(species_cols) <- species_order

# Create italic labels (except 'Others')
species_labels <- setNames(
  ifelse(
    species_order == 'Others',
    'Others',
    paste0('<i>', species_order, '</i>')
  ),
  species_order
)

# Plot top contributors
hm3_strat_contributor_plot <- ggplot(
  hm3_strat_plot,
  aes(x = sample_type, y = mean_abundance, fill = species_plot)
) +
  geom_col(width = 0.8, color = 'black', linewidth = 0.2) +
  facet_wrap(~ pathway_short, scales = 'free_y', ncol = 5) +
  scale_fill_manual(
    values = species_cols,
    labels = species_labels,
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = '(A)',
    x = NULL,
    y = 'Mean Pathway Abundance (CPM)',
    fill = 'Species'
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    strip.text = element_text(face = 'bold', size = 12),
    axis.title.y = element_text(face = 'bold', size = 12, margin = margin(r = 10)),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_markdown(size = 12),
    panel.grid = element_blank()
  )

# Check the plot
hm3_strat_contributor_plot

# Save the plot
ggsave(
  hm3_strat_contributor_plot,
  filename = 'plots/05_functional_profiling/humann3/hm3_strat_contributor_plot.tiff',
  height = 40,
  width = 50,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)

# Save top contributors table
top_contributors_table <- hm3_strat_summary |>
  group_by(pathway, sample_type) |>
  arrange(desc(mean_abundance), .by_group = TRUE) |>
  slice_head(n = 3) |>
  ungroup()

write_tsv(
  top_contributors_table,
  'data/results/05_functional_profiling/humann3/top_contributors_selected_pathways.tsv'
)