# HUMAnN3 Functional Profile - Stratified
# ~/Documents/RStudio/mk_caucasus/scripts/05_functional_profiling/5.1_humann3
# 10/04/2026

# Include All Species Mentioned in taxonomic Composition ------------------

# Load the necessary packages.
library(ggtext)
library(paletteer)
library(tidyverse)
library(ggthemes)

## Load and Check Data ----------------------------------------------------

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

# Load data
hm3_strat_data <- read_tsv(
  'data/05_functional_profiling/humann3/caucasus_stratified_pathabundance_cpm/hm3_caucasus_joined_pathabundance_cpm_stratified.tsv',
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_Abundance-CPM$')) |>
  filter(!str_detect(pathway, 'UNINTEGRATED')) |> # alternative: use base R filter function: filter(!grepl('UNINTEGRATED', pathway))
  
  # split pathway and taxonomy
  separate(
    pathway,
    into = c('pathway', 'taxonomy'),
    sep = '\\|',
    remove = TRUE) |>
  filter(!grepl('unclassified', taxonomy)) |> # removing unclassified taxa in taxonomy column
  
  # extract genus and species
  mutate(
    genus = str_extract(taxonomy, '(?<=g__)[^\\.]+'),
    species = str_extract(taxonomy, '(?<=s__).*'),
    species = str_replace_all(species, '_', ' ')
  ) |>
  
  # drop the raw taxonomy column if not needed
  select(-taxonomy) |>
  
  # collapse everything else to "Others"
  mutate(
    species = if_else(species %in% keep_species, species, 'Others')
  )

# 1) Pivot longer so we can calculate prevalence by sample
hm3_long <- hm3_strat_data |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = 'sample',
    values_to = 'abundance'
  )

# 2) Community-level per sample × pathway (sum across species)
hm3_comm <- hm3_long |>
  group_by(sample, pathway) |>
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = 'drop')

# 3) Pathway prevalence across samples
n_samples <- n_distinct(hm3_comm$sample)

pathway_prevalence <- hm3_comm |>
  mutate(is_present = abundance > 0) |>
  group_by(pathway) |>
  summarise(
    n_present = sum(is_present, na.rm = TRUE),
    prevalence = n_present / n_samples,
    .groups = 'drop'
  ) |>
  mutate(
    presence_group = case_when(
      prevalence == 1 ~ 'core',
      prevalence >= 0.70 & prevalence < 1 ~ 'high',
      prevalence >= 0.50 & prevalence < 0.70 ~ 'moderate',
      prevalence > 0 & prevalence < 0.50 ~ 'low',
      prevalence == 0 ~ 'sample_specific'  # usually won’t occur if pathway exists at all
    )
  )

# 4) Add pathway group back to stratified long data (or to hm3_strat_data)
hm3_strat_data_grouped <- hm3_strat_data |>
  left_join(pathway_prevalence, by = 'pathway')

# 5) Split into separate objects
hm3_core_strat_data <- hm3_strat_data_grouped |> filter(presence_group == 'core')
hm3_high_strat_data <- hm3_strat_data_grouped |> filter(presence_group == 'high')
hm3_moderate_strat_data <- hm3_strat_data_grouped |> filter(presence_group == 'moderate')
hm3_low_strat_data <- hm3_strat_data_grouped |> filter(presence_group == 'low')
hm3_sample_specific_strat_data <- hm3_strat_data_grouped |> filter(presence_group == 'sample_specific')

# Helper function: remove pathways with zero abundance across all samples
remove_all_zero_pathways <- function(df) {
  df |>
    rowwise() |>
    filter(any(c_across(where(is.numeric)) > 0)) |>
    ungroup()
}

hm3_core_strat_data <- remove_all_zero_pathways(hm3_core_strat_data)
hm3_high_strat_data <- remove_all_zero_pathways(hm3_high_strat_data)
hm3_moderate_strat_data <- remove_all_zero_pathways(hm3_moderate_strat_data)
hm3_low_strat_data <- remove_all_zero_pathways(hm3_low_strat_data)
hm3_sample_specific_strat_data <- remove_all_zero_pathways(hm3_sample_specific_strat_data)

# Count the exact number of unique pathways per category
pathway_counts <- tibble(
  category = c('core', 'high', 'moderate', 'low', 'sample_specific'),
  n_pathways = c(
    n_distinct(hm3_core_strat_data$pathway),
    n_distinct(hm3_high_strat_data$pathway),
    n_distinct(hm3_moderate_strat_data$pathway),
    n_distinct(hm3_low_strat_data$pathway),
    n_distinct(hm3_sample_specific_strat_data$pathway)
  )
)

# Print
pathway_counts


# With Addition of Collection Site ----------------------------------------

# Core Metabolic Pathways -------------------------------------------------

# adding pathway_ID and order alphabetically
hm3_core_strat_data <- hm3_core_strat_data |>
  mutate(
    pathway_ID = str_trim(str_extract(pathway, '^[^:]+'))
  ) |>
  arrange(pathway_ID) |>
  mutate(
    pathway_ID = factor(pathway_ID, levels = unique(pathway_ID)),
    pathway    = factor(pathway,    levels = unique(pathway))
  )

# Define sample order
sample_order <- c(
  'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g',
  'k1_m','k2_m','k3_m','k5_m','k6_m','k7_m',
  'k8_m','k9_m','k4_m','k10_m','k11_m','k12_m',
  'k13_m','k14_m','k15_m','k16_m'
)

# Species colors (Classic_20 for keep_species + Others as grey)
classic20 <- paletteer::paletteer_d('ggthemes::Classic_20')
classic20 <- classic20[seq_along(keep_species)]
names(classic20) <- keep_species

species_cols <- c(classic20, 'Others' = '#EDEDED')

# Legend labels with italics (Others stays plain)
species_labels <- c(
  setNames(paste0('<i>', keep_species, '</i>'), keep_species),
  'Others' = 'Others'
)

# Adding collection sites layer ------------------------------------------

kab_samples <- c(
  'k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
  'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g'
)

site_cols <- c(
  'Kabardino-Balkaria' = '#4FACAC',
  'Karachay-Cherkessia' = '#E8403D'
)

# Long CPM table (no relative normalization)
hm3_core_long_cpm <- hm3_core_strat_data |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = 'sample',
    values_to = 'cpm'
  ) |>
  mutate(
    cpm = replace_na(cpm, 0),  # prevents NA bars appearing
    sample = factor(sample, levels = sample_order),
    species = factor(species, levels = c(keep_species, 'Others')),
    pathway_ID = factor(pathway_ID, levels = sort(unique(pathway_ID))),
    collection_site = if_else(as.character(sample) %in% kab_samples,
                              'Kabardino-Balkaria', 'Karachay-Cherkessia')
  ) |>
  filter(!is.na(sample))      # drops any samples not in your sample_order

# Build colored x-axis labels (HTML spans) per sample
sample_levels <- levels(hm3_core_long_cpm$sample)

sample_site <- hm3_core_long_cpm |>
  distinct(sample, collection_site) |>
  mutate(col = site_cols[collection_site]) |>
  arrange(match(sample, sample_levels))

sample_labels <- setNames(
  paste0('<span style="color:', sample_site$col, ';">', sample_site$sample, '</span>'),
  sample_site$sample
)

# Plot the data -----------------------------------------------------------

hm3_core_barplot <- ggplot(
  hm3_core_long_cpm,
  aes(x = sample, y = cpm, fill = species)
) +
  geom_col(width = 0.9, color = 'black', linewidth = 0.2) +
  geom_point(aes(color = collection_site), alpha = 0, size = 5) + # <- dummy layer
  facet_wrap(~ pathway_ID, ncol = 4, scales = 'free_y') +
  
  scale_fill_manual(
    values = species_cols,
    breaks = c(keep_species, 'Others'),
    labels = species_labels,
    drop = FALSE
  ) +
  scale_color_manual(
    values = site_cols,
    name = 'Caucasus Region'
  ) +
  
  scale_x_discrete(labels = sample_labels) +
  guides(
    fill  = guide_legend(order = 1),
    color = guide_legend(order = 2, override.aes = list(alpha = 1))
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  labs(
    x = 'Samples',
    y = 'Copies per million (CPM)',
    fill = 'Contributing Species'
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 12, face = 'bold', margin = margin(t = 10)),
    axis.text = element_text(size = 12),
    axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(size = 12, face = 'bold', margin = margin(r = 10)),
    #strip.background = element_rect(fill = 'white'),
    strip.text = element_text(face = 'bold', size = 12),
    legend.title = element_text(size = 12, face = 'bold'),
    legend.text = ggtext::element_markdown(size = 12)
  )

# Check the plot
hm3_core_barplot

# Save the plot
ggsave(
  hm3_core_barplot,
  filename = 'plots/05_functional_profiling/humann3/hm3_core_barplot_all.tiff',
  height = 50,
  width = 60,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Determine how many species are involved in each pathway -----------------

hm3_core_long2 <- hm3_core_strat_data |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = 'sample',
    values_to = 'cpm'
  ) |>
  filter(cpm > 0)

species_pathway_prevalence <- hm3_core_long2 |>
  distinct(species, pathway) |>
  count(species, name = 'n_core_pathways')