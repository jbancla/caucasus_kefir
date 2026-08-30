# GTDB-Tk Taxonomic Profile
# ~/Documents/RStudio/mk_caucasus/scripts/07_downstream_analysis/7.1_gtdbtk
# 02/10/2025

# Load the necessary packages.
library(tidyverse)
library(paletteer)


# Load and Check Data -----------------------------------------------------

# prefixes that should map to Karachay-Cherkessia (kc)
kc_prefix <- c('k4_m','k10_m','k11_m','k12_m','k13_m','k14_m','k15_m','k16_m')
kc_regex  <- paste0('^(', paste(kc_prefix, collapse = '|'), ')')

# list of samples to remove
remove_prefix <- c('k12_g','s1.2_m','s7.2_m','s9.2_m','sfrankenstein')
remove_regex  <- paste0('^(', paste(remove_prefix, collapse = '|'), ')')

# Load data
gtdb_hqmags_data <- read_tsv(
  'data/07_downstream_analysis/7.1_gtdbtk/hq_mags/classify/gtdbtk.bac120.summary.tsv',
  show_col_types = FALSE
) |>
  select(user_genome, classification, closest_genome_ani) |>
  filter(!str_detect(user_genome, remove_regex)) |>
  mutate(
    sample_type = if_else(str_detect(user_genome, '_g'), 'Kefir Grains', 'Milk Kefir'),
    collection_site = if_else(str_detect(user_genome, kc_regex),
                              'Karachay-Cherkessia', 'Kabardino-Balkaria')
  ) |>
  # split classification into ranks (d__;p__;c__;o__;f__;g__;s__)
  separate_wider_delim(
    classification, delim = ';',
    names = c('domain','phylum','class','order','family','genus','species'),
    too_few = 'align_start'
  ) |>
  # drop x__ prefixes and make them pretty
  mutate(
    across(domain:species, ~ str_remove(., '^[a-z]__')),
    genus   = str_replace_all(genus,   '_', ' '),
    species = str_replace_all(species, '_', ' ')
  )

checkm2_hqmags_data <- read_tsv('data/03_bin_classification/checkm2/hq_mags/quality_report.tsv', show_col_types = FALSE) |> 
  filter(!str_detect(Name, remove_regex)) |> 
  rename(completeness = 'Completeness',
         contamination = 'Contamination')

checkm2_quality <- checkm2_hqmags_data |> 
  select(completeness, contamination)

# Quick sanity check to see if the list and arrangement of user_genome from GTDB-Tk and Name in Checkm2 are the same
stopifnot(nrow(gtdb_hqmags_data) == nrow(checkm2_hqmags_data))
stopifnot(identical(gtdb_hqmags_data$user_genome, checkm2_hqmags_data$Name))
message(all(identical(nrow(gtdb_hqmags_data), nrow(checkm2_hqmags_data)), identical(gtdb_hqmags_data$user_genome, checkm2_hqmags_data$Name)))

# Append completeness and contamination column to gtdb-tk result
gtdb_hqmags_final_data <- bind_cols(gtdb_hqmags_data, checkm2_quality)

# Removing MAGs with <90 completeness
gtdb_hqmags_final_data <- gtdb_hqmags_final_data |> 
  filter(completeness >= 90)

write_tsv(gtdb_hqmags_final_data, 'data/07_downstream_analysis/7.1_gtdbtk/hq_mags/gtdb_hqmags_with_qc_data.tsv')

# Count the number of species in sample_type and collection_site
species_counts <- gtdb_hqmags_final_data |> 
  group_by(species, sample_type, collection_site) |> 
  summarise(n_MAGs = n(), .groups = 'drop') |> 
  arrange(species, sample_type, collection_site)

# Count the total number of species
species_totals <- species_counts |>
  group_by(species) |>
  summarise(total_MAGs = sum(n_MAGs))

# Converting species_counts dataset to matirix for heatmap
species_matrix <- species_counts |>
  pivot_wider(names_from = sample_type,
              values_from = n_MAGs,
              values_fill = 0)


# Stack Barplot -----------------------------------------------------------

# Plot the data
gtdb_species_count_a <- species_counts |>
  mutate(species = fct_reorder(species, n_MAGs, .fun = sum)) |> 
  ggplot(aes(x = species, y = n_MAGs, fill = sample_type)) +
  geom_col(position = 'stack') +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  scale_fill_manual(
    values = c('Milk Kefir' = '#FAB12F', 'Kefir Grains' = '#624E88'),
    name   = 'Sample Type'
  ) +
  coord_flip() +
  facet_wrap(~ collection_site, strip.position = 'right') +
  labs(
    # title = 'Species Across All Samples in each Collection Sites',
    title = '(B)',
    x = NULL, 
    y = 'No. of High Quality MAGs', 
    fill = 'Sample Type') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.y = element_text(face = 'italic', size = 10, color = 'black'),
    axis.text.x = element_text(size = 12, color = 'black'),
    strip.text.y = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Check the plot
gtdb_species_count_a

# Save the plot
ggsave(
  gtdb_species_count_a,
  filename = 'plots/07_downstream_analysis/7.1_gtdbtk/gtdb_species_counts_stack.tiff',
  height = 15,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# MAGs Completeness -------------------------------------------------------

library(ggnewscale)

# 1) desired prefix order
prefix_order <- c('k1_g','k3_g','k5_g','k6_g','k7_g','k9_g', paste0('k', 1:16, '_m'))

# 2) add a prefix column and lock its order
gtdb_ordered <- gtdb_hqmags_final_data |>
  mutate(
    prefix = str_extract(user_genome, '^k\\d+_[gm]'),
    prefix = factor(prefix, levels = prefix_order)
  ) |>
  # 3) order rows by prefix, then (optionally) by completeness descending
  arrange(prefix, desc(completeness), user_genome) |>
  # 4) freeze user_genome factor in that order (reverse for coord_flip)
  mutate(user_genome = factor(user_genome, levels = rev(unique(user_genome))))

# Plot the data
gtdb_lollipop_plot <- gtdb_ordered |>
  mutate(species_clean = stringr::str_replace_all(species, '_', ' ')) |>
  ggplot(aes(x = user_genome, y = completeness)) +
  
  # 1) LINES by sample_type (fixed kefir colors)
  geom_segment(
    aes(xend = user_genome, y = 0, yend = completeness, colour = sample_type),
    linewidth = 0.5
  ) +
  scale_colour_manual(
    values = c('Milk Kefir' = '#FAB12F', 'Kefir Grains' = '#624E88'),
    name   = 'Sample Type'
  ) +
  
  # reset colour scale for points
  ggnewscale::new_scale_colour() +
  
  # 2) POINTS by species using Paletteer (Classic_10)
  geom_point(aes(colour = species_clean), size = 3) +
  paletteer::scale_color_paletteer_d(
    palette = 'ggthemes::Classic_10',
    name    = 'Species',
    labels  = function(x) parse(text = paste0('italic("', x, '")'))
  ) +
  
  # extras
  geom_hline(yintercept = 90, linetype = 'dashed', colour = 'red') +
  coord_flip() +
  # facet_wrap(~ collection_site, ncol = 1, scales = 'free_y') +
  facet_grid(collection_site ~ ., scales = 'free_y', space = 'free_y') + # << stacked + consistent spacing
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    # title = 'MAGs Completeness and Taxonomic Profile (CheckM2 and GTDB-Tk)'),
    title = '(A)',
    x = 'High Quality MAGs', 
    y = 'Completeness (%)') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.x = element_text(size = 12, color = 'black'),
    axis.text.y = element_text(size = 8, color = 'black'),
    strip.text.y = element_text(face = 'bold', size = 10),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Check the plot
gtdb_lollipop_plot

# Save the plot
ggsave(
  gtdb_lollipop_plot,
  filename = 'plots/07_downstream_analysis/7.1_gtdbtk/gtdb_species_lollipop.tiff',
  height = 20,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)