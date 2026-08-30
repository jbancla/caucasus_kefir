# Metaphlan4 Taxonomic Profile
# ~/Documents/RStudio/mk_caucasus/scripts/04_taxonomic_profiling/4.1_metaphlan4
# 06/04/2026

# Load the necessary packages.
library(tidyverse)


# Load and Check Data -----------------------------------------------------

# Load data
mp4_species_data <- read_tsv('data/04_taxonomic_profiling/metaphlan4/mp4_caucasus_merged_species.tsv', show_col_types = FALSE)

# Define sample order
mp4_sample_order <- c('k1_g', 'k3_g', 'k5_g', 'k6_g', 'k7_g', 'k9_g',
                      'k1_m', 'k2_m', 'k3_m', 'k5_m', 'k6_m', 'k7_m', 
                      'k8_m',  'k9_m', 'k4_m', 'k10_m', 'k11_m',  'k12_m', 
                      'k13_m', 'k14_m', 'k15_m', 'k16_m')


# Species Level -----------------------------------------------------------

# Prepare data
mp4_species_long <- mp4_species_data |>
  pivot_longer(
    cols = all_of(mp4_sample_order),
    names_to = 'sample',
    values_to = 'abundance'
  )

# Split into >=0.1 (keep) and <0.1 (pool per sample)
mp4_species_keep <- mp4_species_long |> filter(abundance >= 0.1)

mp4_species_others <- mp4_species_long |>
  filter(abundance < 0.1) |>
  group_by(sample) |>
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = 'drop') |>
  mutate(clade_name = 'Others (< 0.1%)')

# Combine: keep (>=0.1) + Others (<0.1%)
mp4_species_combined <- bind_rows(
  mp4_species_keep,
  mp4_species_others
) |>
  mutate(sample = factor(sample, levels = mp4_sample_order))

# Adding collection sites layer
kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkaria (KB)' = '#4FACAC', 'Karachay-Cherkessia (KC)' = '#E8403D')

mp4_species_combined <- mp4_species_combined |>
  mutate(collection_site = if_else(sample %in% kab_samples,
                                   'Kabardino-Balkaria (KB)', 'Karachay-Cherkessia (KC)'))

# Build colored x-axis labels (HTML spans) per sample
sample_levels <- levels(mp4_species_combined$sample)
sample_site   <- mp4_species_combined |>
  distinct(sample, collection_site) |>
  mutate(col = site_cols[collection_site]) |>
  arrange(match(sample, sample_levels))

axis_labels <- setNames(
  paste0('<span style="color:', sample_site$col, '">', sample_site$sample, '</span>'),
  sample_site$sample
)

# Order species by total abundance (desc)
rest_order <- mp4_species_combined |>
  filter(!clade_name %in% c('Others (< 0.1%)')) |>
  group_by(clade_name) |>
  summarise(total = sum(abundance, na.rm = TRUE), .groups = 'drop') |>
  arrange(desc(total)) |>
  pull(clade_name)

# Stacking order (bars): Others (bottom), species (top)
stack_order <- c('Others (< 0.1%)', rest_order)

# Legend order (entries): species, Others
legend_order <- c(gsub('_', ' ', rest_order), 'Others (< 0.1%)')

mp4_species_combined <- mp4_species_combined |>
  mutate(
    clade_name  = factor(clade_name, levels = stack_order), # for stacking
    clade_label = case_when(                                # for legend text & fill mapping
      clade_name == 'Others (< 0.1%)' ~ as.character(clade_name),
      TRUE ~ gsub('_', ' ', as.character(clade_name))
    ),
    clade_label = factor(clade_label, levels = legend_order),
    sample_group = ifelse(grepl('_g$', as.character(sample)), 'Kefir Grains', 'Milk Kefir')
  )

# Colors: use ggthemes::Classic_20 for species + Others
n_species <- length(rest_order)
# We need n_species + 1 colors (species + 'Others (<0.1%)'); Classic_20 has 20 colors
classic_cols <- as.character(
  paletteer::paletteer_d('ggthemes::Classic_20', n = max(1, n_species + 1))
)
species_cols <- classic_cols[seq_len(n_species)]
others_col  <- classic_cols[n_species + 1]

color_map <- setNames(
  c(species_cols, others_col),
  c(gsub('_', ' ', rest_order), 'Others (< 0.1%)')
)

# compute the largest stacked total across samples
xmax <- mp4_species_combined |>
  group_by(sample) |>
  summarise(total = sum(abundance, na.rm = TRUE), .groups = 'drop') |>
  summarise(max_total = max(total, na.rm = TRUE)) |>
  pull(max_total)

# nice breaks & tiny right headroom
upper  <- xmax * 1.02
xbreaks <- pretty(c(0, upper), n = 8)

# facet-aware dummy data to avoid gaps
dummy_site <- mp4_species_combined |>
  dplyr::distinct(sample, collection_site, sample_group)

# Plot the data
mp4_species_plot <- ggplot(mp4_species_combined, aes(x = sample, y = abundance, fill = clade_label)) +
  geom_col(width = 0.8, position = 'stack', color = 'black', linewidth = 0.2) +
  # facet-aware dummy layer for site legend (no gaps)
  geom_point(
    data = dummy_site,
    aes(x = sample, y = 0, color = collection_site),
    inherit.aes = FALSE, alpha = 0, size = 0, show.legend = TRUE
  ) +
  # facet_wrap(~ sample_group, scales = 'free_x') +
  facet_grid(. ~ sample_group, scales = 'free_x', space = 'free_x') + # << make bar widths consistent
  scale_fill_manual(
    values = color_map, name = 'Species',
    labels = function(x) {
      lab <- sapply(x, function(name) {
        if (name == 'Others (< 0.1%)') paste0("'", name, "'")
        else paste0('italic(\"', name, '\")')
      })
      parse(text = lab)
    }
  ) +
  scale_color_manual(values = site_cols, name = 'Caucasus Region') +
  # scale_y_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(limits = c(0, upper),
                     breaks  = xbreaks,
                     expand  = expansion(mult = c(0, 0))) +
  scale_x_discrete(labels = axis_labels) +
  labs(
    # title = 'Relative Abundance of Caucasus Kefirs at Species Level (MetaPhlAn4)',
    title = '(A)',
    x = 'Samples', y = 'Relative Abundance (%)'
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    axis.text.x   = ggtext::element_markdown(angle = 45, hjust = 1),
    strip.text.x  = element_text(face = 'bold', size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    legend.position = 'right', # keep legends on the right
    legend.box = 'vertical', # stack them vertically
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill  = guide_legend(order = 1, title.position = 'top'), # Species on top
    color = guide_legend(order = 2, title.position = 'top',
                         override.aes = list(alpha = 1, size = 4)) # Collection site below
  )

# Check the plot
mp4_species_plot

# Save the plot
ggsave(
  mp4_species_plot,
  filename = 'plots/04_taxonomic_profiling/metaphlan4/mp4_species_region.tiff',
  height = 20,
  width = 30,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)