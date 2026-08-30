# DRAM MAGs SCFA and Alcohol Conversions
# ~/Documents/RStudio/mk_caucasus/scripts/07_downstream_analysis/7.2_dram
# 21/01/2026

# Load the necessary packages.
library(tidyverse)


# Load and Check Data -----------------------------------------------------

dram_hqmags_data <- read_tsv(
  'data/07_downstream_analysis/7.2_dram/hq_mags/dram_hqmags_with_taxonomy_and_qc_data.tsv',
  show_col_types = FALSE
)

# Keep an ID column for each MAG on the x-axis (change 'genome' if needed)
dram_hqmags_scfa <- dram_hqmags_data |>
  select(genome, species, sample_type, collection_site, 88:100)

metadata_columns <- c('genome', 'species', 'sample_type', 'collection_site')
metafunc_columns <- setdiff(names(dram_hqmags_scfa), metadata_columns)

# Order species by how many MAGs they have (most -> least)
species_order <- dram_hqmags_scfa |>
  distinct(species, genome) |>
  count(species, name = 'n_mags') |>
  arrange(desc(n_mags), species) |>
  pull(species)

# Long format for plotting
dram_hqmags_scfa_long <- dram_hqmags_scfa |>
  mutate(species = factor(species, levels = species_order)) |>
  pivot_longer(
    cols = all_of(metafunc_columns),
    names_to = 'metabolic_function',
    values_to = 'present'
  ) |>
  mutate(
    # Remove DRAM category prefix
    metabolic_function = stringr::str_remove(
      metabolic_function,
      '^SCFA and alcohol conversions:\\s*'
    ),
    # Explicit factor levels control legend order
    presence = factor(
      if_else(present, 'Present', 'Absent'),
      levels = c('Present', 'Absent')
    ),
    metabolic_function = factor(
      metabolic_function,
      levels = rev(unique(metabolic_function))
    )
  )


# Plot the Data -----------------------------------------------------------

scfa_plot <- ggplot(
  dram_hqmags_scfa_long,
  aes(x = genome, y = metabolic_function, fill = presence)
) +
  geom_tile(
    color = 'white',
    linewidth = 0.2,
    width = 0.90,
    height = 0.90
  ) +
  facet_grid(
    . ~ species,
    scales = 'free_x',
    space = 'free_x',
    labeller = labeller(
      species = function(x) sub(' ', '\n', x)
    )
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), position = 'right') + # moves the y-axis text to the right
  scale_fill_manual(
    name = 'Presence/Absence',
    values = c('Present' = '#2BAE66', 'Absent' = '#EDEDED')
  ) +
  labs(
    title = '(D)',
    x = 'No. of HQ MAGs per Species',
    y = 'SCFA and Alcohol Conversion'
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(margin = margin(l = 10)),
    axis.text.y = element_text(size = 11, hjust = 0),
    axis.ticks.y = element_line(),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.07, 'cm'),
    panel.border = element_rect(colour = 'grey40', fill = NA, linewidth = 0.5),
    strip.text.x = element_text(face = 'italic', size = 9, angle = 90),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    legend.position = 'bottom'
  )

# Check the plot
scfa_plot

# Save the plot
ggsave(
  scfa_plot,
  filename = 'plots/07_downstream_analysis/7.2_dram/dram_scfa_plot.tiff',
  height = 10,
  width = 30,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine the Plots of CAZymes Comparisons and SCFA -----------------------

library(gridExtra)

layout_1 <- rbind(
  c(1, 2),
  c(3, 3),
  c(4, 4)
)

cazyme_scfa_combined <- grid.arrange(
  sample_cazyme_plot,
  region_cazyme_plot,
  cazyme_species_nmds_plot,
  scfa_plot,
  layout_matrix = layout_1,
  heights = c(1.2, 1, 0.8),
  widths = c(1, 1)
)

# Save the plot
ggsave(
  cazyme_scfa_combined,
  filename = 'plots/07_downstream_analysis/7.2_dram/cazyme_comparisons/cazyme_scfa_combined.tiff',
  height = 50,
  width = 35,
  units = 'cm',
  bg = 'white',
  dpi = 600,
  compression = 'lzw'
)