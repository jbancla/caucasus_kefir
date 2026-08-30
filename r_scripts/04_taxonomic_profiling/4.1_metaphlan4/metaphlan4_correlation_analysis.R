# MetaPhlAn4 Correlation Analysis
# ~/Documents/RStudio/mk_caucasus/scripts/04_taxonomic_profiling/4.1_metaphlan4
# 17/03/2026

# Load the necessary packages --------------------------------------------
library(tidyverse)


# Load and Check Data -----------------------------------------------------

# Load relative abundance data from MetaPhlAn4
mp4_species_data <- read_tsv('data/04_taxonomic_profiling/metaphlan4_latest/mp4_caucasus_merged_species.tsv', show_col_types = FALSE)

# Convert to long format
mp4_species_long <- mp4_species_data |>
  pivot_longer(-clade_name, names_to = 'sample', values_to = 'abundance')

# Filter milk samples
meta_milk <- mp4_species_long |> 
  filter(str_detect(sample, '_m$')) # only milk samples

# Create pH metadata
ph_data <- tibble(
  sample = c('k1_m', 'k2_m', 'k3_m', 'k4_m', 'k5_m', 'k6_m', 'k7_m', 'k8_m',
             'k9_m', 'k10_m', 'k11_m', 'k12_m', 'k13_m', 'k14_m', 'k15_m', 'k16_m'),
  pH = c(6.3, 4.57, 6.33, 4.46, 5.89, 6.04, 5.56, 5.56,
         4.56, 4.42, 4.55, 4.56, 4.64, 4.58, 4.51, 4.68)
)

# Merge data
meta_ph <- mp4_species_long |>
  inner_join(ph_data, by = 'sample')


# Identify Top Species for Correlation ------------------------------------

# Calculate mean abundance per species
top_species_df <- meta_milk |>
  group_by(clade_name) |>
  summarise(mean_abundance = mean(abundance)) |>
  arrange(desc(mean_abundance))

top_species_df |> head(10)

top_species <- c(
  'Lactococcus_lactis',
  'Enterococcus_durans',
  'Enterococcus_faecalis',
  'Lacticaseibacillus_paracasei',
  'Leuconostoc_mesenteroides'
)

# Filter selected species
mp4_filtered <- meta_ph |>
  filter(clade_name %in% top_species)


# Run Spearman Correlation ------------------------------------------------

mp4_cor_results <- mp4_filtered |>
  group_by(clade_name) |>
  summarise(
    rho = cor(abundance, pH, method = 'spearman'),
    p_value = cor.test(abundance, pH, method = 'spearman')$p.value
  )

mp4_cor_results

# Adjust p-values
mp4_cor_results <- mp4_cor_results |>
  mutate(p_adj = p.adjust(p_value, method = 'BH'))

mp4_cor_results

# Filter only significant species
sig_species <- mp4_cor_results |>
  filter(p_adj < 0.05) |>
  pull(clade_name)

mp4_sig <- mp4_filtered |>
  filter(clade_name %in% sig_species)


# Plot the Significant Species --------------------------------------------

# Manually place correlation labels for significant species
mp4_labels_manual <- tibble(
  clade_name = c(
    'Enterococcus_durans',
    'Lactococcus_lactis'
  ),
  x = c(98, 98),   # adjust these manually
  y = c(4.25, 6.3), # adjust these manually
  label = c(
    'rho = 0.63\nadj. p = 0.022',
    'rho = -0.70\nadj. p = 0.012'
  )
) |>
  mutate(
    clade_name = factor(clade_name, levels = unique(mp4_sig$clade_name))
  )

mp4_cor_plot <- ggplot(mp4_sig, aes(x = abundance, y = pH)) +
  geom_smooth(
    method = 'lm', 
    color = '#FAB12F',
    fill = '#FFF6C0',
    se = TRUE) +
  geom_point(color = '#624E88', size = 3) +
  facet_wrap(~ clade_name, scales = 'free_x',
    labeller = as_labeller(
      function(x) paste0('bolditalic(', gsub('_', '~', x), ')'),
      label_parsed
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100)
  ) +
  geom_text(
    data = mp4_labels_manual,
    mapping = aes(x = x, y = y, label = label),
    hjust = 1, # alignment right (left to right)
    vjust = 1, # alignment right (top to bottom)
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    title = '(A)',
    subtitle = 'Taxonomic',
    x = 'Relative Abundance (%)',
    y = 'pH'
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    plot.subtitle = element_text(face = 'bold', size = 12, hjust = 0.5),
    strip.text.x = element_text(size = 12),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12)
  )

# Check the plot
mp4_cor_plot

# Save the plot
ggsave(
  mp4_cor_plot,
  filename = 'plots/04_taxonomic_profiling/metaphlan4_latest/mp4_correlation_pH.tiff',
  height = 15,
  width = 30,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)