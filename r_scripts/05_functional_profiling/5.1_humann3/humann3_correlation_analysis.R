# HUMAnN3 Correlation Analysis
# ~/Documents/RStudio/mk_caucasus/scripts/05_functional_profiling/5.1_humann3
# 18/03/2026

# Load the necessary packages --------------------------------------------
library(tidyverse)


# Load and Check Data -----------------------------------------------------

# Define sample order
sample_order <- c('k1_g','k3_g','k5_g','k6_g','k7_g','k9_g',
                  'k1_m','k2_m','k3_m','k5_m','k6_m','k7_m',
                  'k8_m','k9_m','k4_m','k10_m','k11_m','k12_m',
                  'k13_m','k14_m','k15_m','k16_m')

# Define Kabardino-Balkaria samples
kab_samples <- c('k1_g','k3_g','k5_g','k6_g','k7_g','k9_g',
                 'k1_m','k2_m','k3_m','k5_m','k6_m','k7_m',
                 'k8_m','k9_m')


# Load unstratified pathway abundance data
hm3_unstrat_data <- read_tsv(
  'data/05_functional_profiling/humann3_caucasus/stratified_pathabundance_cpm/hm3_joined_pathabundance_cpm_unstratified.tsv',
  show_col_types = FALSE
) |>
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

# Create pH metadata for milk kefir
ph_data <- tibble(
  sample = c('k1_m', 'k2_m', 'k3_m', 'k4_m', 'k5_m', 'k6_m', 'k7_m', 'k8_m',
             'k9_m', 'k10_m', 'k11_m', 'k12_m', 'k13_m', 'k14_m', 'k15_m', 'k16_m'),
  pH = c(6.3, 4.57, 6.33, 4.46, 5.89, 6.04, 5.56, 5.56,
         4.56, 4.42, 4.55, 4.56, 4.64, 4.58, 4.51, 4.68)
)


# Correlation between top pathway abundance and pH ------------------------

# Keep only milk kefir samples and merge with pH
hm3_ph <- hm3_unstrat_long |>
  filter(sample_type == 'Milk Kefir') |>
  inner_join(ph_data, by = 'sample')

# Identify top pathways by mean abundance across milk kefir samples
top_pathways_df <- hm3_ph |>
  group_by(pathway) |>
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = 'drop') |>
  arrange(desc(mean_abundance))

top_pathways_df |> head(10)

# Select top pathways manually
top_pathways <- c(
  'LACTOSECAT-PWY: lactose and galactose degradation I',
  'PWY0-1586: peptidoglycan maturation (meso-diaminopimelate containing)',
  'PWY-7221: guanosine ribonucleotides de novo biosynthesis',
  'CALVIN-PWY: Calvin-Benson-Bassham cycle',
  'ANAGLYCOLYSIS-PWY: glycolysis III (from glucose)',
  'PWY-5973: cis-vaccenate biosynthesis',
  'PWY0-1296: purine ribonucleosides degradation',
  'PWY-7663: gondoate biosynthesis (anaerobic)',
  'PWY-6122: 5-aminoimidazole ribonucleotide biosynthesis II',
  'PWY-6277: superpathway of 5-aminoimidazole ribonucleotide biosynthesis'
)

# Filter selected pathways
hm3_filtered <- hm3_ph |>
  filter(pathway %in% top_pathways)

# Run Spearman correlation
hm3_cor_results <- hm3_filtered |>
  group_by(pathway) |>
  summarise(
    rho = cor(abundance, pH, method = 'spearman'),
    p_value = cor.test(abundance, pH, method = 'spearman', exact = FALSE)$p.value,
    .groups = 'drop'
  ) |>
  mutate(
    p_adj = p.adjust(p_value, method = 'BH')
  ) |>
  arrange(p_adj)

hm3_cor_results

# Filter only significant pathways
sig_pathways <- hm3_cor_results |>
  filter(p_adj < 0.05) |>
  pull(pathway)

hm3_sig <- hm3_filtered |>
  filter(pathway %in% sig_pathways) |>
  mutate(
    pathway_short = str_extract(pathway, '^[^:]+')
  )

# Manually place correlation labels for significant pathways
hm3_labels_manual <- tibble(
  pathway = c(
    'PWY0-1586: peptidoglycan maturation (meso-diaminopimelate containing)',
    'PWY0-1296: purine ribonucleosides degradation',
    'PWY-7221: guanosine ribonucleotides de novo biosynthesis',
    'CALVIN-PWY: Calvin-Benson-Bassham cycle'
  ),
  pathway_short = c(
    'PWY0-1586',
    'PWY0-1296',
    'PWY-7221',
    'CALVIN-PWY'
  ),
  x = c(1870, 1780, 1630, 1470),   # adjust these manually
  y = c(6.5, 3.85, 3.85, 6.5),       # adjust these manually
  label = c(
    'rho = -0.77\nadj. p = 0.005',
    'rho = 0.69\nadj. p = 0.017',
    'rho = 0.65\nadj. p = 0.023',
    'rho = -0.61\nadj. p = 0.029'
  )
) |>
  mutate(
    pathway_short = factor(pathway_short, levels = unique(hm3_sig$pathway_short))
  )

# Plot significant pathways
hm3_cor_plot <- ggplot(hm3_sig, aes(x = abundance, y = pH)) +
  geom_smooth(
    method = 'lm',
    color = '#E8403D',
    fill = '#F5CBCB',
    se = TRUE) +
  geom_point(color = '#4FACAC', size = 3) +
  facet_wrap(~ pathway_short, scales = 'free_x') +
  geom_text(
    data = hm3_labels_manual,
    mapping = aes(x = x, y = y, label = label),
    hjust = 1, # alignment right (left to right)
    vjust = 1, # alignment right (top to bottom)
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    title = '(B)',
    subtitle = 'Functional',
    x = 'Pathway Abundance (CPM)',
    y = 'pH'
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    plot.subtitle = element_text(face = 'bold', size = 12, hjust = 0.5),
    strip.text = element_text(face = 'bold', size = 12),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12)
  )

# Check the plot
hm3_cor_plot

# Save the plot
ggsave(
  hm3_cor_plot,
  filename = 'plots/05_functional_profiling/humann3_latest/hm3_pathway_correlation_pH.tiff',
  height = 20,
  width = 30,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine Taxonomic and Functional Correlation Analysis -------------------

library(gridExtra)

tax_func_combined <- grid.arrange(mp4_cor_plot, hm3_cor_plot, nrow = 2, heights = c(1, 1.8))

# Save the plot
ggsave(
  tax_func_combined,
  filename = 'plots/05_functional_profiling/humann3_latest/combined_correlation_analysis_plot.tiff',
  height = 45,
  width = 30,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)