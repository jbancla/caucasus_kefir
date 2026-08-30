# Abricate MAGs Antibiotic Resistance Genes Profile
# ~/Documents/RStudio/mk_caucasus/scripts/07_downstream_analysis/7.3_abricate
# 24/01/2026


# This is the Revised Working Script --------------------------------------

# Load the necessary packages.
library(paletteer)
library(tidyverse)


# Load and Check Data -----------------------------------------------------

# Adding collection sites layer (kab = Kabardino-Balkaria)
kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

# Load antibiotic resistance genes data
arg_hqmags_data <- read_tsv('data/07_downstream_analysis/7.3_abricate/hq_mags/abricate_card_summary.tab', show_col_types = FALSE) |> 
  rename(bin_name = `#FILE`, arg_counts = NUM_FOUND) |> 
  mutate(
    # drop the long prefix
    bin_name = str_replace(bin_name, '^/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/07_downstream_analysis/7.3_abricate/hq_mags/', ''),
    # drop the suffix
    bin_name = str_replace(bin_name, '_card\\.tab$', ''),
    sample_name = str_replace(bin_name, "_bin\\..*$", ""),
    sample_type = case_when(grepl('_g$', sample_name) ~ 'Kefir Grains', TRUE ~ 'Milk Kefir'),
    collection_site = if_else(sample_name %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia')) |> 
  filter(!sample_name %in% c('k12_g', 's1.2_m', 's7.2_m', 's9.2_m', 'sfrankenstein')) |> 
  filter(!bin_name %in% c('k14_m_bin.2.orig', 'k14_m_bin.3.strict')) |> 
  select(-Bifidobacteria_intrinsic_ileS_conferring_resistance_to_mupirocin, -Bifidobacterium_adolescentis_rpoB_conferring_resistance_to_rifampicin,
         -emrR, -tetS)

# Load GDTB-Tk to append species names
gtdb_hqmags_data <- read_tsv('data/07_downstream_analysis/7.1_gtdbtk/hq_mags/gtdb_hqmags_with_qc_data.tsv', show_col_types = FALSE)

gtdb_species <- gtdb_hqmags_data |> select(species)

# Append species column to arg_hqmags_data result
arg_hqmags_final_data <- bind_cols(arg_hqmags_data, gtdb_species)

metadata_columns <- c('bin_name', 'arg_counts', 'species', 'sample_name', 'sample_type', 'collection_site')
arg_columns <- setdiff(names(arg_hqmags_final_data), metadata_columns)

arg_hqmags_final_data <- arg_hqmags_final_data |>
  mutate(
    across(
      all_of(arg_columns),
      ~ {
        x <- as.character(.)
        x[x == '.'] <- '0' # turn dot placeholders into '0'
        x <- readr::parse_number(x) # robust numeric parse
        replace(x, is.na(x), 0) # any leftovers to 0
      }
    )
  )

# Order species by how many MAGs they have (most -> least)
species_order <- arg_hqmags_final_data |>
  distinct(species, bin_name) |>
  count(species, name = 'n_mags') |>
  arrange(desc(n_mags), species) |>
  pull(species)

# Long format for plotting
arg_hqmags_long <- arg_hqmags_final_data |>
  mutate(species = factor(species, levels = species_order)) |>
  pivot_longer(
    cols = all_of(arg_columns),
    names_to = 'arg',
    values_to = 'coverage') |>
  mutate(
    # keep coverage as numeric so the legend is continuous (0 -> 100)
    coverage = pmin(pmax(as.numeric(coverage), 0), 100),
    # order ARGs alphabetically, A at the top
    arg = factor(arg, levels = rev(sort(unique(arg))))
  )


# Plot the Data -----------------------------------------------------------

arg_plot <- ggplot(arg_hqmags_long,
                   aes(x = bin_name, y = arg, fill = coverage)) +
  geom_tile(color = 'white', linewidth = 0.2, width = 0.90, height = 0.90) +
  facet_grid(. ~ species, scales = 'free_x', space = 'free_x', labeller = labeller(species = function(x) sub(' ', '\n', x))) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_paletteer_c(
    'grDevices::RdPu',
    direction = -1, # reverse colors, not values
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    guide = guide_colorbar(
      title.position = 'top',
      title.hjust = 0.5,
      direction = 'vertical'
    )
  ) +
  labs(
    x = 'No. of HQ MAGs per Species',
    y = 'Antibiotic Resistance Genes',
    fill = 'Coverage (%)') +
  theme_minimal(base_size = 10) +
  theme(
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.y = element_text(size = 12, color = 'black'),
    axis.ticks.y = element_line(),
    # panel.grid = element_blank(),
    panel.spacing.x = unit(0.07, 'cm'),
    panel.border = element_rect(colour = 'grey40', fill = NA, linewidth = 0.5),
    strip.text.x = element_text(face = 'italic', size = 7, angle = 90, color = 'black'),
    # legend.position = 'right',
    legend.title = element_text(face = 'bold', size = 12, margin = margin(b = 15)),
    legend.text = element_text(size = 12)
  )

# Check the plot
arg_plot

# Save the plot
ggsave(
  arg_plot,
  filename = 'plots/07_downstream_analysis/7.3_abricate/abricate_arg_plot_v2.tiff',
  height = 11,
  width = 35,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)