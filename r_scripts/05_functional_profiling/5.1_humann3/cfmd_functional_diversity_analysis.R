# cFMD vs Caucasus Diversity Analysis
# ~/Documents/RStudio/mk_caucasus/scripts/05_functional_profiling/5.1_humann3
# 05/02/2026

# Load the necessary packages.
library(ggrepel)
library(patchwork)
library(tidyverse)
library(vegan)


# Load and Check Data -----------------------------------------------------

# Load caucasus and cFMD metadata
hm3_caucasus_metadata <- read_csv(
  'data/05_functional_profiling/humann3/hm3_caucasus_kefir_clean_metadata.csv',
  show_col_types = FALSE
) |>
  filter(sample_type == 'milk kefir')

hm3_cfmd_metadata <- read_tsv(
  'data/05_functional_profiling/humann3_latest/cfmd/cfmd_metadata_latest.tsv',
  show_col_types = FALSE
) |>
  filter(type == 'kefir', subtype == 'milk_kefir') |>
  mutate(subtype = case_when(
    subtype == 'milk_kefir' ~ 'milk kefir'
  )) |>
  rename(sample_type = subtype, country_code = country, sample = sample_id) |>
  select(sample, dataset_name, sample_type, country_code) |>
  filter(!sample %in% c('T12', 'U1', 'U3', 'U6', # samples with 0 metaphlan output
                        'T11', 'U2', 'U5', # super outlier samples (~20000) in metaphlan
                        'U4', # sample with 0 humann output
                        'Y4')) # super outlier sample (~20000) in humann

# Merge both metadata
hm3_merged_mk_metadata <- bind_rows(hm3_caucasus_metadata, hm3_cfmd_metadata)

# Save the merged functional profile so you don't have to run the above commands.
write_tsv(hm3_merged_mk_metadata, 'data/05_functional_profiling/humann3_latest/cfmd/merged_mk_metadata_strat.tsv')

# Get the milk kefir sample names in hm3_cfmd_metadata to filter samples needed in humann3 output in cfmd 
hm3_sample_names <- hm3_cfmd_metadata$sample

# Load caucasus and cFMD functional profiles and select only milk kefir samples
hm3_ancla_strat <- read_tsv(
  'data/05_functional_profiling/humann3_latest/caucasus_stratified_pathabundance_cpm/hm3_caucasus_joined_pathabundance_cpm_stratified.tsv',
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_Abundance-CPM$')) |> # remove suffix _Abundance-CPM in all column names
  filter(!str_detect(pathway, 'UNINTEGRATED')) |> # removing UNINTEGRATED rows in pathway column
  select(-matches('_g$')) |> # removing kefir grains
  filter(rowSums(across(-pathway)) > 0) |> # removes pathways with zero values
  separate(
    pathway,
    into = c('pathway', 'taxonomy'), sep = '\\|', remove = TRUE
  ) |> # split pathway and taxonomy
  filter(!grepl('unclassified', taxonomy)) |> # removing unclassified taxa in taxonomy column
  mutate(
    genus = str_extract(taxonomy, '(?<=g__)[^\\.]+'),
    species = str_extract(taxonomy, '(?<=s__).*'),
    species = str_replace_all(species, '_', ' ')
  ) |> # extract genus and species
  select(-taxonomy, -genus) # drop the raw taxonomy and genus column if not needed

hm3_gonzalezorozco_strat <- read_tsv(
  'data/05_functional_profiling/humann3_latest/cfmd/GonzalezOrozcoB_2023_stratified_pathabundance_cpm/hm3_GonzalezOrozcoB_2023_joined_pathabundance_cpm_stratified.tsv',
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_pathabundance_cpm$')) |> # remove suffix '_pathabundance_cpm' in all column names
  filter(!str_detect(pathway, 'UNINTEGRATED')) |> # removing UNINTEGRATED rows in pathway column
  select(pathway, any_of(hm3_sample_names)) |> # selecting columns present in metadata
  filter(rowSums(across(-pathway)) > 0) |> # removes pathways with zero values
  separate(
    pathway,
    into = c('pathway', 'taxonomy'), sep = '\\|', remove = TRUE
  ) |> # split pathway and taxonomy
  filter(!grepl('unclassified', taxonomy)) |> # removing unclassified taxa in taxonomy column
  mutate(
    genus = str_extract(taxonomy, '(?<=g__)[^\\.]+'),
    species = str_extract(taxonomy, '(?<=s__).*'),
    species = str_replace_all(species, '_', ' ')
  ) |> # extract genus and species
  select(-taxonomy, -genus) # drop the raw taxonomy and genus column if not needed

hm3_leech_strat <- read_tsv(
  'data/05_functional_profiling/humann3_latest/cfmd/LeechJ_2020_stratified_pathabundance_cpm/hm3_LeechJ_2020_joined_pathabundance_cpm_stratified.tsv',
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_pathabundance_cpm$')) |> # remove suffix '_pathabundance_cpm' in all column names
  filter(!str_detect(pathway, 'UNINTEGRATED')) |> # removing UNINTEGRATED rows in pathway column
  select(pathway, any_of(hm3_sample_names)) |> # selecting columns present in metadata
  filter(rowSums(across(-pathway)) > 0) |> # removes pathways with zero values
  separate(
    pathway,
    into = c('pathway', 'taxonomy'), sep = '\\|', remove = TRUE
  ) |> # split pathway and taxonomy
  filter(!grepl('unclassified', taxonomy)) |> # removing unclassified taxa in taxonomy column
  mutate(
    genus = str_extract(taxonomy, '(?<=g__)[^\\.]+'),
    species = str_extract(taxonomy, '(?<=s__).*'),
    species = str_replace_all(species, '_', ' ')
  ) |> # extract genus and species
  select(-taxonomy, -genus) # drop the raw taxonomy and genus column if not needed

hm3_salgado_strat <- read_tsv(
  'data/05_functional_profiling/humann3_latest/cfmd/SalgadoTS_2021_stratified_pathabundance_cpm/hm3_SalgadoTS_2021_joined_pathabundance_cpm_stratified.tsv',
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_Abundance-CPM$')) |> # remove suffix '_Abundance-CPM' in all column names
  filter(!str_detect(pathway, 'UNINTEGRATED')) |> # removing UNINTEGRATED rows in pathway column
  select(pathway, any_of(hm3_sample_names)) |> # selecting columns present in metadata
  filter(rowSums(across(-pathway)) > 0) |> # removes pathways with zero values
  separate(
    pathway,
    into = c('pathway', 'taxonomy'), sep = '\\|', remove = TRUE
  ) |> # split pathway and taxonomy
  filter(!grepl('unclassified', taxonomy)) |> # removing unclassified taxa in taxonomy column
  mutate(
    genus = str_extract(taxonomy, '(?<=g__)[^\\.]+'),
    species = str_extract(taxonomy, '(?<=s__).*'),
    species = str_replace_all(species, '_', ' ')
  ) |> # extract genus and species
  select(-taxonomy, -genus) # drop the raw taxonomy and genus column if not needed

hm3_walshAM_strat <- read_tsv(
  'data/05_functional_profiling/humann3_latest/cfmd/WalshAM_2016_stratified_pathabundance_cpm/hm3_WalshAM_2016_joined_pathabundance_cpm_stratified.tsv',
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_Abundance-CPM$')) |> # remove suffix '_Abundance-CPM' in all column names
  filter(!str_detect(pathway, 'UNINTEGRATED')) |> # removing UNINTEGRATED rows in pathway column
  select(pathway, any_of(hm3_sample_names)) |> # selecting columns present in metadata
  filter(rowSums(across(-pathway)) > 0) |> # removes pathways with zero values
  separate(
    pathway,
    into = c('pathway', 'taxonomy'), sep = '\\|', remove = TRUE
  ) |> # split pathway and taxonomy
  filter(!grepl('unclassified', taxonomy)) |> # removing unclassified taxa in taxonomy column
  mutate(
    genus = str_extract(taxonomy, '(?<=g__)[^\\.]+'),
    species = str_extract(taxonomy, '(?<=s__).*'),
    species = str_replace_all(species, '_', ' ')
  ) |> # extract genus and species
  select(-taxonomy, -genus) # drop the raw taxonomy and genus column if not needed

hm3_walshL_strat <- read_tsv(
  'data/05_functional_profiling/humann3_latest/cfmd/WalshL_xxxx_stratified_pathabundance_cpm/hm3_WalshL_xxxx_joined_pathabundance_cpm_stratified.tsv',
  show_col_types = FALSE) |>
  rename_with(~ str_remove(., '_Abundance-CPM$')) |> # remove suffix '_Abundance-CPM' in all column names
  filter(!str_detect(pathway, 'UNINTEGRATED')) |> # removing UNINTEGRATED rows in pathway column
  select(pathway, any_of(hm3_sample_names)) |> # selecting columns present in metadata
  filter(rowSums(across(-pathway)) > 0) |> # removes pathways with zero values
  separate(
    pathway,
    into = c('pathway', 'taxonomy'), sep = '\\|', remove = TRUE
  ) |> # split pathway and taxonomy
  filter(!grepl('unclassified', taxonomy)) |> # removing unclassified taxa in taxonomy column
  mutate(
    genus = str_extract(taxonomy, '(?<=g__)[^\\.]+'),
    species = str_extract(taxonomy, '(?<=s__).*'),
    species = str_replace_all(species, '_', ' ')
  ) |> # extract genus and species
  select(-taxonomy, -genus) # drop the raw taxonomy and genus column if not needed


# Build a Pathway-only Abundance Matrix (for Bray–Curtis) -----------------

# Function to prepare pathway-only matrix
hm3_prep_pathway_matrix <- function(df) {
  df |>
    select(-species) |> # drop species
    group_by(pathway) |> # sum across species
    summarise(across(where(is.numeric), sum), .groups = 'drop') |>
    column_to_rownames('pathway') |>
    t() |>
    as.data.frame()
}

# Applying function hm3_prep_pathway_matrix to each dataset
hm3_ancla_path_mat <- hm3_prep_pathway_matrix(hm3_ancla_strat)
hm3_gonzalezorozco_path_mat <- hm3_prep_pathway_matrix(hm3_gonzalezorozco_strat)
hm3_leech_path_mat <- hm3_prep_pathway_matrix(hm3_leech_strat)
hm3_salgado_path_mat <- hm3_prep_pathway_matrix(hm3_salgado_strat)
hm3_walshAM_path_mat <- hm3_prep_pathway_matrix(hm3_walshAM_strat)
hm3_walshL_path_mat <- hm3_prep_pathway_matrix(hm3_walshL_strat)

# Combine all of the matrix path profiles
hm3_merged_path_mat <- bind_rows(
  hm3_ancla_path_mat,
  hm3_gonzalezorozco_path_mat,
  hm3_leech_path_mat,
  hm3_salgado_path_mat,
  hm3_walshAM_path_mat,
  hm3_walshL_path_mat
)

hm3_merged_path_mat[is.na(hm3_merged_path_mat)] <- 0

# Save the hm3_merged_path_matrix so you don't have to run the above commands.
write_tsv(hm3_merged_path_mat, 'data/05_functional_profiling/humann3_latest/cfmd/merged_path_matrix.tsv')


# Build a Species Contribution Matrix (for envfit) ------------------------

# Function to prepare species-only matrix
hm3_prep_species_envfit <- function(df) {
  df |>
    group_by(species) |>
    summarise(across(where(is.numeric), sum), .groups = 'drop') |>
    column_to_rownames('species') |>
    t() |>
    as.data.frame()
}

# Applying function hm3_prep_species_envfit to each dataset
hm3_ancla_species_env <- hm3_prep_species_envfit(hm3_ancla_strat)
hm3_gonzalezorozco_species_env <- hm3_prep_species_envfit(hm3_gonzalezorozco_strat)
hm3_leech_species_env <- hm3_prep_species_envfit(hm3_leech_strat)
hm3_salgado_species_env <- hm3_prep_species_envfit(hm3_salgado_strat)
hm3_walshAM_species_env <- hm3_prep_species_envfit(hm3_walshAM_strat)
hm3_walshL_species_env <- hm3_prep_species_envfit(hm3_walshL_strat)

# Merge all species_env (union of species columns); NAs -> 0
hm3_species_env <- bind_rows(
  hm3_ancla_species_env,
  hm3_gonzalezorozco_species_env,
  hm3_leech_species_env,
  hm3_salgado_species_env,
  hm3_walshAM_species_env,
  hm3_walshL_species_env
)

hm3_species_env[is.na(hm3_species_env)] <- 0


# Milk Kefir Beta Diversity and Ordination --------------------------------

# Confirm empty samples before Bray-Curtis
hm3_row_sums <- rowSums(hm3_merged_path_mat)

# Identify the empty samples (sum == 0)
hm3_empty_samples <- names(hm3_row_sums[hm3_row_sums == 0])
length(hm3_empty_samples)
hm3_empty_samples

# Remove empty samples before Bray–Curtis
hm3_merged_path_mat_filt <- hm3_merged_path_mat[rowSums(hm3_merged_path_mat) > 0, , drop = FALSE]

# Keep metadata + species_env aligned to the filtered samples
hm3_merged_mk_metadata_filt <- hm3_merged_mk_metadata |>
  filter(sample %in% rownames(hm3_merged_path_mat_filt))

hm3_species_env_filt <- hm3_species_env[rownames(hm3_merged_path_mat_filt), , drop = FALSE]

# Bray-Curtis dissimilarity
hm3_bray_dist <- vegdist(hm3_merged_path_mat_filt, method = 'bray')


## NMDS Ordination --------------------------------------------------------

# Non-Metric Multidimensional Scaling (NMDS)
hm3_nmds_ctry <- metaMDS(hm3_bray_dist, k = 2, trymax = 100)
# k = 2: sets the number of dimensions (axes) / 2D space
# trymax = 100: sets the max number of random starts (iteration) that the algorithm will try to find a stable solution
# Try up to 100 random configurations to fine the best (lowest stress) solution

# Format results for plotting using ggplot
hm3_nmds_df_ctry <- as.data.frame(hm3_nmds_ctry$points) # converting to data frame
hm3_nmds_df_ctry$sample <- rownames(hm3_nmds_df_ctry) # adding sample column
hm3_nmds_df_ctry <- hm3_nmds_df_ctry |>
  left_join(hm3_merged_mk_metadata |> select(sample, country_code), by = 'sample')
hm3_nmds_df_ctry$shape <- ifelse(grepl('_m$', hm3_nmds_df_ctry$sample), 'Caucasus', 'cFMD')

# Save hm3_nmds_df_ctry
write.csv(hm3_nmds_df_ctry,
  'data/results/05_functional_profiling/humann3_latest/cfmd/mk_nmds_country_strat.csv',
  row.names = TRUE)


# Plotting with NO ellipse ------------------------------------------------

# Create the label text manually based on PERMANOVA result
hm3_ctry_nmds_perm_label_a <- data.frame(
  MDS1 = min(hm3_nmds_df_ctry$MDS1) + 0.02,
  MDS2 = max(hm3_nmds_df_ctry$MDS2) - 0.02,
  label = "bold('PERMANOVA')"
)

hm3_ctry_nmds_r2p_label_a <- data.frame(
  MDS1 = min(hm3_nmds_df_ctry$MDS1) + 0.02,
  MDS2 = max(hm3_nmds_df_ctry$MDS2) - 0.08,
  label = "italic(R^2)*' = 0.3, '*italic(p)*' = 0.001'"
)

hm3_distinct_cols <- c(
  '#E41A1C', # ALA
  '#377EB8', # AUS
  '#4DAF4A', # BEL
  '#984EA3', # BGR
  '#FF4F0F', # CAN
  '#016B61', # CHN
  '#A65628', # DEU
  '#F781BF', # DNK
  '#999999', # ESP
  '#B2D8CE', # FRA
  '#FC8D62', # GBR
  '#8DA0CB', # GRC
  '#F5CBCB', # HRV
  '#A6D854', # IND
  '#1B9E77', # IRL
  '#FFD92F', # ISR
  '#E5C494', # ITA
  '#0BA6DF', # MEX
  '#D95F02', # NZL
  '#7570B3', # PRT
  '#E7298A', # RUS
  '#0A400C', # SGP
  '#E1AA36', # SWE
  '#A6761D', # TUN
  '#640D5F', # TWN
  '#00809D', # USA
  '#8A784E'  # ZAF
)

# Plot the data
hm3_nmds_plot_ctry_a <- ggplot(hm3_nmds_df_ctry, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = country_code, shape = shape), size = 3) +
  scale_shape_manual(values = c('Caucasus' = 17, 'cFMD' = 16)) +
  # scale_fill_viridis_d(option = 'F') + # for ellipses
  # scale_color_viridis_d(option = 'F') + # for points
  scale_fill_manual(values = hm3_distinct_cols, name = 'Countries') + # for ellipses
  scale_color_manual(values = hm3_distinct_cols, name = 'Countries') + # for points
  labs(
    # title = 'Bray-Curtis - NMDS Ordination by Countries (HUMAnN3)', 
    title = '(B)',
    x = 'NMDS 1',
    y = 'NMDS 2',
    color = 'Countries',
    fill = 'Countries',
    shape = 'Datasets'
  ) +
  # label *only* Caucasus samples
  geom_text(
    data = subset(hm3_nmds_df_ctry, shape == 'Caucasus'),
    aes(label = sample),
    hjust = -0.4,
    vjust = -0.05,
    size = 3
  ) +
  geom_text(
    data = hm3_ctry_nmds_perm_label_a,
    aes(x = MDS1, y = MDS2, label = label),
    parse = TRUE, hjust = 0, vjust = 1, size = 5
  ) +
  geom_text(
    data = hm3_ctry_nmds_r2p_label_a,
    aes(x = MDS1, y = MDS2, label = label),
    parse = TRUE, hjust = 0, vjust = 1, size = 5
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill  = guide_legend(order = 1, title.position = 'top'), # Countries on top
    color = guide_legend(order = 2, title.position = 'top',
                         override.aes = list(alpha = 1, size = 4)) # Datasets below
  )

# Check the plot
hm3_nmds_plot_ctry_a

# Save the plot
ggsave(
  hm3_nmds_plot_ctry_a,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/mk_nmds_country_noellipse_strat.png',
  height = 15,
  width = 22,
  units = 'cm',
  dpi = 300
)

# Save the plot
ggsave(
  hm3_nmds_plot_ctry_a,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/mk_nmds_country_noellipse_strat.tiff',
  height = 15,
  width = 22,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


## PERMANOVA --------------------------------------------------------------

hm3_permanova_ctry <- adonis2(hm3_bray_dist ~ country_code, data = hm3_nmds_df_ctry)
# PERMANOVA result is the same even when using mk_pcoa_df_ctry

# Convert to data frame
hm3_permanova_df_ctry <- as.data.frame(hm3_permanova_ctry)

# Write to CSV
write.csv(
  hm3_permanova_df_ctry,
  'data/results/05_functional_profiling/humann3_latest/cfmd/mk_permanova_country_strat.csv',
  row.names = TRUE)


## Vegan - Envfit (species drivers) --------------------------------------

# Optional but recommended for envfit stability on compositional data
hm3_species_env_hell <- decostand(hm3_species_env_filt, method = 'hellinger')

set.seed(1)
hm3_ef_species <- envfit(hm3_nmds_ctry, hm3_species_env_hell, permutations = 999)

# Extract vectors
hm3_sp_vec <- as.data.frame(hm3_ef_species$vectors$arrows)
hm3_sp_vec$r2 <- hm3_ef_species$vectors$r
hm3_sp_vec$p  <- hm3_ef_species$vectors$pvals
hm3_sp_vec$taxon <- rownames(hm3_sp_vec)
hm3_sp_vec$q <- p.adjust(hm3_sp_vec$p, method = 'BH')

# Scale arrows into ~30% of axis span
hm3_x_rng <- range(hm3_nmds_df_ctry$MDS1, na.rm = TRUE)
hm3_y_rng <- range(hm3_nmds_df_ctry$MDS2, na.rm = TRUE)
hm3_x_span <- diff(hm3_x_rng)
hm3_y_span <- diff(hm3_y_rng)

hm3_max_abs_x <- max(abs(hm3_sp_vec$NMDS1), na.rm = TRUE)
hm3_max_abs_y <- max(abs(hm3_sp_vec$NMDS2), na.rm = TRUE)

hm3_mul_x <- (hm3_x_span * 0.30) / hm3_max_abs_x
hm3_mul_y <- (hm3_y_span * 0.30) / hm3_max_abs_y
hm3_arrow_mul <- min(hm3_mul_x, hm3_mul_y)

hm3_sp_vec$Axis.1 <- hm3_sp_vec$NMDS1 * hm3_arrow_mul
hm3_sp_vec$Axis.2 <- hm3_sp_vec$NMDS2 * hm3_arrow_mul

# Keep significant species or fallback to top 12 by r²
hm3_sp_sig <- hm3_sp_vec |> filter(q <= 0.05) |> arrange(desc(r2))
if (nrow(hm3_sp_sig) == 0) {
  message('⚠️ No species significant at q <= 0.05. Showing top 12 by r².')
  hm3_sp_top <- hm3_sp_vec |> arrange(desc(r2)) |> slice_head(n = 12)
} else {
  hm3_sp_top <- hm3_sp_sig |> slice_head(n = 12)
}

# Species label formatter (italicize species; preserve ' group' suffix if present)
format_taxon_label <- function(taxon) {
  clean <- gsub('_', ' ', taxon)
  
  if (grepl(' group$', clean)) {
    main <- sub(' group$', '', clean)
    paste0('italic("', main, '")~plain(" group")')
  } else {
    paste0('italic("', clean, '")')
  }
}

hm3_sp_top$label <- vapply(hm3_sp_top$taxon, format_taxon_label, character(1))

# Plot NMDS with envfit arrows
hm3_nmds_plot_ctry_species <- hm3_nmds_plot_ctry_a +
  geom_segment(
    data = hm3_sp_top,
    aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
    arrow = arrow(length = unit(0.15, 'cm')),
    linewidth = 0.4,
    linetype = 'dashed'
  ) +
  geom_text_repel(
    data = hm3_sp_top,
    aes(x = Axis.1, y = Axis.2, label = label),
    parse = TRUE,
    size = 3,
    color = '#295F98',
    max.overlaps = 50,
    force = 2,
    box.padding = 0.3,
    segment.color = 'grey50'
  )

# Check the plot
hm3_nmds_plot_ctry_species

# Save the plot
ggsave(
  hm3_nmds_plot_ctry_species,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/mk_nmds_country_species.png',
  height = 20,
  width = 25,
  units = 'cm',
  dpi = 300
)

# Save the plot
ggsave(
  hm3_nmds_plot_ctry_species,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/mk_nmds_country_species.tiff',
  height = 20,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)

# Export CSV results

# Mark significant taxa
hm3_sp_vec <- hm3_sp_vec |>
  mutate(significant = ifelse(q <= 0.05, 'yes', 'no'))

# Save full list of taxa tested
write.csv(hm3_sp_vec,
  'data/results/05_functional_profiling/humann3_latest/cfmd/mk_envfit_species_all_NMDS.csv',
  row.names = TRUE)

# Save top taxa (the ones shown in the plot)
write.csv(hm3_sp_top,
  'data/results/05_functional_profiling/humann3_latest/cfmd/mk_envfit_species_top12_NMDS.csv',
  row.names = TRUE)


## Distance of Countries from Caucasus using Bray–Curtis ------------------

library(countrycode)

# --- helper: ISO2 -> flag emoji (e.g., 'es' -> '🇪🇸') ---
hm3_flag_from_iso2 <- function(iso2) {
  vapply(iso2, function(x) {
    if (is.na(x) || nchar(x) != 2) return('')
    u <- utf8ToInt(toupper(x))
    intToUtf8(u + 127397)
  }, character(1))
}

# Bray–Curtis distance matrix
hm3_bray_mat <- as.matrix(hm3_bray_dist)

# Align metadata to bray_mat row order (critical)
hm3_meta_aligned <- hm3_merged_mk_metadata_filt |>
  filter(sample %in% rownames(hm3_bray_mat)) |>
  distinct(sample, country_code) |>
  slice(match(rownames(hm3_bray_mat), sample))

# Identify Caucasus samples (your rule: end with '_m')
hm3_cauc_samples <- hm3_meta_aligned$sample[grepl('_m$', hm3_meta_aligned$sample)]
hm3_cfmd_samples <- setdiff(hm3_meta_aligned$sample, hm3_cauc_samples)

# Mean Bray–Curtis distance of each cFMD sample to all Caucasus samples
hm3_cfmd_mean_to_cauc <- tibble(sample = hm3_cfmd_samples) |>
  mutate(mean_to_cauc = rowMeans(
    hm3_bray_mat[sample, hm3_cauc_samples, drop = FALSE],
    na.rm = TRUE
  ))

# Summarise per country & top 10 closest (lowest dissimilarity)
hm3_top10_countries_bray <- hm3_cfmd_mean_to_cauc |>
  left_join(hm3_meta_aligned, by = 'sample') |>
  group_by(country_code) |>
  summarise(
    mean_distance = mean(mean_to_cauc, na.rm = TRUE),
    n = dplyr::n(),
    .groups = 'drop'
  ) |>
  arrange(mean_distance) |>
  slice_head(n = 10) |>
  mutate(
    country_name = countrycode(country_code, 'iso3c', 'country.name'),
    iso2 = tolower(countrycode(country_code, 'iso3c', 'iso2c')),
    flag = hm3_flag_from_iso2(iso2),
    axis_label = paste0(country_name, '  ', flag)
  )

# Named vector for axis labels
hm3_axis_labs <- setNames(hm3_top10_countries_bray$axis_label, hm3_top10_countries_bray$country_name)

# Ranked bar plot (same style as your NMDS version)
hm3_bray_plot_distance <- ggplot(
  hm3_top10_countries_bray,
  aes(
    x = reorder(country_name, -mean_distance),
    y = mean_distance
  )
) +
  geom_col(width = 0.7, fill = '#4FACAC') +
  geom_text(aes(label = round(mean_distance, 3)),
            hjust = -0.1, size = 4) +
  scale_x_discrete(labels = hm3_axis_labs) +
  coord_flip() +
  labs(
    title = '(D)',
    x = 'Countries',
    y = 'Mean Bray-Curtis dissimilarity to Caucasus samples'
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  expand_limits(y = max(hm3_top10_countries_bray$mean_distance, na.rm = TRUE) * 1.1)

# Check the plot
hm3_bray_plot_distance

# Save the plot
ggsave(
  hm3_bray_plot_distance,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/mk_bray_country_distance.tiff',
  height = 15,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine Taxonomic and Functional Data -----------------------------------

library(gridExtra)

layout_1 <- rbind(
  c(1, 2),
  c(3, 4)
)

cfmd_tax_func_combined_plot <- grid.arrange(
  mp4_nmds_plot_ctry_taxa,
  mp4_bray_plot_distance,
  hm3_nmds_plot_ctry_species,
  hm3_bray_plot_distance,
  layout_matrix = layout_1,
  heights = c(1, 1),
  widths = c(1.4, 1)
)

# Save the plot
ggsave(
  cfmd_tax_func_combined_plot,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/cfmd_tax_func_combined_plot.tiff',
  height = 40,
  width = 45,
  units = 'cm',
  bg = 'white',
  dpi = 600,
  compression = 'lzw'
)


# Combine Beta Diversity and Distance -------------------------------------

hm3_cfmd_combined_plot <- (hm3_nmds_plot_ctry_species | hm3_bray_plot_distance) +
  plot_layout(widths = c(1, 0.7))

# Check the plot
hm3_cfmd_combined_plot

# Save the plot
ggsave(
  hm3_cfmd_combined_plot,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/cfmd_combined_plot_final.png',
  height = 18,
  width = 44,
  units = 'cm',
  dpi = 300
)


# Final Combination of Taxonomy and Function Caucasus vs cFMD -------------

final_cfmd_combined_plot <- (mp4_cfmd_combined_plot / hm3_cfmd_combined_plot)

# Check the plot
final_cfmd_combined_plot

# Save the plot
ggsave(final_cfmd_combined_plot,
       filename = 'plots/05_functional_profiling/humann3_latest/cfmd/final_cfmd_combined_plot.png',
       height = 40,
       width = 45,
       units = 'cm',
       dpi = 300
)


# Dispersion Analysis -----------------------------------------------------

# Functional dispersion by dataset (Caucasus vs cFMD)
hm3_bd_dataset <- betadisper(hm3_bray_dist, group = hm3_nmds_df_ctry$shape)
permutest(hm3_bd_dataset, permutations = 999)

# Plot distances to centroid (within-group dispersion)
# NOTE: betadisper runs a PCoA internally (on distances-to-centroid space),
# so the axes *do* have eigenvalues and % variance explained.

# Extract site scores (PCoA coordinates) and centroids
hm3_sites <- as.data.frame(scores(hm3_bd_dataset, display = 'sites'))
hm3_sites$sample <- rownames(hm3_sites)
hm3_sites$dataset <- hm3_bd_dataset$group

hm3_centroids <- as.data.frame(scores(hm3_bd_dataset, display = 'centroids'))
hm3_centroids$dataset <- rownames(hm3_centroids)

# Add centroid coordinates per sample (for spider segments)
hm3_sites <- hm3_sites |>
  left_join(
    hm3_centroids |> rename(centroid_x = PCoA1, centroid_y = PCoA2),
    by = 'dataset'
  )

# % variance explained by PCoA 1 and 2 (HUMAnN / function)
# (use only positive eigenvalues for the denominator, consistent with PCoA reporting)
hm3_eig <- hm3_bd_dataset$eig
hm3_eig_pos <- hm3_eig[hm3_eig > 0]

hm3_pcoa1_pct <- round(100 * hm3_eig_pos[1] / sum(hm3_eig_pos), 1)
hm3_pcoa2_pct <- round(100 * hm3_eig_pos[2] / sum(hm3_eig_pos), 1)

# Convex hull coordinates per dataset (to mimic base plot outline)
hm3_hulls <- hm3_sites |>
  group_by(dataset) |>
  slice(chull(PCoA1, PCoA2)) |>
  ungroup()

# Colors (as requested)
hm3_cols <- c('Caucasus' = '#E8403D', 'cFMD' = '#4FACAC')

# Create label text based on betadisper permutest result
hm3_bd_label_title <- data.frame(
  PCoA1 = - 0.5,
  PCoA2 = + 0.35,
  label = "bold('Dispersion Test')"
)

hm3_bd_label_stats <- data.frame(
  PCoA1 = - 0.5,
  PCoA2 = + 0.32,
  label = "italic(F)*' = 0.38, '*italic(p)*' = 0.58'"
)

# Build ggplot object
hm3_bd_plot <- ggplot(hm3_sites, aes(x = PCoA1, y = PCoA2)) +
  # Spider segments: sample -> centroid
  geom_segment(
    aes(x = centroid_x, y = centroid_y, xend = PCoA1, yend = PCoA2, color = dataset),
    alpha = 0.35,
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  # Hull outline per dataset
  geom_polygon(
    data = hm3_hulls,
    aes(x = PCoA1, y = PCoA2, group = dataset, color = dataset, fill = dataset),
    alpha = 0.08,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  # Points
  geom_point(aes(color = dataset, shape = dataset), size = 2.5, alpha = 0.9) +
  # Centroids
  geom_point(
    data = hm3_centroids,
    aes(x = PCoA1, y = PCoA2, color = dataset),
    size = 4,
    shape = 4,   # X marker for centroid
    stroke = 1.1,
    show.legend = FALSE
  ) +
  geom_text(
    data = hm3_bd_label_title,
    aes(x = PCoA1, y = PCoA2, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 4.5
  ) +
  geom_text(
    data = hm3_bd_label_stats,
    aes(x = PCoA1, y = PCoA2, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 4.5
  ) +
  # Label only Caucasus milk kefir samples
  geom_text_repel(
    data = subset(hm3_sites, dataset == 'Caucasus'),
    aes(label = sample),
    size = 2.6,
    fontface = 'plain',
    color = '#E8403D',
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = 'grey60',
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = hm3_cols, name = 'Datasets') +
  scale_fill_manual(values = hm3_cols) +
  scale_shape_manual(values = c('Caucasus' = 17, 'cFMD' = 16), name = 'Datasets') +
  labs(
    title = '(B)',
    x = paste0('PCoA 1 (', hm3_pcoa1_pct, '%)'),
    y = paste0('PCoA 2 (', hm3_pcoa2_pct, '%)')
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', size = 12),
    axis.title = element_text(face = 'bold', size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'right'
  )

# Check the plot
hm3_bd_plot

# Save the plot
ggsave(
  hm3_bd_plot,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/hm3_betadispers_plot_final.png',
  height = 15,
  width = 22,
  units = 'cm',
  dpi = 300
)

# Save the plot
ggsave(
  hm3_bd_plot,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/hm3_betadispers_plot_final.tiff',
  height = 15,
  width = 22,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine Beta Disperse Plot ----------------------------------------------

cfmd_tax_func_bd_combined_plot <- grid.arrange(
  mp4_bd_plot,
  hm3_bd_plot,
  ncol = 1,
  heights = c(1, 1)
)

# Save the plot
ggsave(
  cfmd_tax_func_bd_combined_plot,
  filename = 'plots/05_functional_profiling/humann3_latest/cfmd/cfmd_tax_func_bd_combined_plot.tiff',
  height = 30,
  width = 22,
  units = 'cm',
  bg = 'white',
  dpi = 600,
  compression = 'lzw'
)