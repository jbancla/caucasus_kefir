# cFMD vs Caucasus Diversity Analysis
# ~/Documents/RStudio/mk_caucasus/scripts/04_taxonomic_profiling/4.1_metaphlan4
# 28/08/2025

# Load the necessary packages.
library(countrycode)
library(ggrepel)
library(patchwork)
library(tidyverse)
library(vegan)


# Load and Check Data -----------------------------------------------------

# Load caucasus and cFMD metadata
mp4_caucasus_metadata <- read_csv(
  'data/04_taxonomic_profiling/metaphlan4/mp4_caucasus_kefir_clean_metadata.csv',
  show_col_types = FALSE
) |>
  filter(sample_type == 'milk kefir')

mp4_cfmd_metadata <- read_tsv(
  'data/04_taxonomic_profiling/metaphlan4_latest/cfmd/cfmd_metadata_latest.tsv',
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
mp4_merged_mk_metadata <- bind_rows(mp4_caucasus_metadata, mp4_cfmd_metadata)

# Save the merged_taxprofile so you don't have to run the above commands.
write_tsv(mp4_merged_mk_metadata, 'data/04_taxonomic_profiling/metaphlan4_latest/cfmd/merged_mk_metadata.tsv')

# Get the milk kefir sample names in mp4_cfmd_metadata to filter samples needed in metaphlan output in cfmd 
mp4_sample_names <- mp4_cfmd_metadata$sample

# Load caucasus and cFMD taxonomic profiles and select only milk kefir samples
mp4_ancla_taxprofile <- read_tsv(
  'data/04_taxonomic_profiling/metaphlan4_latest/mp4_caucasus_merged_species.tsv',
  show_col_types = FALSE) |>
  select(-matches('_g$')) |>
  filter(rowSums(across(-clade_name)) > 0)

mp4_ancla_taxprofile_tp <- mp4_ancla_taxprofile |>
  column_to_rownames(var = 'clade_name') |>
  t() |>
  as.data.frame()

mp4_gonzalezorozco_taxprofile <- read_tsv(
  'data/04_taxonomic_profiling/metaphlan4_latest/mp4_gonzalezorozco_merged_species.tsv',
  show_col_types = FALSE) |>
  select(clade_name, any_of(mp4_sample_names)) |>
  filter(rowSums(across(-clade_name)) > 0)

mp4_gonzalezorozco_taxprofile_tp <- mp4_gonzalezorozco_taxprofile |>
  column_to_rownames(var = 'clade_name') |>
  t() |>
  as.data.frame()

mp4_leech_taxprofile <- read_tsv(
  'data/04_taxonomic_profiling/metaphlan4_latest/mp4_leech_merged_species.tsv',
  show_col_types = FALSE) |>
  select(clade_name, any_of(mp4_sample_names)) |>
  filter(rowSums(across(-clade_name)) > 0)

mp4_leech_taxprofile_tp <- mp4_leech_taxprofile |>
  column_to_rownames(var = 'clade_name') |>
  t() |>
  as.data.frame()

mp4_salgado_taxprofile <- read_tsv(
  'data/04_taxonomic_profiling/metaphlan4_latest/mp4_salgado_merged_species.tsv',
  show_col_types = FALSE) |>
  select(clade_name, any_of(mp4_sample_names)) |>
  filter(rowSums(across(-clade_name)) > 0)

mp4_salgado_taxprofile_tp <- mp4_salgado_taxprofile |>
  column_to_rownames(var = 'clade_name') |>
  t() |>
  as.data.frame()

mp4_walshAM_taxprofile <- read_tsv(
  'data/04_taxonomic_profiling/metaphlan4_latest/mp4_walshAM_merged_species.tsv',
  show_col_types = FALSE) |>
  select(clade_name, any_of(mp4_sample_names)) |>
  filter(rowSums(across(-clade_name)) > 0)

mp4_walshAM_taxprofile_tp <- mp4_walshAM_taxprofile |>
  column_to_rownames(var = 'clade_name') |>
  t() |>
  as.data.frame()

mp4_walshL_taxprofile <- read_tsv(
  'data/04_taxonomic_profiling/metaphlan4_latest/mp4_walshL_merged_species.tsv',
  show_col_types = FALSE) |>
  select(clade_name, any_of(mp4_sample_names)) |>
  filter(rowSums(across(-clade_name)) > 0)

mp4_walshL_taxprofile_tp <- mp4_walshL_taxprofile |>
  column_to_rownames(var = 'clade_name') |>
  t() |>
  as.data.frame()

mp4_merged_taxprofile_tp <- bind_rows(
  mp4_ancla_taxprofile_tp,
  mp4_gonzalezorozco_taxprofile_tp,
  mp4_leech_taxprofile_tp,
  mp4_salgado_taxprofile_tp,
  mp4_walshAM_taxprofile_tp,
  mp4_walshL_taxprofile_tp
)

mp4_merged_taxprofile_tp[is.na(mp4_merged_taxprofile_tp)] <- 0

# Save the merged_taxprofile so you don't have to run the above commands.
write_tsv(mp4_merged_taxprofile_tp, 'data/04_taxonomic_profiling/metaphlan4_latest/cfmd/merged_taxprofile_tp.tsv')

# Use the mp4_merged_mk_metadata and the mp4_merged_taxprofile_tp for diversity analysis

# Check for empty rows
row_sums <- rowSums(mp4_merged_taxprofile_tp)
sum(row_sums == 0)
row_sums[row_sums == 0]

# Upon checking, there are 4 samples with missing metaphlan output in walshL (T12, U1, U3, and U6)
# Filter again the cfmd_metadata before joining with caucasus_metadata
# Because those 4 samples will cause a problem in Bray-Curtis later.

# Milk Kefir Beta Diversity and  Ordination -------------------------------

# Bray-Curtis dissimilarity
mp4_bray_dist <- vegdist(mp4_merged_taxprofile_tp, method = 'bray')


# Country -----------------------------------------------------------------

## NMDS Ordination --------------------------------------------------------

# Non-Metric Multidimensional Scaling (NMDS)
mp4_nmds_ctry <- metaMDS(mp4_bray_dist, k = 2, trymax = 100)
# k = 2: sets the number of dimensions (axes) / 2D space
# trymax = 100: sets the max number of random starts (iteration) that the algorithm will try to find a stable solution
# Try up to 100 random configurations to fine the best (lowest stress) solution

# Format results for plotting using ggplot
mp4_nmds_df_ctry <- as.data.frame(mp4_nmds_ctry$points) # converting to data frame
mp4_nmds_df_ctry$sample <- rownames(mp4_nmds_df_ctry) # adding sample column
mp4_nmds_df_ctry <- mp4_nmds_df_ctry |>
  left_join(mp4_merged_mk_metadata |> select(sample, country_code), by = 'sample')
mp4_nmds_df_ctry$shape <- ifelse(grepl('_m$', mp4_nmds_df_ctry$sample), 'Caucasus', 'cFMD')

# Save mp4_nmds_df_ctry
write.csv(mp4_nmds_df_ctry,
  'data/results/04_taxonomic_profiling/metaphlan4_latest/cfmd/mk_nmds_country.csv',
  row.names = TRUE)

# Plotting with NO ellipse

mp4_distinct_cols <- c(
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

# Create the label text manually based on PERMANOVA result
mp4_ctry_nmds_perm_label_a <- data.frame(
  MDS1 = min(mp4_nmds_df_ctry$MDS1) + 0.02,
  MDS2 = max(mp4_nmds_df_ctry$MDS2) - 0.06,
  label = "bold('PERMANOVA')"
)

mp4_ctry_nmds_r2p_label_a <- data.frame(
  MDS1 = min(mp4_nmds_df_ctry$MDS1) + 0.02,
  MDS2 = max(mp4_nmds_df_ctry$MDS2) - 0.17,
  label = "italic(R^2)*' = 0.44, '*italic(p)*' = 0.001'"
)

# Plot the data
mp4_nmds_plot_ctry_a <- ggplot(mp4_nmds_df_ctry, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = country_code, shape = shape), size = 3) +
  scale_shape_manual(values = c('Caucasus' = 17, 'cFMD' = 16)) +
  # scale_fill_viridis_d(option = 'F') + # for ellipses
  # scale_color_viridis_d(option = 'F') + # for points
  scale_fill_manual(values = mp4_distinct_cols, name = 'Countries') + # for ellipses
  scale_color_manual(values = mp4_distinct_cols, name = 'Countries') + # for points
  labs(
    # title = 'Bray-Curtis - NMDS Ordination by Countries (MetaPhlAn4)', 
    title = '(A)',
    x = 'NMDS 1',
    y = 'NMDS 2',
    color = 'Countries',
    fill = 'Countries',
    shape = 'Datasets'
  ) +
  # label *only* Caucasus samples
  geom_text(
    data = subset(mp4_nmds_df_ctry, shape == 'Caucasus'),
    aes(label = sample),
    hjust = -0.4,
    vjust = -0.05,
    size = 3
  ) +
  geom_text(
    data = mp4_ctry_nmds_perm_label_a,
    aes(x = MDS1, y = MDS2, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  geom_text(
    data = mp4_ctry_nmds_r2p_label_a,
    aes(x = MDS1, y = MDS2, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 5
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
    fill = guide_legend(order = 1, title.position = 'top'), # Countries on top
    color = guide_legend(order = 2, title.position = 'top',
                         override.aes = list(alpha = 1, size = 4)) # Datasets below
  )

# Check the plot
mp4_nmds_plot_ctry_a

# Save the plot
ggsave(
  mp4_nmds_plot_ctry_a,
  filename = 'plots/04_taxonomic_profiling/metaphlan4_latest/cfmd/mk_nmds_country_noellipse.tiff',
  height = 15,
  width = 22,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


## PERMANOVA --------------------------------------------------------------

mp4_permanova_ctry <- adonis2(mp4_bray_dist ~ country_code, data = mp4_nmds_df_ctry)
# PERMANOVA result is the same even when using mk_pcoa_df_ctry

# Convert to data frame
mp4_permanova_df_ctry <- as.data.frame(mp4_permanova_ctry)

# Write to CSV
write.csv(mp4_permanova_df_ctry,
  'data/results/04_taxonomic_profiling/metaphlan4_latest/cfmd/mk_permanova_country.csv',
  row.names = TRUE)


## Vegan - Envfit ---------------------------------------------------------

# envfit with species (NMDS drivers)

# Hellinger transform your table (optional but good for NMDS stability)
mp4_taxa_mat_hell <- decostand(mp4_merged_taxprofile_tp, method = 'hellinger')

# Fit species to NMDS ordination
set.seed(1)
mp4_ef_taxa <- envfit(mp4_nmds_ctry, mp4_taxa_mat_hell, permutations = 999)

# Extract vectors
mp4_tax_vec <- as.data.frame(mp4_ef_taxa$vectors$arrows)
mp4_tax_vec$r2 <- mp4_ef_taxa$vectors$r
mp4_tax_vec$p  <- mp4_ef_taxa$vectors$pvals
mp4_tax_vec$taxon <- rownames(mp4_tax_vec)
mp4_tax_vec$q <- p.adjust(mp4_tax_vec$p, method = 'BH')

# Axis ranges (for scaling arrows)
mp4_x_rng <- range(mp4_nmds_df_ctry$MDS1, na.rm = TRUE)
mp4_y_rng <- range(mp4_nmds_df_ctry$MDS2, na.rm = TRUE)
mp4_x_span <- diff(mp4_x_rng)
mp4_y_span <- diff(mp4_y_rng)

# Scale arrows into ~30% of axis span
mp4_max_abs_x <- max(abs(mp4_tax_vec$NMDS1), na.rm = TRUE)
mp4_max_abs_y <- max(abs(mp4_tax_vec$NMDS2), na.rm = TRUE)
mp4_mul_x <- (mp4_x_span * 0.30) / mp4_max_abs_x
mp4_mul_y <- (mp4_y_span * 0.30) / mp4_max_abs_y
mp4_arrow_mul <- min(mp4_mul_x, mp4_mul_y)

mp4_tax_vec$Axis.1 <- mp4_tax_vec$NMDS1 * mp4_arrow_mul
mp4_tax_vec$Axis.2 <- mp4_tax_vec$NMDS2 * mp4_arrow_mul

# Keep significant taxa or fallback to top 10 by r²
mp4_tax_sig <- mp4_tax_vec |> filter(q <= 0.05) |> arrange(desc(r2))
if (nrow(mp4_tax_sig) == 0) {
  message('⚠️ No taxa significant at q <= 0.05. Showing top 10 by r².')
  mp4_tax_top <- mp4_tax_vec |> arrange(desc(r2)) |> slice_head(n = 10)
} else {
  mp4_tax_top <- mp4_tax_sig |> slice_head(n = 10)
}

# Function to format names: italicize genus+species, leave suffix plain
format_taxon_label <- function(taxon) {
  clean <- gsub('_', ' ', taxon)
  
  # If it ends with " group"
  if (grepl(' group$', clean)) {
    main <- sub(' group$', '', clean)
    return(paste0('italic("', main, '")~plain(" group")'))
  } else {
    return(paste0('italic("', clean, '")'))
  }
}

# Apply function to all taxa
mp4_tax_top$label <- vapply(mp4_tax_top$taxon, format_taxon_label, character(1))

# Plot the data
mp4_nmds_plot_ctry_taxa <- mp4_nmds_plot_ctry_a +
  geom_segment(
    data = mp4_tax_top,
    aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
    arrow = arrow(length = unit(0.15, 'cm')),
    linewidth = 0.4,
    linetype = 'dashed'
  ) +
  geom_text_repel(
    data = mp4_tax_top,
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
mp4_nmds_plot_ctry_taxa

# Save the plot
ggsave(
  mp4_nmds_plot_ctry_taxa,
  filename = 'plots/04_taxonomic_profiling/metaphlan4_latest/cfmd/mk_nmds_country_taxa.tiff',
  height = 20,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)

# Export CSV results

# Mark significant taxa
mp4_tax_vec <- mp4_tax_vec |>
  mutate(significant = ifelse(q <= 0.05, 'yes', 'no'))

# Save full list of taxa tested
write.csv(mp4_tax_vec,
  'data/results/04_taxonomic_profiling/metaphlan4_latest/cfmd/mk_envfit_species_all_NMDS.csv',
  row.names = TRUE)

# Save top 10 taxa (the ones shown in the plot)
write.csv(mp4_tax_top,
  'data/results/04_taxonomic_profiling/metaphlan4_latest/cfmd/mk_envfit_species_top10_NMDS.csv',
  row.names = TRUE)


## Distance of Countries from Caucasus using Bray–Curtis ------------------

# --- helper: ISO2 -> flag emoji (e.g., 'es' -> '🇪🇸') ---
flag_from_iso2 <- function(iso2) {
  vapply(iso2, function(x) {
    if (is.na(x) || nchar(x) != 2) return('')
    u <- utf8ToInt(toupper(x))
    intToUtf8(u + 127397)
  }, character(1))
}

# Bray–Curtis distance matrix
mp4_bray_mat <- as.matrix(mp4_bray_dist)

# Align metadata to mp4_bray_mat row order (critical)
mp4_meta_aligned <- mp4_merged_mk_metadata |>
  filter(sample %in% rownames(mp4_bray_mat)) |>
  distinct(sample, country_code) |>
  slice(match(rownames(mp4_bray_mat), sample))

# Identify Caucasus samples (your rule: end with '_m')
mp4_cauc_samples <- mp4_meta_aligned$sample[grepl('_m$', mp4_meta_aligned$sample)]
mp4_cfmd_samples <- setdiff(mp4_meta_aligned$sample, mp4_cauc_samples)

# Mean Bray–Curtis distance of each cFMD sample to all Caucasus samples
mp4_cfmd_mean_to_cauc <- tibble(sample = mp4_cfmd_samples) |>
  mutate(mean_to_cauc = rowMeans(
    mp4_bray_mat[sample, mp4_cauc_samples, drop = FALSE],
    na.rm = TRUE
  ))

# Summarise per country & top 10 closest (lowest dissimilarity)
mp4_top10_countries_bray <- mp4_cfmd_mean_to_cauc |>
  left_join(mp4_meta_aligned, by = 'sample') |>
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
    flag = flag_from_iso2(iso2),
    axis_label = paste0(country_name, '  ', flag)
  )

# Named vector for axis labels
mp4_axis_labs <- setNames(mp4_top10_countries_bray$axis_label, mp4_top10_countries_bray$country_name)

# Ranked bar plot (same style as your NMDS version)
mp4_bray_plot_distance <- ggplot(
  mp4_top10_countries_bray,
  aes(
    x = reorder(country_name, -mean_distance),
    y = mean_distance
  )
) +
  geom_col(width = 0.7, fill = '#624E88') +
  geom_text(aes(label = round(mean_distance, 3)),
            hjust = -0.1, size = 4) +
  scale_x_discrete(labels = mp4_axis_labs) +
  coord_flip() +
  labs(
    title = '(C)',
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
  expand_limits(y = max(mp4_top10_countries_bray$mean_distance, na.rm = TRUE) * 1.1)

# Check the plot
mp4_bray_plot_distance

# Save the plot
ggsave(
  mp4_bray_plot_distance,
  filename = 'plots/04_taxonomic_profiling/metaphlan4_latest/cfmd/mk_bray_country_distance.tiff',
  height = 15,
  width = 25,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine Beta Diversity and Distance -------------------------------------

mp4_cfmd_combined_plot <- (mp4_nmds_plot_ctry_taxa | mp4_bray_plot_distance) +
  plot_layout(widths = c(1, 0.7))

# Check the plot
mp4_cfmd_combined_plot

# Save the plot
ggsave(
  mp4_cfmd_combined_plot,
  filename = 'plots/04_taxonomic_profiling/metaphlan4_latest/cfmd/cfmd_combined_plot_final.tiff',
  height = 18,
  width = 44,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Dispersion Analysis -----------------------------------------------------

# Taxonomy dispersion by dataset (Caucasus vs cFMD)
mp4_bd_dataset <- betadisper(mp4_bray_dist, group = mp4_nmds_df_ctry$shape)
permutest(mp4_bd_dataset, permutations = 999)

# Plot distances to centroid (within-group dispersion)
# NOTE: betadisper runs a PCoA internally (on distances-to-centroid space),
# so the axes *do* have eigenvalues and % variance explained.

# Extract site scores (PCoA coordinates) and centroids
mp4_sites <- as.data.frame(scores(mp4_bd_dataset, display = 'sites'))
mp4_sites$sample <- rownames(mp4_sites)
mp4_sites$dataset <- mp4_bd_dataset$group

mp4_centroids <- as.data.frame(scores(mp4_bd_dataset, display = 'centroids'))
mp4_centroids$dataset <- rownames(mp4_centroids)

# Add centroid coordinates per sample (for spider segments)
mp4_sites <- mp4_sites |>
  left_join(
    mp4_centroids |> rename(centroid_x = PCoA1, centroid_y = PCoA2),
    by = 'dataset'
  )

# % variance explained by PCoA 1 and 2 (MetaPhlAn / taxonomy)
# (use only positive eigenvalues for the denominator, consistent with PCoA reporting)
mp4_eig <- mp4_bd_dataset$eig
mp4_eig_pos <- mp4_eig[mp4_eig > 0]

mp4_pcoa1_pct <- round(100 * mp4_eig_pos[1] / sum(mp4_eig_pos), 1)
mp4_pcoa2_pct <- round(100 * mp4_eig_pos[2] / sum(mp4_eig_pos), 1)

# Convex hull coordinates per dataset (to mimic base plot outline)
mp4_hulls <- mp4_sites |>
  group_by(dataset) |>
  slice(chull(PCoA1, PCoA2)) |>
  ungroup()

# Colors (as requested)
mp4_cols <- c('Caucasus' = '#E8403D', 'cFMD' = '#4FACAC')

# Create label text based on betadisper permutest result
mp4_bd_label_title <- data.frame(
  PCoA1 = -0.45,
  PCoA2 = 0.65,
  label = "bold('Dispersion Test')"
)

mp4_bd_label_stats <- data.frame(
  PCoA1 = -0.45,
  PCoA2 = 0.61,
  label = "italic(F)*' = 6.7, '*italic(p)*' = 0.016'"
)

# Build ggplot object
mp4_bd_plot <- ggplot(mp4_sites, aes(x = PCoA1, y = PCoA2)) +
  # Spider segments: sample -> centroid
  geom_segment(
    aes(x = centroid_x, y = centroid_y, xend = PCoA1, yend = PCoA2, color = dataset),
    alpha = 0.35,
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  # Hull outline per dataset
  geom_polygon(
    data = mp4_hulls,
    aes(x = PCoA1, y = PCoA2, group = dataset, color = dataset, fill = dataset),
    alpha = 0.08,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  # Points
  geom_point(aes(color = dataset, shape = dataset), size = 2.5, alpha = 0.9) +
  # Centroids
  geom_point(
    data = mp4_centroids,
    aes(x = PCoA1, y = PCoA2, color = dataset),
    size = 4,
    shape = 4,   # X marker for centroid
    stroke = 1.1,
    show.legend = FALSE
  ) +
  geom_text(
    data = mp4_bd_label_title,
    aes(x = PCoA1, y = PCoA2, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 4.5
  ) +
  geom_text(
    data = mp4_bd_label_stats,
    aes(x = PCoA1, y = PCoA2, label = label),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 4.5
  ) +
  # Label only Caucasus milk kefir samples
  geom_text_repel(
    data = subset(mp4_sites, dataset == 'Caucasus'),
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
  scale_color_manual(values = mp4_cols, name = 'Datasets') +
  scale_fill_manual(values = mp4_cols) +
  scale_shape_manual(values = c('Caucasus' = 17, 'cFMD' = 16), name = 'Datasets') +
  labs(
    title = '(A)',
    x = paste0('PCoA 1 (', mp4_pcoa1_pct, '%)'),
    y = paste0('PCoA 2 (', mp4_pcoa2_pct, '%)')
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
mp4_bd_plot

# Save the plot
ggsave(
  mp4_bd_plot,
  filename = 'plots/04_taxonomic_profiling/metaphlan4_latest/cfmd/mp4_betadispers_plot_final.tiff',
  height = 15,
  width = 22,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)