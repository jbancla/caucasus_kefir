############################################################################################################################################################
Humann (HMP Unified Metabolic Analysis Network) script for functional profiling of metagenomic data
############################################################################################################################################################
cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/scripts/05_functional_profiling
nano 5.1_humann.sh
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=humann3
#SBATCH --error=5.1_humann.err
#SBATCH --output=5.1_humann.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=25
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input, output, combined, humann results, taxonomic profile, nucleotide, and protein database directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/1.2_fasnicar"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/05_functional_profiling/5.1_humann_caucasus"
combined_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/05_functional_profiling/5.1_humann/5.1.1_combined_files"
humann_results="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/05_functional_profiling/5.1_humann_caucasus/5.1.2_humann_results"
tax_profile="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/04_taxonomic_profiling/4.1_metaphlan_caucasus_bowtie"
nucleotide_db="/data/databases_food/humann_db/chocophlan"
protein_db="/data/databases_food/humann_db/uniref"

# Create necessary directories.
mkdir -p "$out_dir" "$humann_results"

# Load HUMAnN3 tool.
module load humann/3.8

# Find all directory of interest from the in_dir.
for folder in $(find "$in_dir" -maxdepth 1 -type d | sed '1d' | sort -u); do

    # Find all the files *.fastq.gz from the folder.
    for filename in $(find "$folder" -type f -name "*_R?.fastq.bz2" | awk -F / '{print $NF}' | sed 's/_R..fastq.bz2//g' | sort -u); do

    # Make an output directory that contains the folder_base and the filename.
    mkdir -p "$humann_results/$filename"

	# Run humann script.
	humann --input "$combined_dir"/${filename}.fastq \
	    --output "$humann_results/$filename" \
	    --nucleotide-database "$nucleotide_db" \
	    --protein-database "$protein_db" \
	    --taxonomic-profile "$tax_profile"/${filename}_profile.txt \
	    --threads "$SLURM_CPUS_PER_TASK"

	echo "Processed sample $filename successfully!"

    done # End of filename loop.

done # End of folder loop.

# Unload HUMAnN3 tool.
module unload humann/3.8

############################################################################################################################################################
Humann script for normalizing pathabundance.tsv (Counts per Million)
############################################################################################################################################################
cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/scripts/05_functional_profiling
nano 5.2_humann_caucasus_pa_normalized_cpm.sh
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=hm3_cpm_caucasus
#SBATCH --error=5.2_humann_caucasus_pa_normalized_cpm.err
#SBATCH --output=5.2_humann_caucasus_pa_normalized_cpm.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input and output directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/05_functional_profiling/5.1_humann_caucasus/5.1.2_humann_results"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/05_functional_profiling/5.1_humann_caucasus/5.1.2_humann_results"

# Load the functional profiling tool.
module load humann/3.8

# Find all directory of interest from the in_dir.
for folder in $(find "$in_dir" -maxdepth 1 -type d | sed '1d' | sort -u); do

# Print only the base name from the in_dir.
folder_base=$(echo "$folder" | awk -F / '{print $NF}')

    # Find all the files *.fastq.gz from the folder.
    for filename in $(find "$folder" -type f -name "*_pathabundance.tsv" | awk -F / '{print $NF}' | sed 's/_pathabundance.tsv//g' | sort -u); do

    # Find all pathabundance.tsv file from the input directory.
    tsv_file="$folder/${filename}_pathabundance.tsv"

    # Check if tsv_file exist before proceeding to humann normalization cpm.
    if [ -f "$tsv_file" ]; then
	# If tsv_file is found from the folder, perform humann normalization cpm.
	echo "File $tsv_file found. Performing humann normalization cpm..."

        # Run humann script.
        humann_renorm_table --input "$tsv_file" \
            --output "$out_dir/$folder_base"/${filename}_pathabundance_cpm.tsv \
            --units cpm --update-snames

    else
	# If tsv_file is missing, skip this sample.
	echo "File $tsv_file NOT found. Skipping..."
    fi

	echo "Generated $out_dir/$folder_base/${filename}_pathabundance_cpm.tsv successfully!"

    done # End of filename loop.

done # End of folder loop.

# Unload the functional profiling tool.
module unload humann/3.8

############################################################################################################################################################
Script for joining all pathabundance.tsv files
############################################################################################################################################################

# Define input and output directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/05_functional_profiling/5.1_humann_caucasus/5.1.2_humann_results"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/05_functional_profiling/5.1_humann_caucasus/5.1.3_humann_joined"
pathabundance_outdir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/05_functional_profiling/5.1_humann_caucasus/5.1.3_humann_joined/pathabundance_cpm"

mkdir -p "$pathabundance_outdir"

# Copy all *_cpm.tsv files into one directory
find "$in_dir" -type f -name '*_pathabundance_cpm.tsv' -exec cp {} "$pathabundance_outdir" \;

# Load the functional profiling tool.
module load humann/3.8

    humann_join_tables \
        --input "$pathabundance_outdir" \
        --file_name pathabundance_cpm \
        --output "$out_dir"/hm3_caucasus_joined_pathabundance_cpm.tsv

    # Run SEPARATELY stratify_table script (go first to the input directory).
    humann_split_stratified_table --input "$out_dir"/hm3_caucasus_joined_pathabundance_cpm.tsv --output "$out_dir"/caucasus_stratified_pathabundance_cpm

# Unload the functional profiling tool.
module unload humann/3.8

############################################################################################################################################################