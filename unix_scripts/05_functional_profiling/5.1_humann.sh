############################################################################################################################################################
Humann (HMP Unified Metabolic Analysis Network) script for functional profiling of metagenomic data
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
- END -
