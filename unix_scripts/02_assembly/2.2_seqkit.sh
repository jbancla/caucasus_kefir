######################################################################################################################################################################
SeqKit script for filtering assembled contigs longer than 1000bp
######################################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=seqkit
#SBATCH --error=2.2_seqkit.err
#SBATCH --output=2.2_seqkit.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input and output directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/02_assembly/2.1_metaspades"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/02_assembly/2.1_metaspades"

# Load the seqkit tool.
module load seqkit/2.8.0

# Find all directory of interest from the in_dir.
for folder in $(find "$in_dir" -type d -maxdepth 1 | sed '1d' | sort -u); do

# Print only the base name from the in_dir.
folder_base=$(echo "$folder" | awk -F / '{print $NF}')

    # Find all the sample names from the folder.
    for filename in $(find "$folder" -type f ! -name "first_pe_contigs.fasta" -name "*contigs.fasta" -maxdepth 1 | awk -F / '{print $NF}' | sed 's/_contigs.fasta//g' | sort -u); do

    # Find all contigs.fasta file from the input directory.
    assembled_contigs="$folder/${filename}_contigs.fasta"

    # Check if assembled_contigs exist before proceeding to seqkit filtering.
    if [ -f "$assembled_contigs" ]; then
	# If assembled_contigs is found from the folder, perform seqkit filtering.
	echo "File $assembled_contigs found. Performing seqkit filtering..."

        # Run seqkit script.
        seqkit seq -m 1000 "$assembled_contigs" > "$out_dir/$folder_base"/${filename}_contigs_filtered.fasta
	    seqkit stats -T "$out_dir/$folder_base"/${filename}_contigs_filtered.fasta > "$out_dir/$folder_base"/${filename}_contigs_filtered_stats.tsv

    else
	# If assembled_contigs is missing, skip this sample.
	echo "File $assembled_contigs NOT found. Skipping..."
    fi

	echo "Generated $out_dir/$folder_base/${filename}_contigs_filtered.fasta successfully!"

    done # End of filename loop.

done # End of folder loop.

# Unload the seqkit tool.
module unload seqkit/2.8.0

######################################################################################################################################################################
- END -
