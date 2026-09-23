############################################################################################################################################################
fastqc script for quality control check
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=fastqc
#SBATCH --error=1.1_fastqc.err
#SBATCH --output=1.1_fastqc.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=25
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input and output directory.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/00_download"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/1.1_fastqc"

# Load the quality control analysis tool.
module load fastqc/0.12.1

# Find all directory of interest from the in_dir.
for folder in $(find "$in_dir" -type d | sed '1d' | sort -u); do

# Print only the base name from the in_dir.
folder_base=$(echo "$folder" | awk -F / '{print $NF}')

    # Find all the files *.fastq.gz from the folder.
    for filename in $(find "$folder" -type f -name "*.fastq.gz" | awk -F / '{print $NF}' | sed 's/_R._001.fastq.gz//g' | sort -u); do

    # Make an output directory that contains the modified file name.
    mkdir -p "$out_dir/$folder_base/$filename"

	# Find all .fastq.gz files in the output directory.
	for fastqfile in "$in_dir/$folder_base"/${filename}_R?_001.fastq.gz; do

	echo "File processing $fastqfile"

	# Run fastqc script.
	fastqc "$fastqfile" --threads "$SLURM_CPUS_PER_TASK" -o "$out_dir/$folder_base/$filename"

	done # End of fastqfile loop.

    done # End of filename loop.

done # End of folder loop.

# Unload the quality control analysis tool.
module unload fastqc/0.12.1

#######################################################################################################################################################################
- END -
