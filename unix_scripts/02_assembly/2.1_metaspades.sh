############################################################################################################################################################
MetaSPAdes script for metagenomics assembly of clean reads
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=metaspades
#SBATCH --error=2.1_metaspades.err
#SBATCH --output=2.1_metaspades.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=25
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input and output directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/clean_reads"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/02_assembly/2.1_metaspades"

# Create necessary directories.
mkdir -p "$out_dir"

# Load the assembly tool.
module load spades/3.15.5

# Find all directory of interest from the in_dir.
for folder in $(find "$in_dir" -type d | sed '1d' | sort -u); do

    # Find all the interleaved files from the folder.
    for filename in $(find "$folder" -type f -name "*fastq.gz" | awk -F / '{print $NF}' | sed 's/_R..fastq.gz//g' | sort -u); do

    # Make an output directory that contains the folder_base and the filename.
     mkdir -p "$out_dir/$filename"

    # Find all fastq.gz files in the input directory.
    file1="$folder/${filename}_R1.fastq.gz"
    file2="$folder/${filename}_R2.fastq.gz"

    # # Check if fastq_file exist before proceeding to metaspades assembly.
    if [ -f "$file1" ] && [ -f "$file2" ]; then
	# If both file1 and file2 are found from the folder, perform metaspades assembly.
	echo "File $file1 and $file2 found. Performing metaspades assembly..."

        # Run metaspades script.
        metaspades.py -1 "$file1" -2 "$file2" -t "$SLURM_CPUS_PER_TASK" -o "$out_dir/$filename"

    else
	# If either file is missing, skip this sample.
	echo "File $file1 and/or $file2 NOT found. Skipping..."
    fi

	# Rename the output files to include the sample name.
	mv "$out_dir/$filename"/contigs.fasta "$out_dir/$filename"/${filename}_contigs.fasta
	mv "$out_dir/$filename"/scaffolds.fasta "$out_dir/$filename"/${filename}_scaffolds.fasta

	echo "Generated $out_dir/$filename/${filename}_contigs.fasta successfully!"

    done # End of filename loop.

done # End of folder loop.

# Unload the assembly tool.
module unload spades/3.15.5

############################################################################################################################################################
- END -
