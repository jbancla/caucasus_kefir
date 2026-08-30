############################################################################################################################################################
Fasnicar script for filtering and trimming raw reads
############################################################################################################################################################
cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/scripts/01_preprocessing
nano 1.2_fasnicar.sh
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=fasnicar
#SBATCH --error=1.2_fasnicar.err
#SBATCH --output=1.2_fasnicar.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input, output, and index directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/00_download"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/1.2_fasnicar"
index_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/databases/built_index/phix_cow_human"

# Load the filtering-trimming tool.
module load fasnicar/0.2.4

# Find all directory of interest from the in_dir.
for folder in $(find "$in_dir" -type d | sed '1d' | sort -u | sed -n '1,10p'); do

# Print only the base name from the in_dir.
folder_base=$(echo "$folder" | awk -F / '{print $NF}')

    # Find all the files *.fastq.gz from the in_dir.
    for filename in $(find "$folder" -type f -name "*.fastq.gz" | awk -F / '{print $NF}' | sed 's/_R._001.fastq.gz//g' | sort -u); do

    # Make an output directory that contains the filename.
    mkdir -p "$out_dir/$folder_base"

    # Copy the modified directory to the output directory.
    cp "$folder"/${filename}_* "$out_dir/$folder_base"

        # Run fasnicar script.
        preprocess.new.py -e .fastq.gz -f _R1 -r _R2 -x "$index_dir" --verbose -n "$SLURM_CPUS_PER_TASK" -i "$out_dir/$folder_base"

    done # End of filename loop.

done # End of folder loop.

# Unload the filtering-trimming tool.
module load fasnicar/0.2.4

#######################################################################################################################################################################