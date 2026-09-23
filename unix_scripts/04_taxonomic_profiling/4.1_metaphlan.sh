############################################################################################################################################################
MetaPhlAn script for taxonomic profiling
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=metaphlan4
#SBATCH --error=4.1_metaphlan.err
#SBATCH --output=4.1_metaphlan.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=30
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input, output, interleaved, and temporary directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/1.2_fasnicar"
temp_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/temp_dir"
interleaved_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/04_taxonomic_profiling/4.1_metaphlan/4.1.3_interleaved_files"
metaphlan_outdir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/04_taxonomic_profiling/4.1_metaphlan"
metaphlan_dbdir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/databases/metaphlan_db_22102024"

# Create necessary directories.
mkdir -p "$temp_dir"
mkdir -p "$interleaved_dir"

# Find all directory of interest from the in_dir.
for folder in $(find "$in_dir" -type d -maxdepth 1 | sed '1d' | sort -u); do

# Print only the base name from the in_dir.
folder_base=$(echo "$folder" | awk -F / '{print $NF}')

    for filename in $(find "$folder" -type f -name "*R?.fastq.bz2" | sed 's/_R..fastq.*//g' | sort -u | awk -F / '{print $NF}'); do

        find "$folder" -type f \( -name "*${filename}*" -a -name "*_R*.fastq.bz2" \) -exec cp {} "$temp_dir" \;

        echo "Processing $folder."

        # Check for .bz2 files.
        find "$temp_dir" -type f -name "*.bz2" -exec bash -c '
            for file; do
                bunzip2 "$file"
            done
        ' bash {} +

        # Check for .gz files.
        find "$temp_dir" -type f -name "*.gz" -exec bash -c '
            for file; do
                gunzip "$file"
            done
        ' bash {} +

        echo "Uncompressed files currently located at $temp_dir."

        files=("$temp_dir"/*.fastq*)

        if [ ${#files[@]} -eq 2 ]; then
            # Save the file paths as separate variables.
            file1="${files[0]}"
            file2="${files[1]}"
        else
            echo "Error: There are not exactly two FASTQ files in $temp_dir."
        fi

        echo "Now processing $file1 and $file2."

        # Concatenate R1 and R2 for MetaPhlan
        cat "$file1" "$file2" > "$interleaved_dir"/${filename}_interleaved.fastq

	mkdir -p "$metaphlan_outdir/4.1.1_results"
	mkdir -p "$metaphlan_outdir/4.1.2_bowtie"

	results_dir="$metaphlan_outdir/4.1.1_results"
	bowtie_dir="$metaphlan_outdir/4.1.2_bowtie"

        # Load the metaphlan tool.
        module load metaphlan/4.1.1

        # Run the metaphlan script.
        metaphlan "$interleaved_dir"/${filename}_interleaved.fastq --input_type fastq \
            --output_file "$results_dir"/${filename}_profile.txt \
	        --bowtie2out "$bowtie_dir"/${filename}.bowtie2out.txt \
	        --bowtie2db "$metaphlan_dbdir" --index mpa_vOct22_CHOCOPhlAnSGB_202403 \
	        --add_viruses --unclassified_estimation --nproc "$SLURM_CPUS_PER_TASK"

        # Unload the metaphlan tool.
        module unload metaphlan/4.1.1

        # Clean up temporary files.
        rm "$temp_dir"/*

    done # End of the filename loop.

done # End of the folder loop.

rm -r "$temp_dir"

############################################################################################################################################################
- END -
