############################################################################################################################################################
CoverM script for read alignment statistics in metagenomics
############################################################################################################################################################
cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/scripts/07_downstream_analysis
nano 7.4_coverm_hq_mags.sh 
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=coverm_HQ
#SBATCH --error=7.4_coverm_hq_mags.err
#SBATCH --output=7.4_coverm_hq_mags.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=25
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define clean reads, MAGs, and output directories.
clean_reads="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/clean_reads"
mag_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/mags/hq_mags"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/07_downstream_analysis/7.4_coverm"

# Create necessary directories.
mkdir -p "$out_dir"/hq_mags

# Load coverm tool.
module load anaconda/3.10
source activate base
conda activate coverm

for sample_dir in "${clean_reads}"/*; do
    [[ -d "${sample_dir}" ]] || continue

    sample_name=$(basename "${sample_dir}")

    r1="${sample_dir}/${sample_name}_R1.fastq.gz"
    r2="${sample_dir}/${sample_name}_R2.fastq.gz"

    # Skip if read files are missing
    if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
        echo "[SKIP] Missing reads for ${sample_name}: ${r1} or ${r2}"
        continue
    fi

    echo "[RUN] ${sample_name}"

    # Run coverm script.
    coverm genome \
	    -1 "${r1}" \
	    -2 "${r2}" \
	    --genome-fasta-directory "${mag_dir}" \
	    --threads "$SLURM_CPUS_PER_TASK" \
	    --min-read-percent-identity 95 \
	    --min-read-aligned-percent 75 \
	    --methods mean covered_bases length rpkm \
	    --genome-fasta-extension fa \
	    -o "${out_dir}/hq_mags/${sample_name}_mag_coverage.tsv"

done # End of sample_dir loop.

echo "[DONE] Outputs in: ${out_dir}/hq_mags"

# Unload coverm tool.
conda deactivate
module unload anaconda/3.10

############################################################################################################################################################