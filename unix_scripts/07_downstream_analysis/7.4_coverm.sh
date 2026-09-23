############################################################################################################################################################
CoverM script for read alignment statistics in metagenomics
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=covermHQbam
#SBATCH --error=7.7_coverm_hq_mags_relab_bam_cache.err
#SBATCH --output=7.7_coverm_hq_mags_relab_bam_cache.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# STEP 2 (rerun): Perform coverage calculation using coverm, this time caching BAM files.

# Define clean reads, MAGs, and output directories.
clean_reads="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/clean_reads"
mag_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/mags/hq_mags"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/07_downstream_analysis/7.7_coverm/hq_mags_relab_bam_cache"

# Create necessary directories.
mkdir -p "$out_dir"/hq_mags_relab
mkdir -p "$out_dir"/bam_cache

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

    # Create a per-sample subdirectory for the cached BAM so file names don't collide across samples.
    sample_bam_dir="${out_dir}/bam_cache/${sample_name}"
    mkdir -p "${sample_bam_dir}"

    # Run coverm script (identical parameters to 7.6, plus BAM caching).
    coverm genome \
    -1 "${r1}" \
    -2 "${r2}" \
    --genome-fasta-directory "${mag_dir}" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --methods mean covered_bases rpkm relative_abundance \
    --genome-fasta-extension fa \
    --bam-file-cache-directory "${sample_bam_dir}" \
    -o "${out_dir}/hq_mags_relab/${sample_name}_mag_coverage_relab.tsv"

done # End of sample_dir loop.

echo "[DONE] Relative abundance outputs in: ${out_dir}/hq_mags_relab"
echo "[DONE] Cached BAM files in: ${out_dir}/bam_cache/<sample_name>/"

# Unload coverm tool.
conda deactivate
module unload anaconda/3.10

############################################################################################################################################################
- END -
