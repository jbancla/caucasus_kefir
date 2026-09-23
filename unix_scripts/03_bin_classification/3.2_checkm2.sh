############################################################################################################################################################
CheckM2 script for assessing the quality of reassembled bins
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=checkm2_HQ
#SBATCH --error=3.2_checkm2_HQ_MAGS.err
#SBATCH --output=3.2_checkm2_HQ_MAGS.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input, output, and checkm2 database directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/hq_mags"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/3.2_checkm2"
checkm2_dbdir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/databases/checkm2"

# Create necessary directories.
mkdir -p "$out_dir"/hq_mags

# Load the checkm2 tool.
module load anaconda/3.10
source activate base
conda activate checkm2

    echo "Started checkm2 analysis."

    # Run checkm2 script.
    checkm2 predict --threads "$SLURM_CPUS_PER_TASK" -x fa \
        --input "$in_dir" --output-directory "${out_dir}/hq_mags" \
        --database_path "$checkm2_dbdir"/uniref100.KO.1.dmnd

    echo "Done with checkm2 analysis."

# Unload the checkm2 tool.
conda deactivate
module unload anaconda/3.10

############################################################################################################################################################
- END -
