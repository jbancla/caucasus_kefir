############################################################################################################################################################
GTDBTk script for taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy (GTDB)
############################################################################################################################################################
cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/scripts/07_downstream_analysis
nano 7.1_gtdbtk_hq_mags.sh 
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=gtdbtk_HQ
#SBATCH --error=7.1_gtdbtk_hq_mags.err
#SBATCH --output=7.1_gtdbtk_hq_mags.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input, output, and gtdb database directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/hq_mags"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/07_downstream_analysis/7.1_gtdbtk"
gtdb_path="/data/databases_food/GTDB/release226"
mash_db="${gtdb_path}/gtdbtk_ref_sketch.msh"

# Create necessary directories.
mkdir -p "$out_dir"/hq_mags

# Load gtdbtk tool.
module load gtdbtk/2.4.0
source activate gtdbtk-2.4.0
export GTDBTK_DATA_PATH="$gtdb_path"

echo "Started running gtdbtk."

    # Run gtdbtk script.
    gtdbtk classify_wf \
        --genome_dir "$in_dir" \
        --out_dir "$out_dir"/hq_mags \
        --cpus 20 \
        --extension fa \
        --mash_db "$mash_db"

echo "Done running gtdbtk."

conda deactivate
module unload gtdbtk/2.4.0

############################################################################################################################################################