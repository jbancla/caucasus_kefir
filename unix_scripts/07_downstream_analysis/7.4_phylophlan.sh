############################################################################################################################################################
PhyloPhlAn script for reconstructing the tree of life in prokaryotes
############################################################################################################################################################
cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/07_downstream_analysis/7.4_phylophlan/hq_mags
nano 7.4_phylophlan_hq_mags.sh
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=phylophlan_HQ
#SBATCH --error=7.4_phylophlan_hq_mags.err
#SBATCH --output=7.4_phylophlan_hq_mags.out
#SBATCH -p GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=30
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Load phylophlan tool
module load phylophlan/3.1.1

    phylophlan \
        -i input_genomes \
        -d phylophlan \
        --diversity high \
        --accurate \
        -f tol.cfg \
        -o output_tol \
        --genome_extension .fa \
        --verbose \
        --nproc "$SLURM_CPUS_PER_TASK"

# Unload phylophlan tool
module unload phylophlan/3.1.1

############################################################################################################################################################