############################################################################################################################################################
DRAM script for the annotation and curation of function for microbial and viral genomes
############################################################################################################################################################
cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/scripts/07_downstream_analysis
nano 7.2_dram_hq_mags.sh 
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=dram_HQ
#SBATCH --error=7.2_dram_hq_mags.err
#SBATCH --output=7.2_dram_hq_mags.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input and output directories.
bin_names="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/mags/hq_mags_filenames.txt"
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/mags/hq_mags"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/07_downstream_analysis/7.2_dram"

# Create necessary directories.
mkdir -p "$out_dir"/hq_mags

annotation_outdir="$out_dir/hq_mags/annotation"
distilled_outdir="$out_dir/hq_mags/distilled"

mkdir -p "$annotation_outdir"
mkdir -p "$distilled_outdir"

# Load DRAM tool.
module load dram/1.5.0
source activate DRAM

for x in $(cat "$bin_names"); do

    echo "Started annotating bin $x"

    # Annotation

    DRAM.py annotate -i "$in_dir"/"$x".fa -o "$annotation_outdir/$x" --threads "$SLURM_CPUS_PER_TASK" 

    echo "Done annotating bin $x"

    echo "Started distilling bin $x"

    # Distill

    DRAM.py distill -i "$annotation_outdir/$x"/annotations.tsv -o "$distilled_outdir/$x" \
        --trna_path "$annotation_outdir/$x"/trnas.tsv --rrna_path "$annotation_outdir/$x"/rrnas.tsv

    echo "Done distilling bin $x"

done # End of x loop.

# Unload DRAM tool.
module unload dram/1.5.0

############################################################################################################################################################