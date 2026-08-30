############################################################################################################################################################
Abricate script for mass screening of contigs for antibiotic resistance and virulence genes
############################################################################################################################################################
cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/scripts/07_downstream_analysis
nano 7.3_abricate_hq_mags.sh 
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=abricate_HQ
#SBATCH --error=7.3_abricate_hq_mags.err
#SBATCH --output=7.3_abricate_hq_mags.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define input and output directories.
bin_names="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/mags/hq_mags_filenames.txt"
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/mags/hq_mags"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/07_downstream_analysis/7.3_abricate"

# Create necessary directories.
mkdir -p "$out_dir"/hq_mags

# Load abricate tool.
module load anaconda/3.10
source activate base
conda activate abricate

# Antibiotic Resistance

for x in $(cat "$bin_names"); do

    echo "Started screening for antibiotic resistance genes in $x"

    abricate --db card --quiet "$in_dir"/"$x".fa > "$out_dir"/hq_mags/"$x"_card.tab

    echo "Done screening in $x"

done # End of x loop.

    abricate --summary "$out_dir"/hq_mags/*_card.tab > "$out_dir"/hq_mags/abricate_card_summary.tab

# Virulence Factors

for x in $(cat "$bin_names"); do

    echo "Started screening for virulence factors genes in $x"

    abricate --db vfdb --quiet "$in_dir"/"$x".fa > "$out_dir"/hq_mags/"$x"_vfdb.tab

    echo "Done screening in $x"

done # End of x loop.

    abricate --summary "$out_dir"/hq_mags/*_vfdb.tab > "$out_dir"/hq_mags/abricate_vfdb_summary.tab

# Unload abricate tool.
conda deactivate
module unload anaconda/3.10

############################################################################################################################################################