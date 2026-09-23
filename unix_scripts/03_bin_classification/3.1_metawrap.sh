############################################################################################################################################################
MetaWRAP script for binning, bin refinement, and reassembly of assembled contigs
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=metawrap_9005
#SBATCH --error=3.1_metawrap_metaspades_9005.err
#SBATCH --output=3.1_metawrap_metaspades_9005.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=25
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Define samples, input, and output directory.
samples="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/filenames.txt"
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/02_assembly/2.1_metaspades"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification/3.1_metawrap_metaspades_9005"

# Load the metawrap tool.
module load metawrap/1.3.2env

cd /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/03_bin_classification
mkdir 3.1_metawrap_metaspades_9005

for x in $(cat "$samples"); do

    cp /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/fastq_dir_v2/"$x"/${x}_1.fastq \
    /data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/fastq_dir_v2/"$x"/${x}_2.fastq .
    mkdir 3.1_metawrap_metaspades_9005/"$x"

# STEP 1: Initial Binning

    echo "Started binning on ${x}"

	# Run metawrap initial binning script.
	metawrap binning -o ./3.1_metawrap_metaspades_9005/"$x"/INITIAL_BINNING -t 25 -m 200 \
	    -a "$in_dir"/"$x"/${x}_contigs_filtered.fasta --metabat2 --maxbin2 --concoct ./*fastq

    echo "Done binning on ${x}"

# STEP 2: Bin Refinement

    echo "Started bin refinement on ${x}"

	# Run metawrap bin refinement script.
	metawrap bin_refinement -o ./3.1_metawrap_metaspades_9005/"$x"/BIN_REFINEMENT --quick -m 200 -t 25 \
	    -A ./3.1_metawrap_metaspades_9005/"$x"/INITIAL_BINNING/metabat2_bins \
	    -B ./3.1_metawrap_metaspades_9005/"$x"/INITIAL_BINNING/maxbin2_bins \
	    -C ./3.1_metawrap_metaspades_9005/"$x"/INITIAL_BINNING/concoct_bins -c 90 -x 5

    echo "Done bin refinement on ${x}"

# STEP 3: Reassembly

    echo "Started bin reassembly on ${x}"

	# Run metawrap bin reassembly script.
	metawrap reassemble_bins -t 25 -m 200 -o ./3.1_metawrap_metaspades_9005/"$x"/BIN_REASSEMBLY \
	    -1 "$x"_1.fastq -2 "$x"_2.fastq -c 90 -x 5 -b ./3.1_metawrap_metaspades_9005/"$x"/BIN_REFINEMENT/metawrap_90_5_bins

    echo "Done bin reassembly on ${x}"

    rm "$x"_1.fastq "$x"_2.fastq

done

# Unload the metawrap tool.
module unload metawrap/1.3.2env

############################################################################################################################################################
- END -
