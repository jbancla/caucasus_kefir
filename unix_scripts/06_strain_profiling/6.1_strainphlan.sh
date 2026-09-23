############################################################################################################################################################
StrainPhlAn script for performing strain-level population genomics on large metagenomics datasets
############################################################################################################################################################
# STEP 1
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=strainphlan_step1
#SBATCH --error=6.1_strainphlan_s1.err
#SBATCH --output=6.1_strainphlan_s1.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=35
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# STEP 1: Run metaphlan to obtain sam output files.

# Define input, output, and metaphlan database directories.
strain_names="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/strain_filenames.txt"
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/01_preprocessing/fastq_dir"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan"
metaphlan_dbdir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/databases/metaphlan_db"

# Create necessary directories.
mkdir -p "$out_dir"/6.1.1_interleaved
mkdir -p "$out_dir"/6.1.1_sams
mkdir -p "$out_dir"/6.1.1_bowtie2
mkdir -p "$out_dir"/6.1.1_profiles

for x in $(cat "$strain_names"); do

    # Create first an interleaved file using the reformat.sh script.

    # Load bbmap tool.
    module load bbmap/39.06

    echo "Started running reformat.sh on sample $x."

    # Run reformat.sh script.
    reformat.sh in="$in_dir/$x"/${x}_R1.fastq \
        in2="$in_dir/$x"/${x}_R2.fastq \
        out="$out_dir"/6.1.1_interleaved/${x}_interleaved.fastq

    echo "Done running reformat.sh on sample $x."

    # Use the interleaved file as input to run metaphlan.

    # Load metaphlan tool.
    module load metaphlan/4.1.1

    echo "Started running metaphlan on sample $x."

    # Run metaphlan script.
    metaphlan "$out_dir"/6.1.1_interleaved/${x}_interleaved.fastq \
        --input_type fastq \
        -s "$out_dir"/6.1.1_sams/${x}.sam.bz2 \
        --bowtie2db "$metaphlan_dbdir" \
        --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
        --bowtie2out "$out_dir"/6.1.1_bowtie2/${x}.bowtie2.bz2 \
        -o "$out_dir"/6.1.1_profiles/${x}_profile.tsv \
        --nproc "$SLURM_CPUS_PER_TASK"

    echo "Done running metaphlan on sample $x."

done # End of x loop.

# Unload bbmap and metaphlan tool.
module unload bbmap/39.06
module unload metaphlan/4.1.1

############################################################################################################################################################
# STEP 2
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=strainphlan_step2
#SBATCH --error=6.1_strainphlan_s2.err
#SBATCH --output=6.1_strainphlan_s2.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=35
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# STEP 2: Run samples-to-markers to generate a marker file for each sample.

# Strainphlan will be looking for the database you use in STEP 1, so set the latest mpa database in HPC. Use this:
# module load metaphlan/4.1.1
# export METAPHLAN_DB_DIR='/data/Food/analysis/R1838_DOMINO/joseph_ancla/databases/metaphlan_db'

# Define input and output directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan/6.1.1_sams"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan"

# Create necessary directories.
mkdir -p "$out_dir"/6.1.2_consensus_markers

# Load metaphlan and strainphlan tool.
module load metaphlan/4.1.1
module load strainphlan/4.0.6

    echo "Started running sample2markers on sam files."

    # Run sample2markers.py script.
    sample2markers.py -i "$in_dir"/*.sam.bz2 \
        -o "$out_dir"/6.1.2_consensus_markers \
        --nproc "$SLURM_CPUS_PER_TASK"

    echo "Done running sample2markers on sam files."

# Unload metaphlan and strainphlan tool.
module unload metaphlan/4.1.1
module unload strainphlan/4.0.6

############################################################################################################################################################
# STEP 3
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=strainphlan_step3
#SBATCH --error=6.1_strainphlan_s3.err
#SBATCH --output=6.1_strainphlan_s3.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=35
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# STEP 3: Extract clade markers from MetaPhlAn database.
# Run the script inside the 6.1.3_clade_markers directory.

# Define input, output, and metaphlan database directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan/6.1.2_consensus_markers"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan"
mpadb_pkl="/data/Food/analysis/R1838_DOMINO/joseph_ancla/databases/metaphlan_db"

# Create necessary directories.
mkdir -p "$out_dir"/6.1.3_clade_markers

cd "$out_dir"/6.1.3_clade_markers

# Load strainphlan tool.
module load strainphlan/4.0.6

    echo "Started extracting clade markers on json files."

    # Run strainphlan script.
    strainphlan -s "$in_dir"/*.json.bz2 \
        -d "$mpadb_pkl"/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl \
        --print_clades_only \
        -o . \
        --nproc "$SLURM_CPUS_PER_TASK" > clades.txt

    echo "Done extracting clade markers on json files."

# Unload strainphlan tool.
module unload strainphlan/4.0.6

# This will produce "print_clades_only.tsv" that contains your samples clades.
# I don't know what happen to clades.txt but it contains the log of what happened to the script instead of the output file. 
# Well as long as it's working then it's good to go!

############################################################################################################################################################
# STEP 3b
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=strainphlan_step3b
#SBATCH --error=6.1_strainphlan_s3b.err
#SBATCH --output=6.1_strainphlan_s3b.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=35
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# Continuation of STEP 3, make sure you are in 6.1.3_clade_markers directory, then run:
# cut -f1 print_clades_only.tsv | tail -n +2 | sort -u > clade_list.txt

# Define input, output, and metaphlan database directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan/6.1.3_clade_markers"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan"
mpadb_pkl="/data/Food/analysis/R1838_DOMINO/joseph_ancla/databases/metaphlan_db"

# Load strainphlan tool.
module load strainphlan/4.0.6

for x in $(cat "$in_dir"/clade_list.txt); do

    # Run extract_markers.py script.
    extract_markers.py -c "$x" \
        -d "$mpadb_pkl"/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl \
        -o "$out_dir"/6.1.3_clade_markers

done # End of x loop.

# Unload strainphlan tool.
module unload strainphlan/4.0.6

############################################################################################################################################################
# STEP 4
############################################################################################################################################################

#!/bin/sh
#SBATCH --job-name=strainphlan_step4
#SBATCH --error=6.1_strainphlan_s4.err
#SBATCH --output=6.1_strainphlan_s4.out
#SBATCH -p Priority,Background,GPU
#SBATCH -n 1
#SBATCH --cpus-per-task=35
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Joseph.Ancla@teagasc.ie

# STEP 4: Build the multiple sequence alignment and the phylogenetic tree.

# Define input, output, clade, and metaphlan database directories.
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan/6.1.2_consensus_markers"
clade_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan/6.1.3_clade_markers"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan"
mpadb_pkl="/data/Food/analysis/R1838_DOMINO/joseph_ancla/databases/metaphlan_db"

# Create necessary directories.
mkdir -p "$out_dir"/6.1.4_output

# Load strainphlan tool.
module load strainphlan/4.0.6

for x in $(cat "$clade_dir"/clade_list.txt); do

    echo "Started building multiple sequence alignment on json files."

    mkdir -p "$out_dir"/6.1.4_output/"$x"

    # Run strainphlan script.
    strainphlan -s "$in_dir"/*.json.bz2 \
        -m "$clade_dir"/${x}.fna \
        -o "$out_dir"/6.1.4_output/"$x" \
        -c "$x" \
        -d "$mpadb_pkl"/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl \
        --mutation_rates \
        --phylophlan_mode accurate \
        --nproc "$SLURM_CPUS_PER_TASK"

    echo "Done building multiple sequence alignment on json files."

done # End of x loop.

# Unload strainphlan tool.
module unload strainphlan/4.0.6

############################################################################################################################################################
Extraction of pairwise phylogenetic distances
############################################################################################################################################################

# Define input and output directories.
clade_names="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan/6.1.3_clade_markers"
in_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan/6.1.4_output"
out_dir="/data/Food/analysis/R1838_DOMINO/joseph_ancla/domino_files/mk_caucasus/06_strain_profiling/6.1_strainphlan"

# Create necessary directories.
mkdir -p "$out_dir"/6.1.5_phylogenetic_distances

# Load python tool.
conda activate python

for x in $(cat "$clade_names"/clade_list.txt); do

    echo "Started calculating phylogenetic distances on clade $x."

    # Run tree_pairwisedists.py script.
    python tree_pairwisedists.py -n "$in_dir/$x"/RAxML_bestTree.${x}.StrainPhlAn4.tre "$out_dir"/6.1.5_phylogenetic_distances/${x}_phylogenetic_distance.tsv

    echo "Done calculating phylogenetic distances on clade $x."

done # End of x loop.
    
# Unload python tool.
conda deactivate

############################################################################################################################################################
- END -
