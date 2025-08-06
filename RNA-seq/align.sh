#!/bin/bash
#SBATCH -J align.sh                                       
#SBATCH -p general        # Partition name (check with your HPC)
#SBATCH -o align.txt
#SBATCH -e align.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com     # Your email for notifications
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --array=0-17      # Job array to run 18 parallel tasks
#SBATCH -A <project_number>     # Your HPC project or account ID

# Set the working directory to your project directory
cd /path/to/project/root

# Environment Setup
module load miniconda
conda activate RNA
module load star/2.7.11a    # Load STAR module

# List of SRA accessions corresponding to RNA-seq samples
accession=("SRR21949253" "SRR21949254" "SRR21949255" "SRR21949256" "SRR21949257" "SRR21949258" "SRR21949259" "SRR21949260" "SRR21949261" "SRR21949262" "SRR21949263" "SRR21949264" "SRR21949265" "SRR21949266" "SRR21949267" "SRR21949268" "SRR21949269" "SRR21949270")

# Get accession ID for the array job
acc=${accession[$SLURM_ARRAY_TASK_ID]}

# Run STAR to align paired-end reads
STAR --runThreadN 10 \
     --genomeDir /path/to/project/root/genome \
     --readFilesIn /path/to/project/root/trim/${acc}_1_val_1.fq /path/to/project/root/trim/${acc}_2_val_2.fq \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix /path/to/project/root/align/${acc}/


