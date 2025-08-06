#!/bin/bash
#SBATCH -J qc.sh
#SBATCH -p general        # Partition name (check with your HPC)
#SBATCH -o qc.txt
#SBATCH -e qc.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com     # Your email for notifications
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH -A <project_number>     # Your HPC project or account ID

# Set the working directory to your project directory
cd /path/to/project/root

# Create a directory to store FASTQC output files
mkdir fastqc

# Environment Setup
module load miniconda
conda activate RNA
module load fastqc/0.12.1     # Load FastQC module

# Run FastQC on all files
fastqc -o fastqc *.fastq 
