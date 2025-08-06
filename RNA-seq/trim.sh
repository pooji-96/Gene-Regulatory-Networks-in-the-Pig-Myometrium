#!/bin/bash
#SBATCH -J trim.sh
#SBATCH -p general        # Partition name (check with your HPC)
#SBATCH -o trim.txt
#SBATCH -e trim.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com     # Your email for notifications
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH -A <project_number>     # Your HPC project or account ID

# Set the working directory to your project directory
cd /path/to/project/root

# Create a directory to store trimmed files
mkdir trim

# Environment Setup
module load python/3.11.4
pip install cutadapt     # Cutadapt needs python so this step is done outside the virtual env RNA
module load trimgalore/0.6.10

# List of SRA accessions corresponding to RNA-seq samples
accession=("SRR21949253" "SRR21949254" "SRR21949255" "SRR21949256" "SRR21949257" "SRR21949258" "SRR21949259" "SRR21949260" "SRR21949261" "SRR21949262" "SRR21949263" "SRR21949264" "SRR21949265" "SRR21949266" "SRR21949267" "SRR21949268" "SRR21949269" "SRR21949270")

# Trimming
for acc in "${accession[@]}"; do
    trim_galore --illumina --phred33 --paired --fastqc "${acc}_1.fastq" "${acc}_2.fastq" -o trim
done

