#!/bin/bash
#SBATCH -J counts.sh
#SBATCH -p general        # Partition name (check with your HPC)
#SBATCH -o counts.txt
#SBATCH -e counts.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com     # Your email for notifications
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH -A <project_number>     # Your HPC project or account ID

# Set the working directory to your project directory
cd /path/to/project/root

# Create a directory to store count matrix files
mkdir -p Counts

# Environment Setup
module load miniconda
conda activate RNA
module load subread/2.0.6     # Load Subread module

# List of SRA accessions corresponding to RNA-seq samples
accession=("SRR21949253" "SRR21949254" "SRR21949255" "SRR21949256" "SRR21949257" "SRR21949258" "SRR21949259" "SRR21949260" "SRR21949261" "SRR21949262" "SRR21949263" "SRR21949264" "SRR21949265" "SRR21949266" "SRR21949267" "SRR21949268" "SRR21949269" "SRR21949270")

# Path to the GTF annotation file
gtf_file="/path/to/project/root/genome/Sus_scrofa.Sscrofa11.1.113.chr.gtf"

# Run featureCounts in paired-end mode for each sample
for acc in "${accession[@]}"; do
	sorted_bam="align/$acc/$acc.sorted.bam"
	output="Counts/$acc_counts.txt"

	featureCounts -p -a "$gtf_file" -o "$output" "$sorted_bam"
done



