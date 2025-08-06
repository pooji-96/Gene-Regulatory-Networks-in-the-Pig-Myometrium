#!/bin/bash
#SBATCH -J bamindex.sh
#SBATCH -p general        # Partition name (check with your HPC)
#SBATCH -o bamindex.txt
#SBATCH -e bamindex.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com     # Your email for notifications
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH -A <project_number>     # Your HPC project or account ID

# Set the working directory to your project directory
cd /path/to/project/root

# Environment Setup
module load miniconda
conda activate RNA
module load samtools/1.17     # Load Samtools module

# List of SRA accessions corresponding to RNA-seq samples
accession=("SRR21949253" "SRR21949254" "SRR21949255" "SRR21949256" "SRR21949257" "SRR21949258" "SRR21949259" "SRR21949260" "SRR21949261" "SRR21949262" "SRR21949263" "SRR21949264" "SRR21949265" "SRR21949266" "SRR21949267" "SRR21949268" "SRR21949269" "SRR21949270")

# Sort, index, and validate BAM files for each sample
for acc in "${accession[@]}"; do
	bam_file="align/$acc/Aligned.sortedByCoord.out.bam"
	sorted_bam="align/$acc/$acc.sorted.bam"

	samtools sort "$bam_file" -o "$sorted_bam"
	samtools index "$sorted_bam"
    samtools quickcheck "$sorted_bam"

done


