#!/bin/bash
#SBATCH -J index.sh
#SBATCH -p general        # Partition name (check with your HPC)
#SBATCH -o index.txt
#SBATCH -e index.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com     # Your email for notifications
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH -A <project_number>     # Your HPC project or account ID

# Set the working directory to your project directory
cd /path/to/project/root

# Environment Setup
module load miniconda
conda activate RNA
module load star/2.7.11a      # Load STAR module

# Download Reference Genome GTF file
mkdir genome && cd genome
wget https://ftp.ensembl.org/pub/current_gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.113.chr.gtf.gz
gunzip Sus_scrofa.Sscrofa11.1.113.chr.gtf.gz

wget https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
gunzip Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz

# Indexing
STAR --runMode genomeGenerate \
	--genomeDir /path/to/project/root/genome \
	--genomeFastaFiles /path/to/project/root/genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	--sjdbGTFfile /path/to/project/root/genome/Sus_scrofa.Sscrofa11.1.113.chr.gtf \
	--runThreadN 10
