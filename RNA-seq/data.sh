#!/bin/bash
# ----------------------------------------- SLURM Job Settings ----------------------------------------- #
#SBATCH -J data.sh
#SBATCH -p general        # Partition name (check with your HPC)
#SBATCH -o data.txt
#SBATCH -e data.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com     # Your email for notifications
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH -A <project_number>     # Your HPC project or account ID

# Set the working directory to your project directory
cd /path/to/project/root

# Environment Setup
module load miniconda
conda activate RNA
module load sra-toolkit/3.0.5     # Load SRA Toolkit module

# Download RNA-seq data from NCBI SRA
fasterq-dump --split-files SRR21949253 SRR21949254 SRR21949255 SRR21949256 SRR21949257 SRR21949258 SRR21949259 SRR21949260 SRR21949261 SRR21949262 SRR21949263 SRR21949264 SRR21949265 SRR21949266 SRR21949267 SRR21949268 SRR21949269 SRR21949270

