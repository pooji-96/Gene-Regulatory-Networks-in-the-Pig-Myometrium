# RNA-seq data analysis

## Description:
   This project involves the analysis of RNA-seq data for the uterine circular muscle (CM) and longitudinal muscle (LM) layers at the mesometrial (M) side on days 12 (pre-implantation stage) and 15 (implantation stage) of pregnancy and day 15 of the estrous cycle in pigs. The RNA-seq data is downloaded and processed using bash scripts and Unix commands within a conda environment. Quality control is performed with FastQC, followed by adapter trimming with TrimGalore. High-quality clean reads are aligned with STAR and quantification is done using featureCounts.

### Required files:
The sequencing data for this analysis is sourced from the NCBI Sequence Read Archive (SRA) (Project ID: PRJNA891690) and downloaded using the SRA Toolkit.

### Software Dependencies  
**Prerequisite:** Access to a Slurm-based HPC system with the required software tools available or installed manually (if applicable)  
`miniconda`  
`sratoolkit`  
`fastqc`  
`python`  
`cutadapt`  
`trim_galore`  
`star`  
`subread` 


### Directory Structure  
```plaintext
/path/to/project/root/
├── *.sh                 # Scripts
├── *_1.fastq.gz         # Raw paired-end reads (R1)
├── *_2.fastq.gz         # Raw paired-end reads (R2)
├── genome/              # Reference genome and annotation file
├── fastqc/              # FASTQC reports
├── trim/                # Trimmed FASTQ files
├── align/SRR*           # STAR alignment outputs (sorted BAM files)
├── Counts/              # Gene count output files
├── *.txt                # Output logs or summaries
└── *.err                # Error logs from script runs
```


### Execution:
1. Change path to your project directory
```bash
    cd /path/to/project/root
```

2. Create and activate conda environment
```bash
    module load miniconda
    conda create --name RNA
    conda activate RNA
```

⚠️**Note:** Before running, update all file paths and Slurm settings (`#SBATCH` lines) in the `.sh` scripts as needed.

3. Download SRA dataset
```bash
    sbatch data.sh
```

4. Load FastQC and run Quality Control
```bash
    sbatch qc.sh
```

5. Load TrimGalore and Perform Adapter Trimming
```bash
    sbatch trim.sh
```

6. Download Reference genome and annotations & Genome Indexing with STAR
```bash
    sbatch index.sh
```

7. Read alignment with STAR
```bash
    sbatch align.sh
```

8. Sorting and Indexing BAM Files for Quantification
```bash
    sbatch bamindex.sh
```

9. Quantification with featureCounts
```bash
    sbatch counts.sh
```


### Output
Raw FASTQ files: SRR*/  
Quality reports: fastqc/  
Trimmed FASTQ files: trim/  
Aligned & Sorted BAM files: align/SRR*/    
Gene count tables: Counts/  

---

> Author: Poojitha Kolli  
> Completed: December 17th 2024

