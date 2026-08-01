# WGS_CRKP_HYBRID Pipeline

## 📖 Overview
This repository contains a highly scalable, custom bioinformatics pipeline developed to determine the role of increased *bla*OXA-48 copy numbers and OmpK36 loss in carbapenem resistance in *Klebsiella pneumoniae*. 

The workflow was developed as part of an internship project at **Precision Gene Joint Stock Company (Precigene)** in collaboration with the **School of Biotechnology, International University, VNU-HCM**. The analytical pipeline was applied to investigate genomic and proteogenomic data originally described by Meekes et al. (Antimicrobial Agents and Chemotherapy, 2025).

## ⚙️ Bioinformatics Pipeline Modules
The pipeline integrates seven sequential modules operating entirely within a Linux command-line environment:

*   **Module 1 (Pre-processing & QC):** Short-read quality filtering using FastQC and fastp; long-read assessment and filtering using NanoPlot and Filtlong.
*   **Module 2 (Hybrid Assembly & Validation):** *De novo* hybrid assembly via Unicycler; assembly quality evaluation using QUAST, BUSCO, and CheckM2.
*   **Module 3 & 4 (High-resolution Annotation & Functional Screening):** Protein-coding annotation by Bakta; non-coding RNA identification via Infernal against Rfam; epidemiological typing and AMR/virulence screening using Kleborate, AMRFinderPlus, StarAMR, VFDB, and MOB-suite/PlasmidFinder.
*   **Module 5 (Deep-Dive Mechanisms):** Plasmid Copy Number (PCN) quantification (using BWA-MEM and Samtools), porin integrity analysis, and Insertion Sequence (IS) element detection via ISEScan.
*   **Module 6 (Variant Calling & Effect):** SNP/Indel calling and functional annotation using FreeBayes/Snippy and SnpEff (with a specific focus on regulatory mutations like *copA*).
*   **Module 7 (Phylogenomics & Pangenome):** Core-genome alignment and maximum-likelihood phylogenetic tree construction using Panaroo, Gubbins, and IQ-TREE.

## 📂 Repository Structure
*   **`code/`**: The core directory containing all scripts.
    *   `HYBRIDS_PIPELINE/`: Sequential Shell scripts (`00` to `10`) for automated data processing, QC, hybrid assembly, annotation, and variant calling.
    *   `DATA_VISUALIZATION/`: R scripts used to generate the statistical figures (e.g., PCN amplification, resistome profiling heatmaps, porin status) presented in the final report.
*   **`Genomic_Pipeline_Report_Final.docx.pdf`**: A comprehensive final report detailing the methodology, results, and discussion of the WGS_CRKP_HYBRID pipeline. It provides in-depth assessments of the combined impact of increased *bla*OXA-48 or *bla*CTX-M-15 plasmid copy numbers and OmpK36 porin inactivation on carbapenem resistance.
*   **`environment.yml`**: Conda environment configuration file to reproduce the exact software dependencies used in this study.
*   **`sample_metadata.tsv`** & **`accessions.csv`**: Detailed metadata and accession numbers for the 15 clinical isolates and 4 laboratory-derived mutants analyzed
*   **`pairs.tsv`**: Configuration file defining native-mutant isolate pairs for comparative mutational and PCN analysis.

## 🚀 Usage
To ensure full reproducibility, clone this repository and set up the Conda environment using the provided configuration file:

```bash
# Clone the repository
git clone [https://github.com/nguyenanhduc11032005-creator/WGS_CRKP_HYBRID_pipeline.git](https://github.com/nguyenanhduc11032005-creator/WGS_CRKP_HYBRID_pipeline.git)
cd WGS_CRKP_HYBRID_pipeline

# Create and activate the conda environment
conda env create -f environment.yml
conda activate [replace_with_your_env_name]
