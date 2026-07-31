# Reference-Guided Genome Assembly and Annotation of the Northern Illinois *Culex pipiens pipiens* Mosquito

## Overview

This repository contains the computational workflow used for my MS thesis project to process Illumina sequencing data, generate a reference-guided consensus genome assembly, evaluate assembly quality, and perform structural genome annotation using MAKER.

The workflow includes:

1. Raw sequencing read quality assessment
2. Adapter/barcode trimming and post-trimming quality evaluation
3. Reference-guided read mapping
4. Extraction and characterization of unmapped reads
5. Consensus genome generation
6. Assembly completeness assessment
7. Genome quality visualization
8. Repeat identification and masking
9. Evidence-supported genome annotation using MAKER

The pipeline was designed to produce a high-quality genome assembly and annotation for downstream genomic analyses.

---

## Project Information

**Organism:** *Culex pipiens pipiens*  
**Project type:** Bioinformatics MS Thesis  
**Primary objective:** Generate and annotate a genome assembly for the *Culex pipiens pipiens* mosquito species from Illumina sequencing data 

---
Reference Genome NCBI: idCulPipi1.1, genome accession GCA_963924435.1

Publication: https://doi.org/10.12688/wellcomeopenres.23767.1

---

The workflow was performed using Conda-managed environments to maintain software compatibility and reproducibility.

| Software | Version |
|----------|---------|
| Conda | 26.3.2 |
| Trimmomatic | 0.40 |
| FastQC | 0.12.1 |
| Bowtie2 | 2.5.5 |
| SAMtools | 1.19.2 |
| Seqtk | 1.5-r133 |
| BLAST+ | 2.17.0 |
| GATK | 4.6.2.0 |
| BUSCO | 6.1.0 |
| BlobToolKit | 4.4.5 |
| MAKER | 3.01.04 |
| AUGUSTUS | 3.5.0 |
| RepeatMasker | 4.2.3 |
| RepeatModeler | 2.0.8 |

---

# Conda Environments

Software packages were installed into separate Conda environments to avoid dependency conflicts.

| Environment | Purpose |
|-------------|---------|
| `blast_env` | BLAST+ |
| `busco` | BUSCO genome completeness assessment |
| `gatk_env` | GATK analysis |
| `btk` | BlobToolKit visualization |
| `maker` | MAKER genome annotation |

---
