# Reference-Guided Genome Assembly and Annotation of the Northern Illinois *Culex pipiens pipiens* Mosquito

## Overview

This repository contains the computational workflow developed for my MS thesis project. The workflow process Illumina sequencing data, generate a reference-guided consensus genome, evaluates genome completeness, and perform structural and functional genome annotation.

The workflow code and commands are documented in `code.md`.

The workflow includes:

1. Raw sequencing read quality assessment
2. Barcode trimming and post-trimming quality evaluation
3. Reference-guided read mapping
4. Identification and characterization of unmapped reads
5. Reference-guided consensus genome generation
6. Genome completeness assessment
7. Genome quality visualization
8. Repeat identification and masking
9. Evidence-supported genome prediction using BRAKER4
10. Functional annotation using eggNOG
11. Genome and annotation statistics
12. Transfer of gene annotations to the consensus genome using Liftoff
13. Extraction and comparison of insecticide-resistance-associated genes

The resulting consensus genome and annotations are intended to support downstream genomic analyses of *Culex pipiens pipiens*.

---

## Project Information

**Organism:** *Culex pipiens pipiens*  
**Project type:** Bioinformatics MS Thesis  
**Primary objective:** Generate and annotate a reference-guided consensus genome for a Northern Illinois *Culex pipiens pipiens* mosquito using Illumina sequencing data.

---
**Reference genome:** Culex pipiens, idCulPipi1.1

**NCBI assembly accession:** GCA_963924435.1

**Associated publication:** https://doi.org/10.12688/wellcomeopenres.23767.1

---

The workflow used Conda-managed environments for software installation and dependency management. The table below lists the principal software tools and versions recorded for the analysis.

| Software          | Version     |
| ----------------- | ----------- |
| Conda             | 26.3.2      |
| FastQC            | 0.12.1      |
| Trimmomatic       | 0.40        |
| Bowtie2           | 2.5.4       |
| SAMtools          | 1.19.2      |
| BLAST+            | 2.17.0      |
| Biopython         | 1.87        |
| Seqtk             | 1.5-r133    |
| seqkit            | 2.13.0      |
| gffread           | 0.12.7      |
| Jellyfish         | 2.2.10      |
| GenomeScope       | 2.0         |
| BUSCO             | 6.1.0       |
| BlobToolKit       | 4.4.5       |
| RepeatModeler     | 2.0.9       |
| RepeatMasker      | 4.2.3       |
| BRAKER4 pipeline  | 0.5.0-beta  |
| AUGUSTUS          | 3.5.0       |
| Liftoff           | 1.5.2       |
| eggNOG-mapper     | 2.1.13      |

Additional Pipeline Dependencies

Some tools invoke other software internally. Important versions recorded by the BRAKER4 run or reported by eggNOG-mapper include:

| Software       | Version                | Associated workflow |
| -------------- | ---------------------- | ------------------- |
| DIAMOND        | 2.0.15                 | BRAKER4             |
| DIAMOND        | 2.2.8                  | eggNOG-mapper       |
| ProtHint       | 2.6.0                  | BRAKER4             |
| miniprot       | 0.12-r237              | BRAKER4             |
| Spaln          | 2.3.3f                 | BRAKER4             |
| TSEBRA         | Commit `23ef205`       | BRAKER4             |
| GeneMark-EP+   | 4.*; commit `aacc025`  | BRAKER4             |
| GeneMark-ES    | 4.*; commit `aacc025`  | BRAKER4             |
| compleasm      | 0.2.8                  | BRAKER4             |
| MMseqs2        | 18.8cc5c               | eggNOG-mapper       |
| Apptainer      | 1.5.4                  | BRAKER4 container execution |
| Snakemake      | 8.18.2                 | BRAKER4 environment |
| Snakemake      | 7.19.1                 | BlobToolKit environment |
| RMBlast        | 2.14.1+                | RepeatMasker|
| AGAT           | 1.4.1                  | BRAKER4 dependency |
| AUGUSTUS       | 3.5.0                  | BRAKER4|
| BUSCO          | 6.1.0                  | Genome completeness assessment; BRAKER4 dependency |

---

# Conda Environments

Separate Conda environments were used to reduce dependency conflicts between software packages.

| Environment      | Purpose                                  |
| ---------------- | ---------------------------------------- |
| `blast_env`      | Protein database construction and BLAST searches |
| `busco`          | Genome completeness assessment           |
| `btk`            | Genome quality visualization             |
| `repeatmod_env`  | Repeat library generation                |
| `repeatmask_env` | Repeat masking                           |
| `braker4_env`    | BRAKER4 gene prediction                  |
| `eggnog_env`     | Functional annotation                    |

---
