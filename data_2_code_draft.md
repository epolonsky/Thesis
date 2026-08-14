# Bioinformatics Workflow

This file shows the computational workflow used to process sequencing data, generate a consensus genome assembly, assess assembly quality, and annotate the genome using MAKER.

---

# done in the data_2 directory
cd ~/data_2

---

# Software

| Software      | Version  |
| ------------- | -------- |
| Conda         | 26.3.2   |
| Trimmomatic   | 0.40     |
| FastQC        | 0.12.1   |
| Bowtie2       | 2.5.5    |
| SAMtools      | 1.19.2   |
| Seqtk         | 1.5-r133 |
| BLAST+        | 2.17.0   |
| Jellyfish     | 2.2.10   |
| Liftoff       | 1.5.1    |
| BUSCO         | 6.1.0    |
| BlobToolKit   | 4.4.5    |
| MAKER         | 3.01.04  |
| AUGUSTUS      | 3.5.0    |
| RepeatMasker  | 4.2.3    |
| RepeatModeler | 2.0.8    |
| AGAT          | 1.7.0    |
| seqkit        | 2.13.0   |
| genomescope   | 2.0      | 

---

# Conda Environments

These software packages were installed in separate Conda environments to maintain compatible software dependencies.

| Environment | Purpose          |
| ----------- | ---------------- |
| `blast_env` | BLAST            |
| `busco`     | BUSCO            |
| `gatk_env`  | GATK             |
| `btk`       | BlobToolKit      |
| `maker`     | MAKER annotation |

---

# Workflow

## 1. Quality assessment of raw sequencing reads

Assess raw Illumina read quality using FastQC.

```bash
mkdir fastqc_output
fastqc -o fastqc_output F7_PP_DNA_S15_R1_001.fastq.gz F7_PP_DNA_S15_R2_001.fastq.gz
```

---

## 2. Trim sequencing reads

Trim the first 20 bp from each read to remove barcode sequence.

Java memory was increased using the `-Xms` and `-Xmx` options.

```bash
mkdir trimmomatic_output
trimmomatic PE -threads 16 -Xms1024m -Xmx8g -phred33 F7_PP_DNA_S15_R1_001.fastq.gz F7_PP_DNA_S15_R2_001.fastq.gz R1_trimmed_paired.fastq.gz R1_trimmed_unpaired.fastq.gz R2_trimmed_paired.fastq.gz R2_trimmed_unpaired.fastq.gz HEADCROP:20 CROP:124 SLIDINGWINDOW:4:20 MINLEN:50
```

---

## 3. Quality assessment after trimming

```bash
fastqc -o fastqc_output R1_trimmed_paired.fastq.gz  R2_trimmed_paired.fastq.gz
```

---

## 4. Map reads to the Culex pipiens reference genome

Reads were mapped against the NCBI chromosome-level assembly:

*Culex pipiens* idCulPipi1.1  
Accession: GCA_963924435.1

Build a Bowtie2 index and align paired-end reads to the reference genome.

```bash
mkdir bowtie2_index
bowtie2-build -f GCA_963924435.1_idCulPipi1.1_genomic.fna bowtie2_index/ref_genome_index
```

```bash
mkdir -p bowtie2_output/unmapped_reads
bowtie2 -x bowtie2_index/ref_genome_index -1 R1_trimmed_paired.fastq.gz -2 R2_trimmed_paired.fastq.gz --un-conc-gz bowtie2_output/unmapped_reads/unmapped_reads -S bowtie2_output/aligned.sam -p 4 > bowtie2_output/bowtie2.log 2>&1
```

---

## 5. Convert aligned SAM files to sorted and indexed BAM

```bash
samtools view -bS bowtie2_output/aligned.sam
samtools sort bowtie2_output/aligned.bam -o bowtie2_output/aligned_sorted.bam
samtools index bowtie2_output/aligned_sorted.bam
```

---

## 6. Calculate sequencing coverage

```bash
samtools coverage bowtie2_output/aligned.sam > aligned_sam_coverage.txt
```

---

## 7. Generate a reference-guided consensus genome 

Generate a reference-guided consensus sequence from reads aligned to the reference genome.

```bash
samtools consensus -f fasta bowtie2_output/aligned_sorted.bam -o consensus.fasta
```

---

## 8. Assess assembly completeness

Evaluate genome completeness using BUSCO.

```bash
busco -i consensus.fasta -l diptera_odb10 -o busco_out -m genome 
```

---


# Genome Survey Analysis

##  K-mer counting using Jellyfish

Genome characteristics were estimated from trimmed Illumina reads using a k-mer based approach.  
21-mer frequencies were generated using Jellyfish with canonical k-mers (`-C`). The k-mer size was selected as **k=21**, which is commonly used for short-read genome size estimation.

The hash size (`-s`) was set to **500M**.

```bash
mkdir jellyfish_output
jellyfish count -C -m 21 -s 500M -t 20 -o jellyfish_output/F7_k21_trimmed.jf <(zcat R1_trimmed_paired.fastq.gz) <(zcat R2_trimmed_paired.fastq.gz) > jellyfish_output/jellyfish.log 2>&1 
```

## Generate k-mer histogram

The Jellyfish k-mer database was converted into a k-mer frequency histogram for downstream genome profiling.

```bash
jellyfish histo -t 20 jellyfish_output/F7_k21_trimmed.jf > jellyfish_output/F7_k21_trimmed.histo
```

## genomescope 2 genome survey

Run the genome survey analysis using genomescope
```bash
genomescope2 -i jellyfish_output/F7_k21_trimmed.histo -o genomescope_output -k 21 -p 2 -m 1000000
```

