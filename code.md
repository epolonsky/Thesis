# Bioinformatics Workflow

This file shows the computational workflow used to process illumina short read sequencing data, generate a consensus genome assembly, assess assembly quality, and annotate the genome using Liftoff.

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
| `btk`       | BlobToolKit      |

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

Trim the first 20 bp from each read to remove the barcode sequence.

Java memory was increased using the `-Xms` and `-Xmx` options.

```bash
mkdir trimmomatic_output
trimmomatic PE -threads 2 -Xms1024m -Xmx8g -phred33 F7_PP_DNA_S15_R1_001.fastq.gz F7_PP_DNA_S15_R2_001.fastq.gz trimmomatic_output/R1_paired.fastq.gz trimmomatic_output/R1_unpaired.fastq.gz trimmomatic_output/R2_paired.fastq.gz trimmomatic_output/R2_unpaired.fastq.gz HEADCROP:20
```

---

## 3. Quality assessment after trimming

```bash
fastqc -o fastqc_output R1_paired.fastq.gz R2_paired.fastq.gz
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
mkdir -p bowtie2_output
bowtie2 -x bowtie2_index/ref_genome_index -1 trimmomatic_output/R1_paired.fastq.gz -2 trimmomatic_output/R2_paired.fastq.gz --un-conc-gz bowtie2_output/unmapped_reads -S bowtie2_output/aligned.sam -p 4 > bowtie2_output/bowtie2.log 2>&1
```

---

## 5. Convert aligned SAM files to sorted and indexed BAM

```bash
samtools view -bS bowtie2_output/aligned.sam -o aligned.bam
samtools sort bowtie2_output/aligned.bam -o bowtie2_output/aligned_sorted.bam
samtools index bowtie2_output/aligned_sorted.bam
```

---

## 6. Calculate sequencing coverage

```bash
samtools coverage bowtie2_output/aligned_sorted.bam > aligned_sorted_coverage_output.txt
```

---

## 7. Generate a reference-guided consensus genome 

Generate a reference-guided consensus sequence from reads aligned to the reference genome.

```bash
samtools consensus -m simple --min-depth 1 --call-fract 0.5 bowtie2_output/aligned_sorted.bam > consensus_simple.fasta
```

---

## 8. Assess assembly completeness

Evaluate genome completeness using BUSCO.

```bash
busco -i consensus_simple.fasta -l diptera_odb10 -o busco_out -m genome 
```

---

## 9.  Generate BlobToolKit assembly quality visualization snailplot using BUSCO results

```bash
blobtools create --fasta consensus_simple.fasta my_snail_dataset

blobtools add --busco full_table.tsv my_snail_dataset

blobtools view --view snail --plot --out ./ snail_plot
```

---

# Liftoff annotation transfer

## 10. Generate an annotation using Liftoff

Liftoff is used to project existing annotations from the *Culex pipiens pallens* reference genome assembly onto the `consensus.fasta` genome assembly. The resulting GFF3 annotation file contains the transferred positions of genes and other genomic features mapped onto the consensus assembly.

```bash
mkdir liftoff_simple
liftoff_simple -g GCF_016801865.2_pallens_genomic.gff -o liftoff_simple/consensus_simple_lifted_annotation.gff3 consensus_simple.fasta GCF_016801865.2_pallens_genomic.fna
```

---

# Insecticide resistant gene extraction

The CDS of the genes for *Culex pipiens pallens* and *Culex quinquefasciatus* were extracted from NCBI accession numbers found through literature search.

For *Culex pipiens pipiens* the genes were found by using the gene id locus of the pallens gene and extracting it from the annotation.

## 11. Extract CDS sequences from the Liftoff annotation

The CDS sequences from the Liftoff transferred annotation were extracted from the consensus_simple.fasta assembly using gffread. This produces a FASTA file containing the coding sequences (CDS) corresponding to the transferred gene annotations.

```bash
gffread liftoff_simple/consensus_simple_lifted_annotation.gff3 -g consensus_simple.fasta -x liftoff_simple/consensus_simple_CDS.fasta
```

## 12. Identify and extract genes of interest

Three insecticide resistance-associated genes were selected for comparison.

For Culex pipiens pipiens, the corresponding genes were identified using the gene locus IDs associated with the Culex pipiens pallens annotation transferred by Liftoff.

The following gene loci were extracted:

Gene	Gene locus
Ace-1	LOC120414010
Voltage-gated sodium channel	LOC120419138
GABA receptor	LOC120412863

A separate GFF3 file was created for each insecticide resistance-associated gene by filtering the Liftoff annotation for the corresponding gene locus.

```bash
mkdir pipiens_ins_res_genes
grep 'gene=LOC120414010' liftoff_simple/consensus_simple_lifted_annotation.gff3 > pipiens_ins_res_genes/pipiens_ace1.gff3
grep 'gene=LOC120419138' liftoff_simple/consensus_simple_lifted_annotation.gff3 > pipiens_ins_res_genes/pipiens_voltage_gated_sodium_channel.gff3
grep 'gene=LOC120412863' liftoff_simple/consensus_simple_lifted_annotation.gff3 > pipiens_ins_res_genes/pipiens_gaba.gff3
```

## 13. Extract the CDS sequences for each gene

The gene-specific GFF3 annotations were then used with gffread to extract the corresponding CDS sequences from the consensus_simple.fasta assembly

```bash
gffread pipiens_ins_res_genes/pipiens_ace1.gff3 -g consensus_simple.fasta -x pipiens_ins_res_genes/pipiens_ace1.fasta
gffread pipiens_ins_res_genes/pipiens_voltage_gated_sodium_channel.gff3 -g consensus_simple.fasta -x pipiens_ins_res_genes/pipiens_voltage_gated_sodium_channel.fasta
gffread pipiens_ins_res_genes/pipiens_gaba.gff3 -g consensus_simple.fasta -x pipiens_ins_res_genes/pipiens_gaba.fasta
```

## 15. Inspect the extracted FASTA sequences

The FASTA headers were inspected to identify the transcript/RNA accession corresponding to each gene of interest

```bash
grep '^>' pipiens_ins_res_genes/pipiens_ace1.fasta
grep '^>' pipiens_ins_res_genes/pipiens_voltage_gated_sodium_channel.fasta
grep '^>' pipiens_ins_res_genes/pipiens_gaba.fasta
```
The relevant transcript accessions were then used to extract the desired CDS sequence from each gene-specific FASTA file

## 16. Extract the selected CDS sequences

seqkit grep was used to select the specific transcript accession for each insecticide resistance-associated gene.

```bash
seqkit grep -n -p 'rna-XM_052710739.1' pipiens_ins_res_genes/pipiens_ace1.fasta > pipiens_ins_res_gene_cds/pipiens_ace1_CDS.fasta
seqkit grep -n -p 'rna-XM_039573487.2' pipiens_ins_res_genes/pipiens_gaba.fasta > pipiens_ins_res_gene_cds/pipiens_gaba_CDS.fasta
seqkit grep -n -p 'rna-XM_052707444.1' pipiens_ins_res_genes/pipiens_voltage_gated_sodium_channel.fasta > pipiens_ins_res_gene_cds/pipiens_voltage_gated_sodium_channel_CDS.fasta
```

These files contain the selected Culex pipiens pipiens CDS sequences for the three genes used in downstream comparisons.

## 17. Create pairwise multi-FASTA files

To compare the insecticide resistance-associated genes between species, the Culex pipiens pipiens CDS sequences were combined with the corresponding CDS sequences from Culex pipiens pallens and Culex quinquefasciatus.

Each resulting FASTA file contains two sequences for the same gene, allowing direct pairwise sequence comparison.

The following comparisons were generated:

- C. p. pallens vs. C. p. pipiens
- C. quinquefasciatus vs. C. p. pipiens

```
mkdir ins_res_gene_multi_fasta_files

# C. p. pallens vs. C. p. pipiens
cat pallens_ins_res_gene_cds/pallens_ace1.fasta pipiens_ins_res_gene_cds/pipiens_ace1_CDS.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_ace1.fasta
cat pallens_ins_res_gene_cds/pallens_gaba.fasta pipiens_ins_res_gene_cds/pipiens_gaba_CDS.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_gaba.fasta
cat pallens_ins_res_gene_cds/pallens_voltage_gated_sodium_channel.fasta pipiens_ins_res_gene_cds/pipiens_voltage_gated_sodium_channel_CDS.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_voltage_gated_sodium_channel.fasta

# C. quinquefasciatus vs. C. p. pipiens
cat quinx_ins_res_gene_cds/quinx_ace1.fasta pipiens_ins_res_gene_cds/pipiens_ace1_CDS.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_ace1.fasta
cat quinx_ins_res_gene_cds/quinx_gaba.fasta pipiens_ins_res_gene_cds/pipiens_gaba_CDS.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_gaba.fasta
cat quinx_ins_res_gene_cds/quinx_voltage_gated_sodium_channel.fasta pipiens_ins_res_gene_cds/pipiens_voltage_gated_sodium_channel_CDS.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_voltage_gated_sodium_channel.fasta
```
---

# Annotation statistics

## Generate annotation stats with AGAT

Annotation statistics were generated from the Liftoff-generated GFF3 annotation using AGAT:

```bash
agat_sp_statistics.pl --gff consensus_lifted_annotation.gff3 -o liftoff_simple/liftoff_simple_annotation_gff_statistics.txt
```

The AGAT output file contains statistics for each annotation feature type (mRNA, lncRNA, rRNA, tRNA, snRNA, snoRNA, and transcript). For protein-coding gene annotation statistics, values were taken from the following section:

------------------------------------- mrna ---------------------------
mrna have isoforms! Here are the statistics without isoforms shortest isoforms excluded)

This section was used because it reports one representative transcript per gene.

The following annotation statistics were obtained from the AGAT output:

| Statistic | AGAT output field |
|-----------|-------------------|
| Number of protein-coding genes | `Number of gene` |
| Number of exons per gene | `mean exons per mrna` |
| Mean gene length (bp) | `mean gene length (bp)` |
| Mean exon length (bp) | `mean exon length (bp)` |
| Number of CDSs per gene | `mean cdss per mrna` |
| Mean CDS length (bp) | `mean cds length (bp)` |
| Number of introns per gene | `mean introns in cdss per mrna` |
| Mean intron length (bp) | `mean intron in cds length (bp)` |
| Total gene length (bp) | `Total gene length (bp)` |
| Total exon length (bp) | `Total exon length (bp)` |
| Total CDS length (bp) | `Total cds length (bp)` |
| Total intron length (bp) | `Total intron length per cds (bp)` |

## Generate predicted protein sequences

Protein sequences were extracted from the Liftoff annotation using gffread:

```bash
gffread consensus_simple_lifted_annotation.gff3 -g consensus_simple.fasta -y consensus_simple_proteins.fasta
```

## Number of predicted protein sequences

The number of predicted protein sequences was calculated by counting protein FASTA records:

```bash
grep -c "^>" proteins.fa
```

## Mean protein length (aa)

Mean protein length was calculated from the extracted protein FASTA file:

```bash
seqkit fx2tab -nl consensus_simple_proteins.fasta | \
awk '{sum+=length($2); n++} END {print sum/n}'
```

## Genome size

Genome assembly size was calculated from the genome FASTA:

```bash
seqkit stats simple_consensus.fasta
```

The sum_len value was used as the total genome size.

## Genome, exon, CDS, and intron rations

Feature density ratios were calculated as:

```bash
feature length / genome assembly size × 100
```

using the following AGAT values:

Gene ratio:
Total gene length (Total gene length (bp)) / genome size
Exon ratio:
Total exon length (Total exon length (bp)) / genome size
CDS ratio:
Total CDS length (Total cds length (bp)) / genome size
Intron ratio:
Total intron length (Total intron length per cds (bp)) / genome size

---

# Genome Survey Analysis

##  K-mer counting using Jellyfish

Genome characteristics were estimated from trimmed Illumina reads using a k-mer based approach.  
21-mer frequencies were generated using Jellyfish with canonical k-mers (`-C`). The k-mer size was selected as **k=21**, which is commonly used for short-read genome size estimation.

The hash size (`-s`) was set to **500M**.

```bash
mkdir jellyfish_output
jellyfish count -C -m 21 -s 500M -t 20 -o jellyfish_output/F7_k21_trimmed.jf <(zcat R1_paired.fastq.gz) <(zcat R2_paired.fastq.gz) > jellyfish_output/jellyfish.log 2>&1 
```

## Generate k-mer histogram

The Jellyfish k-mer database was converted into a k-mer frequency histogram for downstream genome profiling.

```bash
jellyfish histo -t 20 jellyfish_output/F7_k21_trimmed.jf > jellyfish_output/F7_k21_trimmed.histo
```

## GenomeScope 2.0 genome survey analysis (online web version)

GenomeScope 2.0 was used to estimate genome characteristics from the k-mer frequency distribution.

The analysis was performed using a diploid model (`p=2`) with the following parameters:

| Parameter | Value | Description |
|-----------|-------|-------------|
| k-mer size (`k`) | 21 | Length of k-mer used for genome profiling |
| Ploidy (`p`) | 2 | Diploid genome model |
| Maximum k-mer coverage (`m`) | 1,000,000 | Maximum coverage threshold included in model fitting |

GenomeScope was run using the Jellyfish-generated histogram as input.

