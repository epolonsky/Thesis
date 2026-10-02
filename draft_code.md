Swedish Genome Annotation codes

# Install and run repeatmodeler 2

```bash
cd ~/data
conda create -n repeatmod_env -c conda-forge -c bioconda repeatmodeler=2.0.9
conda activate repeatmod_env

BuildDatabase -name idCulPipi1_genome_db GCA_963924435.1_idCulPipi1.1_genomic.fna
nohup RepeatModeler -database idCulPipi1_genome_db -threads 8 -dir RM_idCulPipi1 & # change RM to repeat_modled when done runnning
```

# Install and run repeatmasker

```bash
cd ~/data
conda create -n repeatmask_env -c conda-forge -c bioconda repeatmasker=4.2.4
conda activate repeatmask_env

nohup RepeatMasker -pa 4 -gff -xsmall -lib repeat_modled_idCulPipi1/consensi.fa -dir repeat_masked_idCulPipi1 GCA_963924435.1_idCulPipi1.1_genomic.fna &
```













# Bioinformatics Workflow

This file shows the computational workflow used to process illumina short read sequencing data, generate a consensus genome assembly, assess assembly quality, and annotate the genome using Liftoff.

## Assess raw Illumina read quality with FastQC

```bash
mkdir fastqc_output
fastqc -o fastqc_output F7_PP_DNA_S15_R1_001.fastq.gz F7_PP_DNA_S15_R2_001.fastq.gz
```

## Trim sequencing reads with Trimmomatic

Trim the first 20 bp from each read to remove the barcode sequence.

Java memory was increased using the `-Xms` and `-Xmx` options.

```bash
mkdir trimmomatic_output
trimmomatic PE -threads 2 -Xms1024m -Xmx8g -phred33 F7_PP_DNA_S15_R1_001.fastq.gz F7_PP_DNA_S15_R2_001.fastq.gz trimmomatic_output/R1_paired.fastq.gz trimmomatic_output/R1_unpaired.fastq.gz trimmomatic_output/R2_paired.fastq.gz trimmomatic_output/R2_unpaired.fastq.gz HEADCROP:20
```

 ## Asses trimmed read quality with FastQC

```bash
fastqc -o fastqc_output R1_paired.fastq.gz R2_paired.fastq.gz
```

## Map reads to the Culex pipiens reference genome with Bowtie2

Reads were mapped against the NCBI chromosome-level assembly:

*Culex pipiens* idCulPipi1.1  
Accession: GCA_963924435.1

```bash
mkdir bowtie2_index
bowtie2-build -f GCA_963924435.1_idCulPipi1.1_genomic.fna bowtie2_index/ref_genome_index
```

```bash
mkdir -p bowtie2_output
bowtie2 -x bowtie2_index/ref_genome_index -1 trimmomatic_output/R1_paired.fastq.gz -2 trimmomatic_output/R2_paired.fastq.gz --un-conc-gz bowtie2_output/unmapped_reads -S bowtie2_output/aligned.sam -p 4 > bowtie2_output/bowtie2.log 2>&1
```

## Convert aligned SAM files to sorted and indexed BAM with Samtools

```bash
samtools view -bS bowtie2_output/aligned.sam -o aligned.bam
samtools sort bowtie2_output/aligned.bam -o bowtie2_output/aligned_sorted.bam
samtools index bowtie2_output/aligned_sorted.bam
```

## Calculate sequencing coverage with Samtools

```bash
samtools coverage bowtie2_output/aligned_sorted.bam > aligned_sorted_coverage_output.txt
```

## Generate a reference-guided consensus genome from reads aligned to the reference genome with Samtools consensus

```bash
samtools consensus -m simple --min-depth 1 --call-fract 0.5 bowtie2_output/aligned_sorted.bam > consensus_simple.fasta
```

## Evaluate genome completeness with BUSCO

```bash
busco -i consensus_simple.fasta -l diptera_odb10 -o busco_out -m genome 
```

##  Generate snailplot from BUSCO results with BlobToolKit

```bash
cd busco_simple_out

blobtools create --fasta ~/data/consensus_simple.fasta my_snail_dataset

blobtools add --busco run_diptera_odb10/full_table.tsv my_snail_dataset

blobtools view --view snail --plot --out ./ my_snail_dataset
```

---

## Transfer an existing *Culex pipiens pallens* reference genome annotation to my consensus sequence with Liftoff

```bash
mkdir liftoff_simple
liftoff_simple -g GCF_016801865.2_pallens_genomic.gff -o liftoff_simple/consensus_simple_lifted_annotation.gff3 consensus_simple.fasta GCF_016801865.2_pallens_genomic.fna
```

---

# Insecticide resistant gene extraction

The CDS of the genes for *Culex pipiens pallens* and *Culex quinquefasciatus* were extracted from NCBI accession numbers found through literature search.

For *Culex pipiens pipiens* the genes were found by using the gene id locus of the pallens gene and extracting it from the annotation.

## Extract CDS sequences from the Liftoff annotation with gffread

```bash
gffread liftoff_simple/consensus_simple_lifted_annotation.gff3 -g consensus_simple.fasta -x liftoff_simple/consensus_simple_CDS.fasta
```

## Identify and extract genes of interest with Samtools

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
samtools faidx consensus_simple.fasta OZ004313.1:67389111-67610679 > pipiens_ins_res_genes/pipiens_ace1.fasta
samtools faidx consensus_simple.fasta OZ004312.1:191453602-191579832 > pipiens_ins_res_genes/pipiens_sodium_channel.fasta
samtools faidx consensus_simple.fasta OZ004313.1:172172965-172310683 > pipiens_ins_res_genes/pipiens_rdlbagaba.fasta
```

## Extract the CDS sequences for each gene with seqkit

seqkit grep was used to select the specific transcript accession for each insecticide resistance-associated gene.

```bash
mkdir pipiens_ins_res_genes/pipiens_cds

seqkit grep -p 'rna-XM_052710737.1' liftoff_simple/consensus_simple_CDS.fasta > pipiens_ins_res_genes/pipiens_cds/pipiens_ace1_cds.fasta

seqkit grep -p 'rna-XM_052707395.1' liftoff_simple/consensus_simple_CDS.fasta > pipiens_ins_res_genes/pipiens_cds/pipiens_sodium_channel_cds.fasta

seqkit grep -p 'rna-XM_039573482.2' liftoff_simple/consensus_simple_CDS.fasta > pipiens_ins_res_genes/pipiens_cds/pipiens_rdlbagaba_cds.fasta
```

These files contain the selected Culex pipiens pipiens CDS sequences for the three genes used in downstream comparisons.

## Create pairwise multi-FASTA files for the gene DNA sequences

To compare the insecticide resistance-associated genes between species, the Culex pipiens pipiens CDS sequences were combined with the corresponding CDS sequences from Culex pipiens pallens and Culex quinquefasciatus.

Each resulting FASTA file contains two sequences for the same gene, allowing direct pairwise sequence comparison.

The following comparisons were generated:

- C. p. pallens vs. C. p. pipiens
- C. quinquefasciatus vs. C. p. pipiens

```bash
cd ~/data

mkdir ins_res_gene_multi_fasta_files

# C. p. pallens vs. C. p. pipiens
cat quinx_pallens_genes/quinx_pallens_cds/pallens_ace1_cds.fasta pipiens_ins_res_genes/pipiens_cds/pipiens_ace1_cds.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_ace1.fasta

cat quinx_pallens_genes/quinx_pallens_cds/pallens_rdl_gaba_cds.fasta pipiens_ins_res_genes/pipiens_cds/pipiens_rdl_gaba_cds.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_rdl_gaba.fasta

cat quinx_pallens_genes/quinx_pallens_cds/pallens_sodium_channel_cds.fasta pipiens_ins_res_genes/pipiens_cds/pipiens_sodium_channel_cds.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_sodium_channel.fasta

# C. quinquefasciatus vs. C. p. pipiens
cat quinx_pallens_genes/quinx_pallens_cds/quinx_ace1_cds.fasta pipiens_ins_res_genes/pipiens_cds/pipiens_ace1_cds.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_ace1.fasta

cat quinx_pallens_genes/quinx_pallens_cds/quinx_rdl_gaba_cds.fasta pipiens_ins_res_genes/pipiens_cds/pipiens_rdl_gaba_cds.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_rdl_gaba.fasta

cat quinx_pallens_genes/quinx_pallens_cds/quinx_sodium_channel_cds.fasta pipiens_ins_res_genes/pipiens_cds/pipiens_sodium_channel_cds.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_sodium_channel.fasta
```

## Translate extracted CDS sequences from *C. p. pipiens* to their respective amino acid sequence using biopython

```bash
cd pipiens_ins_res_genes

python - <<'PY'
from Bio import SeqIO
from pathlib import Path

input_dir = Path("pipiens_cds")
output_dir = Path("pipiens_prot")
output_dir.mkdir(exist_ok=True)

for cds_file in input_dir.glob("*.fasta"):
    protein_file = output_dir / cds_file.name.replace("_cds.fasta", "_prot.fasta")

    with open(protein_file, "w") as out:
        for record in SeqIO.parse(cds_file, "fasta"):
            record.seq = record.seq.translate()
            SeqIO.write(record, out, "fasta")

    print(f"Translated: {cds_file} -> {protein_file}")
PY
```

## Create pairwise multi-FASTA files for the protein sequences

```bash
cd ~/data

# C. p. pallens vs. C. p. pipiens
cat quinx_pallens_genes/quinx_pallens_prot/pallens_ace1_prot.fasta pipiens_ins_res_genes/pipiens_prot/pipiens_ace1_prot.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_ace1_prot.fasta

cat quinx_pallens_genes/quinx_pallens_prot/pallens_rdl_gaba_prot.fasta pipiens_ins_res_genes/pipiens_prot/pipiens_rdl_gaba_prot.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_rdl_gaba_prot.fasta

cat quinx_pallens_genes/quinx_pallens_prot/pallens_sodium_channel_prot.fasta pipiens_ins_res_genes/pipiens_prot/pipiens_sodium_channel_prot.fasta > ins_res_gene_multi_fasta_files/pallens_vs_pipiens_sodium_channel_prot.fasta

# C. quinquefasciatus vs. C. p. pipiens
cat quinx_pallens_genes/quinx_pallens_prot/quinx_ace1_prot.fasta pipiens_ins_res_genes/pipiens_prot/pipiens_ace1_prot.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_ace1_prot.fasta

cat quinx_pallens_genes/quinx_pallens_prot/quinx_rdl_gaba_prot.fasta pipiens_ins_res_genes/pipiens_prot/pipiens_rdl_gaba_prot.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_rdl_gaba_prot.fasta

cat quinx_pallens_genes/quinx_pallens_prot/quinx_sodium_channel_prot.fasta pipiens_ins_res_genes/pipiens_prot/pipiens_sodium_channel_prot.fasta > ins_res_gene_multi_fasta_files/quinx_vs_pipiens_sodium_channel_prot.fasta
```

# Annotation statistics

## Generate annotation stats with AGAT

```bash
agat_sp_statistics.pl --gff consensus_simple_lifted_annotation.gff3 -o liftoff_simple/liftoff_simple_annotation_gff_statistics.txt
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


## Calculate the number of predicted protein sequences

```bash
grep -c "^>" liftoff_simple/consensus_simple_CDS.fasta
```

## Calculate the mean protein length (aa)

```bash
seqkit fx2tab -nl consensus_simple_proteins.fasta | \
awk '{sum+=length($2); n++} END {print sum/n}'
```

## Calculate the genome assembly size

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

## Get the K-mer counts with Jellyfish

Genome characteristics were estimated from trimmed Illumina reads using a k-mer based approach.  
21-mer frequencies were generated using Jellyfish with canonical k-mers (`-C`). The k-mer size was selected as **k=21**, which is commonly used for short-read genome size estimation.

The hash size (`-s`) was set to **500M**.

```bash
mkdir jellyfish_output
jellyfish count -C -m 21 -s 500M -t 20 -o jellyfish_output/F7_k21_trimmed.jf <(zcat R1_paired.fastq.gz) <(zcat R2_paired.fastq.gz) > jellyfish_output/jellyfish.log 2>&1 
```

## Generate a k-mer histogram with Jellyfish

```bash
jellyfish histo -t 20 jellyfish_output/F7_k21_trimmed.jf > jellyfish_output/F7_k21_trimmed.histo
```

## Use GenomeScope 2.0 to complete a genome survey analysis (online web version)

The analysis was performed using a diploid model (`p=2`) with the following parameters:

| Parameter | Value | Description |
|-----------|-------|-------------|
| k-mer size (`k`) | 21 | Length of k-mer used for genome profiling |
| Ploidy (`p`) | 2 | Diploid genome model |
| Maximum k-mer coverage (`m`) | 1,000,000 | Maximum coverage threshold included in model fitting |

GenomeScope was run using the Jellyfish-generated histogram as input.

