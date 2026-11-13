# Bioinformatics Workflow

This file shows the computational workflow used to process illumina short read sequencing data, generate a consensus genome assembly, assess assembly quality, annotate the existing *Culex pipiens pipiens* genome using BRAKER4 and eggNOG, and transfering that annotation to the consensus genome assembly.

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

## Assess trimmed read quality with FastQC

```bash
fastqc -o fastqc_output trimmomatic_output/R1_paired.fastq.gz trimmomatic_output/R2_paired.fastq.gz
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
samtools view -bS bowtie2_output/aligned.sam -o bowtie2_output/aligned.bam

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
conda activate busco
mkdir busco_simple_out
busco -i consensus_simple.fasta -l diptera_odb10 -o busco_simple_out -m genome 
```

##  Generate snailplot from BUSCO results with BlobToolKit

```bash
cd busco_simple_out

conda activate btk

blobtools create --fasta ~/data/consensus_simple.fasta my_snail_dataset

blobtools add --busco run_diptera_odb10/full_table.tsv my_snail_dataset

blobtools view --view snail --plot --out ./ my_snail_dataset
```

## Install and run repeatmodeler 2

```bash
cd ~/data
conda create -n repeatmod_env -c conda-forge -c bioconda repeatmodeler=2.0.9
conda activate repeatmod_env

BuildDatabase -name idCulPipi1_genome_db GCA_963924435.1_idCulPipi1.1_genomic.fna
nohup RepeatModeler -database idCulPipi1_genome_db -threads 8 -dir repeat_modled_idCulPipi1 & 
```

## Install and run repeatmasker

```bash
cd ~/data
conda create -n repeatmask_env -c conda-forge -c bioconda repeatmasker=4.2.4
conda activate repeatmask_env

nohup RepeatMasker -pa 4 -gff -xsmall -lib repeat_modled_idCulPipi1/consensi.fa -dir repeat_masked_idCulPipi1 GCA_963924435.1_idCulPipi1.1_genomic.fna &
```

## Install and run Braker for gene prediction

```bash
conda create -n braker4_env -c conda-forge python=3.11 apptainer snakemake=8.18.2 pandas
conda activate braker4_env
git clone https://github.com/Gaius-Augustus/BRAKER4
cd ~/data/BRAKER4
cp config.ini.example config.ini
```

Change what is in the samples.csv file using nano to:

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
idCulPipi1,/home/elena/data/GCA_963924435.1_idCulPipi1.1_genomic.fna,/home/elena/data/repeat_masked_idCulPipi1/GCA_963924435.1_idCulPipi1.1_genomic.fna.masked,/home/elena/data/quinx/GCF_015732765.1_quinx_protein.faa:/home/elena/data/pallens/GCF_016801865.2_pallens_protein.faa,,,,,,,,,diptera_odb12,
```

Install the busco lineage

```bash
conda activate busco
cd ~/data/BRAKER4
busco --download diptera_odb12 --download_path shared_data/busco_downloads
```

Run Braker through snakemake

```bash
conda activate braker4_env
cd ~/data/BRAKER4
nohup snakemake --cores 8 --use-singularity --singularity-prefix .singularity_cache --singularity-args "-B /home" --latency-wait 120 --restart-times 3 &
```

## Install and run eggNOG mapper for functional annotation

```bash
cd ~/data
git clone https://github.com/eggnogdb/eggnog-mapper
conda create -n eggnog_env -c conda-forge -c bioconda python=3.11
conda activate eggnog_env
pip install eggnog-mapper
cd ~/data/eggnog-mapper
nohup python download_eggnog_data.py -y --data_dir ~/data/eggnog-mapper/data &

nohup python emapper.py -i ~/data/BRAKER4/output/idCulPipi1/results/braker.aa.gz -o idCulPipi1_eggnog --data_dir ~/data/eggnog-mapper/data --cpu 8 &
```

## Combine the gene prediction and functional annotation outputs into one annotation file with eggNOG decorate gff

```bash
conda activate eggnog_env
cd ~/data/eggnog-mapper

gzip -dc ~/data/BRAKER4/output/idCulPipi1/results/braker.gff3.gz > ~/data/BRAKER4/output/idCulPipi1/results/braker.gff3
gzip -dc ~/data/BRAKER4/output/idCulPipi1/results/braker.aa.gz > ~/data/BRAKER4/output/idCulPipi1/results/braker.aa

sed -i 's/all_orthologies, annot_orthologs, _confidence) = annotation/all_orthologies, annot_orthologs, _confidence, _tax_ceiling, _donor_taxid, _donor_lineage) = annotation/' ~/data/eggnog-mapper/eggnogmapper/deco/decoration.py

nohup python emapper.py -i ~/data/BRAKER4/output/idCulPipi1/results/braker.aa --itype proteins -m no_search --annotate_hits_table idCulPipi1_eggnog.emapper.seed_orthologs --decorate_gff ~/data/BRAKER4/output/idCulPipi1/results/braker.gff3 -o idCulPipi1_decorated --data_dir ~/data/eggnog-mapper/data --cpu 8 &
```

## Get annotation statistics with AGAT

```bash
conda activate agat_env

# idCulPipi1 renefrence BRAKER4 annotation statistics
agat_sp_statistics.pl --gff ~/data/BRAKER4/output/idCulPipi1/results/braker.gff3 --gs ~/data/consensus_simple.fasta -o ~/data/BRAKER4/output/idCulPipi1/results/annotation_statistics.txt

seqkit stats ~/data/BRAKER4/output/idCulPipi1/results/braker.aa

# liftoff consensus annotation statistics
agat_sp_statistics.pl --gff ~/data/liftoff_annotation/consensus_simple_idCulPipi1_lifted_annotation.gff3 --gs ~/data/consensus_simple.fasta -o ~/data/liftoff_annotation/annotation_statistics.txt

seqkit stats ~/data/liftoff_annotation/consensus_simple_proteins_idCulPipi1.faa
```

## Transfer the annotation to my consensus sequence with Liftoff 

```bash
cd ~/data
mkdir liftoff_annotation
liftoff -g eggnog-mapper/idCulPipi1_decorated.emapper.decorated.gff -o liftoff_annotation/consensus_simple_idCulPipi1_lifted_annotation.gff3 consensus_simple.fasta GCA_963924435.1_idCulPipi1.1_genomic.fna 
```

## BLAST reference insectiside resistant gene protein files against BRAKER4 proteins fasta file

```bash
cd ~/data
mkdir blast_ins_res_genes
cd blast_ins_res_genes

cat ~/data/quinx_pallens_genes/quinx_pallens_prot/quinx_sodium_channel_prot.fasta ~/data/quinx_pallens_genes/quinx_pallens_prot/quinx_ace1_prot.fasta ~/data/quinx_pallens_genes/quinx_pallens_prot/quinx_rdl_gaba_prot.fasta ~/data/quinx_pallens_genes/quinx_pallens_prot/pallens_sodium_channel_prot.fasta ~/data/quinx_pallens_genes/quinx_pallens_prot/pallens_ace1_prot.fasta ~/data/quinx_pallens_genes/quinx_pallens_prot/pallens_rdl_gaba_prot.fasta > ~/data/blast_ins_res_genes/reference_proteins.faa

makeblastdb -in ~/data/BRAKER4/output/idCulPipi1/results/braker.aa -dbtype prot -out braker_proteins_db

mkdir output

blastp -query reference_proteins.faa -db braker_proteins_db  -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qstart qend sstart send' -out output/target_gene_hits.tsv
```

# Insecticide resistant gene extraction

The gene sequences along with their CDS and protein sequences for *Culex pipiens pallens* and *Culex quinquefasciatus* were downloaded from NCBI using accession numbers found through literature search.

For *Culex pipiens pipiens* the BRAKER4 protein id was used to extract the protein sequences from the liftoff annotation.

```bash
cd ~/data
mkdir pipiens_ins_res_genes_idCulPipi1_lifted
mkdir pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_prot
mkdir pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_cds

gffread liftoff_annotation/consensus_simple_idCulPipi1_lifted_annotation.gff3 -g consensus_simple.fasta -y liftoff_annotation/consensus_simple_proteins_idCulPipi1.faa -x liftoff_annotation/consensus_simple_CDS_idCulPipi1.fasta
```
```bash
# VGSC: g13896.t4

seqkit grep -p 'g13896.t4' liftoff_annotation/consensus_simple_proteins_idCulPipi1.faa > pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_prot/pipiens_sodium_channel_prot.fasta

# RDL: g22203.t1 
seqkit grep -p 'g22203.t1' liftoff_annotation/consensus_simple_proteins_idCulPipi1.faa >  pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_prot/pipiens_rdl_gaba_prot.fasta

# ACE1: g17576.t1
seqkit grep -p 'g17576.t1' liftoff_annotation/consensus_simple_proteins_idCulPipi1.faa >  pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_prot/pipiens_ace1_prot.fasta
```

```bash
# VGSC CDS: g13896.t4
seqkit grep -p 'g13896.t4' liftoff_annotation/consensus_simple_CDS_idCulPipi1.fasta > pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_cds/pipiens_sodium_channel_cds.fasta

# RDL CDS: g22203.t1 
seqkit grep -p 'g22203.t1'  liftoff_annotation/consensus_simple_CDS_idCulPipi1.fasta > pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_cds/pipiens_rdl_gaba_cds.fasta

# ACE1 CDS: g17576.t1
seqkit grep -p 'g17576.t1' liftoff_annotation/consensus_simple_CDS_idCulPipi1.fasta >  pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_cds/pipiens_ace1_cds.fasta
```

## BLAST reference rdl cds against consensus genome

```bash
cd ~/data
mkdir blast_rdl_to_consensus_simple
cd blast_rdl_to_consensus_simple

cat ~/data/quinx_pallens_genes/quinx_pallens_cds/quinx_rdl_gaba_cds.fasta ~/data/quinx_pallens_genes/quinx_pallens_cds/pallens_rdl_gaba_cds.fasta > ~/data/blast_rdl_to_consensus_simple/reference_rdl_cds.fasta

makeblastdb -in ~/data/consensus_simple.fasta -dbtype nucl -out consensus_simple_db

mkdir output

blastn -query reference_rdl_cds.fasta -db consensus_simple_db -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qstart qend sstart send' -out output/rdl_gene_hits.tsv

# Sort RDL hits by reference CDS position to identify matching fragments
awk '$1 ~ /XM_038263969/' output/rdl_gene_hits.tsv | sort -k9,9n
```

## Extract the sequence for rdl from the BLAST results

```bash
samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172174780-172303616' > output/rdl_genomic_region.fasta
```

## Find where the rdl gene is in the consensus sequence with Exonerate

```bash
cd ~/data/blast_rdl_to_consensus_simple

exonerate --model coding2genome --showalignment yes --showvulgar yes --showtargetgff yes --showsugar no reference_rdl_cds.fasta output/rdl_genomic_region.fasta > output/rdl_exonerate_results.txt
```

```bash
# Compare reference alignments and scores
grep -E 'Query:|Raw score:|Query range:|Target range:' output/rdl_exonerate_results.txt

# View positive-strand exon predictions
grep -E $'\texon\t' output/rdl_exonerate_results.txt | awk -F '\t' '$7 == "+"'

# Convert Exonerate exon coordinates to chromosome coordinates
grep -E $'\texon\t' output/rdl_exonerate_results.txt | awk -F '\t' '$7 == "+" {print "OZ004313.1:" 172174780+$4-1 "-" 172174780+$5-1}' | sort -u -t: -k2,2n
```

```bash
{
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172174780-172174904'
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172226606-172226703'
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172233935-172234002'
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172255260-172255486'
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172257954-172258036'
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172278616-172278753'
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172295112-172295396'
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172298418-172298908'
  samtools faidx ~/data/consensus_simple.fasta 'OZ004313.1:172303452-172303616'
} | awk '
BEGIN {print ">RDL_consensus_CDS"}
!/^>/ {printf "%s", $0}
END {printf "\n"}
' > output/rdl_consensus_cds.fasta
```

```bash
# Translate the cds to protein with seqkit
seqkit translate output/rdl_consensus_cds.fasta > output/rdl_consensus_prot.fasta
```

## Download reference sequences from NCBI with known mutation and make files for alignment

```bash
mkdir ref_mutation
mkdir ref_mutation/files_for_aligment
cd ~/data

cat ref_mutation/quinx_ace1_mut_ref.fasta pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_prot/pipiens_ace1_prot.fasta > ref_mutation/files_for_aligment/quinx_pipiens_ace1_alignment_input.fasta

cat ref_mutation/quinx_rdl_mut_ref.fasta pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_prot/pipiens_rdl_gaba_prot.fasta > ref_mutation/files_for_aligment/quinx_pipiens_rdl_gaba_alignment_input.fasta

cat ref_mutation/quinx_sodium_channel_mut_ref.fasta pipiens_ins_res_genes_idCulPipi1_lifted/pipiens_prot/pipiens_sodium_channel_prot.fasta > ref_mutation/files_for_aligment/quinx_pipiens_sodium_channel_alignment_input.fasta

cat ref_mutation/quinx_rdl_mut_nucleotide.fasta blast_rdl_to_consensus_simple/output/rdl_consensus_cds.fasta >  ref_mutation/files_for_aligment/quinx_pipiens_consensus_nucleotide_rdl_gaba_alignment_input.fasta
```

-
-
-
-
-
-

---

---

# Code Draft Section

## Annotation statistics

## Generate annotation stats with AGAT

```bash
agat_sp_statistics.pl --gff consensus_simple_lifted_annotation.gff3 -o liftoff_simple/liftoff_simple_annotation_gff_statistics.txt
```

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

