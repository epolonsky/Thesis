# Bioinformatics Workflow

This file shows the computational workflow used to process sequencing data, generate a consensus genome assembly, assess assembly quality, and annotate the genome using MAKER.

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
| GATK          | 4.6.2.0  |
| BUSCO         | 6.1.0    |
| BlobToolKit   | 4.4.5    |
| MAKER         | 3.01.04  |
| AUGUSTUS      | 3.5.0    |
| RepeatMasker  | 4.2.3    |
| RepeatModeler | 2.0.8    |

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
fastqc F7_PP_DNA_S15_R1_001.fastq.gz F7_PP_DNA_S15_R2_001.fastq.gz
```

---

## 2. Trim sequencing reads

Trim the first 20 bp from each read to remove barcode sequence.

Java memory was increased using the `-Xms` and `-Xmx` options.

```bash
trimmomatic PE -threads 2 -Xms1024m -Xmx8g -phred33 F7_PP_DNA_S15_R1_001.fastq.gz F7_PP_DNA_S15_R2_001.fastq.gz R1_paired.fastq.gz R1_unpaired.fastq.gz R2_paired.fastq.gz R2_unpaired.fastq.gz HEADCROP:20
```

---

## 3. Quality assessment after trimming

```bash
fastqc R1_paired.fastq.gz R2_paired.fastq.gz
```

---

## 4. Map reads to the Culex pipiens reference genome

Reads were mapped against the NCBI chromosome-level assembly:

*Culex pipiens* idCulPipi1.1  
Accession: GCA_963924435.1

Build a Bowtie2 index and align paired-end reads to the reference genome.

```bash
bowtie2-build -f GCA_963924435.1_idCulPipi1.1_genomic.fna ref_genome_index
```

```bash
bowtie2 -x ref_genome_index -1 R1_paired.fastq.gz -2 R2_paired.fastq.gz --un-conc unmapped_reads.fastq.gz -S aligned.sam -p 4 > bowtie2.log 2>&1
```

---

## 5. Extract unmapped reads

Convert unmapped reads into FASTA format for BLAST analysis.

```bash
samtools view -h -f 4 aligned.sam > unaligned.sam

samtools fastq unaligned.sam > unaligned.fastq

seqtk seq -a unaligned.fastq > unaligned.fasta
```

---

## 6. Identify unmapped sequences

BLAST unmapped reads against the NCBI nucleotide database.

```bash
blastn -db nt -query unaligned.fasta -out unaligned_blast.out
```

---

## 7. Convert aligned SAM files to sorted and indexed BAM

```bash
samtools view -bS aligned.sam | samtools sort -o aligned_sorted.bam

samtools index aligned_sorted.bam
```

---

## 8. Calculate sequencing coverage

```bash
samtools coverage aligned_sorted.bam > aligned_sorted_coverage_output.txt
```

---

## 9. Generate a reference-guided consensus genome 

Generate a reference-guided consensus sequence from reads aligned to the reference genome.

```bash
samtools consensus -f fasta aligned_sorted.bam -o consensus.fasta
```

---

## 10. Assess assembly completeness

Evaluate genome completeness using BUSCO.

```bash
busco -i consensus.fasta -l diptera_odb10 -o busco_out -m genome 
```

---

## 11.  Generate BlobToolKit assembly quality visualization snailplot using BUSCO results

```bash
blobtools create --fasta consensus.fasta my_snail_dataset

blobtools add --busco full_table.tsv my_snail_dataset

blobtools view --view snail --plot --out ./ snail_plot
```

---

# Genome Annotation with MAKER

## 12. Install annotation software

GeneMark-ES was installed manually because it is not available through Conda.

Other annotation software:

```bash
conda install -c bioconda repeatmodeler repeatmasker maker augustus snap busco
```

---

## 13. Remove short scaffolds

Scaffolds shorter than 5 kb were removed because GeneMark-ES did not successfully process very small scaffolds.

```bash
seqkit seq -m 5000 consensus.fasta > consensus_min5000.fasta
```

---

## 14. Download annotation evidence

Protein and RNA evidence from *Culex pipiens pallens* were downloaded because annotations were not available for the *Culex pipiens pipiens* assembly.

```bash
datasets download genome accession GCF_016801865.2 --include gff3,rna,protein
```

---

## 15. Configure MAKER for initial annotation

Before running MAKER, the `maker_opts.ctl` configuration file was edited to specify the input genome assembly, evidence files, and initial gene prediction settings.

The reference-guided consensus genome was used as the input assembly:

```bash
genome=/home/elena/data/consensus_min5000.fasta
```

RNA and protein evidence from the related Culex pipiens pallens assembly were provided as homology evidence:

```bash
protein=/home/elena/data/GCF_016801865.2_pallens_protein.faa
est=/home/elena/data/GCF_016801865.2_pallens_rna.fna
```

The MAKER configuration was also set to use soft-masking during annotation:

```bash
softmask=1
```

The initial MAKER configuration specified the GeneMark-ES executable for the initial MAKER run and GeneMark-ES model generation:

```bash
gmhmm=/home/elena/software/gmes_linux_64_4/gmhmme3
```

---

## 16. Train GeneMark-ES

GeneMark-ES was trained independently on the filtered consensus assembly to generate a species-specific hidden Markov model for gene prediction.

```bash
gmes_petap.pl --sequence consensus_min5000.fasta --ES
```

The trained model generated by GeneMark-ES (`output/gmhmm.mod`) was subsequently specified in `maker_opts.ctl`.

---

## 17. Update MAKER with trained GeneMark model

After GeneMark-ES training completed, the `maker_opts.ctl` file was updated to use the newly generated species-specific model:

```bash
gmhmm=/home/elena/data/output/gmhmm.mod
```

---

## 18. Run MAKER with the trained GeneMark model

MAKER was rerun using the GeneMark-ES trained HMM model (`gmhmm.mod`) to generate gene predictions incorporating RNA evidence, protein homology evidence, and ab initio gene prediction.

```bash
maker -base consensus_min5000 maker_opts.ctl maker_bopts.ctl maker_exe.ctl
```

---

## 19. MAKER output files

The final MAKER annotation produced:

- `consensus_min5000.all.gff` — final genome annotation file
- `consensus_min5000.all.maker.proteins.fasta` — predicted protein sequences
- `consensus_min5000.all.maker.transcripts.fasta` — predicted transcript sequences

These files were used for downstream annotation analysis.

---

## 20. Train additional ab initio predictors (future refinement)

Future MAKER iterations will incorporate trained SNAP and AUGUSTUS models generated from the initial MAKER annotation. These additional predictors will be used alongside GeneMark-ES to improve gene model prediction.

---

# GATK Pileup

## 21. Add read groups

```bash
gatk AddOrReplaceReadGroups -I aligned_sorted.bam -O aligned_sorted_RG.bam --RGID F7_PP_DNA_S15 --RGLB lib1 --RGPL ILLUMINA --RGPU S15 --RGSM F7_PP_DNA --CREATE_INDEX true
```

---

## 22. Index the reference genome

```bash
samtools faidx GCA_963924435.1_idCulPipi1.1_genomic.fna
```

---

## 23. Create the sequence dictionary

```bash
gatk CreateSequenceDictionary -R GCA_963924435.1_idCulPipi1.1_genomic.fna
```

---

## 24. Generate pileup

Generate a pileup file describing the aligned reads at each genomic position.

```bash
gatk Pileup -R GCA_963924435.1_idCulPipi1.1_genomic.fna -I aligned_sorted_RG.bam -O pileup_out.txt
```
