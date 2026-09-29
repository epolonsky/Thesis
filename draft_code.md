Swedish Genome Annotation codes

# Install and run repeatmodeler 2

```bash
cd ~/data
conda create -n repeatmod_env -c conda-forge -c bioconda repeatmodeler=2.0.9
conda activate repeatmod_env

BuildDatabase -name idCulPipi1_genome_db GCA_963924435.1_idCulPipi1.1_genomic.fna
nohup RepeatModeler -database idCulPipi1_genome_db -threads 8 &
```

# Install and run repeatmasker

```bash
cd ~/data
conda create -n repeatmask_env -c conda-forge -c bioconda repeatmasker=4.2.4
conda activate repeatmask_env

nohup RepeatMasker -pa 4 -gff -xsmall -lib idCulPipi1_genome_db-families.fa GCA_963924435.1_idCulPipi1.1_genomic.fna &
```
