# Ebola-VP40 
## Differential expression analysis of cells infected with EBOV, RESTV, or mocks infected with culture medium

#### Step 1: Download the FASTQ files from NCBI BioProject PRJNA1040271
#### Step 2: Check the quality of the FASTQ files
#### Step 3: Human genome data download
    wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
    gunzip Homo_sapiens.GRCh38.109.gtf.gz
#### Step 4: Building human genome index
    mkdir HumanGenome_STARIndex
    STAR --runThreadN 20 --runMode genomeGenerate --genomeDir HumanGenome_STARIndex --genomeFastaFiles HumanGenome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.109.gtf --sjdbOverhang 150
#### Step 5: Aligning FASTQ reads to the human genome index
    python StarAlign.py Input.txt HumanGenome_STARIndex FASTQ Alignment

    Input.txt= Tab delimited file listing paired-end FASTQ files
    HumanGenome_STARIndex = Folder with human genome index
    FASTQ = Folder with FASTQ files in .gz format
    Alignment = Folder to save the output files
#### Step 6: To generate read counts
    Rscript ReadCount.R
