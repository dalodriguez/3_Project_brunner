#  README

## 1. Data sources

This task is based on publicly available single cell RNA sequencing data from the following study:

Title: Fasting induces metabolic switches and spatial redistributions of lipid processing and neuronal interactions in tanycytes
DOI: 10.1038/s41467-024-50913-w

The study evaluate the role of tanycytes, ependymal cells lining the third ventricle (3V) in the mediobasal hypothalamus of the brain, in the control of metabolism. The workflow in this repository replicate the main steps of the data analysis described in the paper, focused on feeding (wild type) samples (SRR28910565 SRR28910566 SRR28910567 SRR28910568), originally sequenced with Illumina HiSeq 4000.

The subsampled FASTQs are stored in `data/Fastq_subsampled` and are used as the inputs for the workflow.

The code to generate the response to the questions are in `workflow/3_Answer_Questions.sh`

## 2. How to download

The following samples were retrieved from SRA

```bash
SRR_FILES=(SRR28910565 SRR28910566 SRR28910567 SRR28910568)
prefetch ${SRR_FILES[@]}
```


## 3. Pre-processing / subsampling

Sample were miniaturized using the code in 1_Subsampling, not needed for Taiga submission. 

###  [1] – Fastq dump
First, raw .SRA samples were converted to fastq and stored in `data/Fastq`

```bash
for SRR in "${SRR_FILES[@]}"; do
  echo "Processing $SRR"
  fasterq-dump Raw/${SRR}/${SRR}.sra \
    --split-files \
    --include-technical \
    -O Fastq \
    --threads $THREADS
  pigz -p $THREADS Fastq/${SRR}_1.fastq
  pigz -p $THREADS Fastq/${SRR}_2.fastq
  pigz -p $THREADS Fastq/${SRR}_3.fastq
done
```
###  [2] – Fastq subsampling
Then, Fastq samples were miniaturized using seqtk to maintain only the 10% of reads. Subsamples are stored in `data/Fastq_subsampled`

```bash
for SRR in "${SRR_FILES[@]}"; do
  seqtk sample -s100 Fastq/${SRR}_1.fastq.gz 0.1 | gzip > "../Fastq_subsampled/${SRR}_1_sub.fastq.gz"
  seqtk sample -s100 Fastq/${SRR}_2.fastq.gz 0.1 | gzip > "../Fastq_subsampled/${SRR}_2_sub.fastq.gz"
  seqtk sample -s100 Fastq/${SRR}_3.fastq.gz 0.1 | gzip > "../Fastq_subsampled/${SRR}_3_sub.fastq.gz"
```
###  [3] – Reference genome subsampling

The mouse GRCm39 reference genome was download and subsampled to contain only the sequence of chromosome 10. The subsampled genome is stored in `data/Reference_ch10`.

```bash
wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa 10 > ../Reference_ch10/chr10.fa

wget https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
gunzip Mus_musculus.GRCm39.112.gtf.gz
awk '$1 == "10" || $1 ~ /^#/' Mus_musculus.GRCm39.112.gtf > ../Reference_ch10/chr10.gtf
```

## 4. How the workflow works

The workflow files are stored in `workflow/`. The files needed are `workflow/2_Workflow.sh` and `workflow/3_Workflow_Seurat.R.sh`


### Step 1 - Building subset reference genome (chromosome 10) 

**Purpose:** To build the reference genome 
**Tools:** `cellranger`
**Inputs:** Need sequence .fa and gene .gft data (from `data/Reference_ch10/`)
**Outputs:** Reference genome folder (in `data/GRCm39`)
**Command:**

```bash
cellranger mkref \
  --genome=GRCm39 \
  --fasta=Reference_ch10/chr10.fa \
  --genes=Reference_ch10/chr10.gtf
```

---
### Step 2 - Cellranger pipeline

**Purpose:** Perform alignment and generate count matrix
**Tools:** `cellranger`
**Inputs:** Need Fastq files and reference genome folder (from `data/Fastq_sumsampled/` and `data/GRCm39`)
**Outputs:** A folder per sample (SRR28910565 SRR28910566 SRR28910567 SRR28910568)
**Command:**

```bash
for SRR in "${SRR_FILES[@]}"; do
    cellranger count \
    --id=$SRR \
    --transcriptome=GRCm39 \
    --fastqs=Fastq_subsampled \
    --sample=$SRR \
    --create-bam=false
done
```

Individual matrices are then aggregated in `data/Combined/`

The `data/aggr.csv` is required, containing the list of SRR files and the path to molecule_info.h5 files. 

```bash
cellranger aggr --id=Combined --csv=aggr.csv --normalize=mapped
```


### Step 3 – Downstream Analysis

**Purpose:** Generate a seurat object, perform quality control
**Tools:** `R4.4` `Seurat_5.3.0` `DropletUtils_1.28.1` `scuttle_1.18.0`
**Inputs:** Filtered feature bc matrix from the aggregated cellranger analysis (from `data\Combined\outs\count\filtered_feature_bc_matrix`)
**Outputs:** RDS seurat object
**Command:**

The full R pipeline is in the `workflow/3_Workflow_Seurat.R.sh` file

