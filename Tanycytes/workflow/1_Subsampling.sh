#!/bin/bash
DONE=(SRR28910565)

SRR_FILES=(SRR28910565 SRR28910566 SRR28910567 SRR28910568)
THREADS=32

mkdir -p ./data/Raw
cd ./data/Raw

prefetch ${SRR_FILES[@]}

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

  seqtk sample -s100 Fastq/${SRR}_1.fastq.gz 0.1 | gzip > "../Fastq_subsampled/${SRR}_1_sub.fastq.gz"
  seqtk sample -s100 Fastq/${SRR}_2.fastq.gz 0.1 | gzip > "../Fastq_subsampled/${SRR}_2_sub.fastq.gz"
  seqtk sample -s100 Fastq/${SRR}_3.fastq.gz 0.1 | gzip > "../Fastq_subsampled/${SRR}_3_sub.fastq.gz"
done

cd Fastq_subsampled
for prefix in "${SRR_FILES[@]}"; do
  echo "Renaming files for $prefix"
  mv "${prefix}_1_sub.fastq.gz" "${prefix}_S1_L001_I1_001.fastq.gz"  
  mv "${prefix}_2_sub.fastq.gz" "${prefix}_S1_L001_R1_001.fastq.gz"
  mv "${prefix}_3_sub.fastq.gz" "${prefix}_S1_L001_R2_001.fastq.gz"
done


wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa 10 > ../Reference_ch10/chr10.fa

wget https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
gunzip Mus_musculus.GRCm39.112.gtf.gz
awk '$1 == "10" || $1 ~ /^#/' Mus_musculus.GRCm39.112.gtf > ../Reference_ch10/chr10.gtf


