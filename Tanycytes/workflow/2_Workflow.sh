#!/bin/bash
SRR_FILES=(SRR28910565 SRR28910566 SRR28910567 SRR28910568)

#export PATH=$PATH:data/cellranger-9.0.1
#cellranger --version

# BUILD REFERENCE GRCm39 GENOME FOR CHROMOSOME 10
cellranger mkref \
  --genome=GRCm39 \
  --fasta=Reference_ch10/chr10.fa \
  --genes=Reference_ch10/chr10.gtf

# GENERATE ALIGNMENT & COUNT MATRICES PER SAMPLE
for SRR in "${SRR_FILES[@]}"; do
    cellranger count \
    --id=$SRR \
    --transcriptome=GRCm39 \
    --fastqs=Fastq_subsampled \
    --sample=$SRR \
    --create-bam=false
done

# GENERATE AGGREGATED COUNT MATRIX
cellranger aggr --id=Combined --csv=aggr.csv --normalize=mapped
