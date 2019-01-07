#!/bin/bash

# workflow and command lines used to process the dataset for CID4404 
# raw bcl files from the NextSeq500 downloaded from Illumina basespace
# the following processing was performed using cellranger version 2.1.1

# sequences were demultiplexed using the mkfastq command
qsub -cwd -b y -j y -pe smp 32 -l h_vmem=300G -P TumourProgression \
	-N mkfastq -V cellranger mkfastq --run=./10X_clinical+PDX-24082017_bclfiles/Files/ \
	--samplesheet=./input_samplesheet.csv --lanes=1,2,3,4 --delete-undetermined

# alignment, filtering, barcode counting and UMI counting was performed using the count command
qsub -cwd -b y -j y -pe smp 32 -l h_vmem=300G -P TumourProgression \
	-N mkfastqCS2 -V cellranger count --id=count_4404_primary_GRCh38 \
	--sample=4404_primary --transcriptome=./refdata-cellranger-GRCh38-1.2.0/ \
	--fastqs=./10X_clinical+PDX-24082017/H5M3VBGX3_1234/outs/fastq_path/Clinical_10X/4404_primary/ \
	--indices=SI-GA-B2