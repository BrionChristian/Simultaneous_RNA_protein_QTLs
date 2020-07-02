#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=4,mem=32gb
#PBS -o /LOG_LOC/logs/
#PBS -e /LOG_LOC/logs/
# using info from bioanalyzer chip on average library size; guessing at SD
# need to run one file at a time

#inputFiles="keyRNA.txt" # provided in this code directory

inputFolder="INPUT_FOLDER_WITH_FASTQ_READS"
outputFolder1="OUTPUT_FOLDER/trimmomatic/"
outputFolder2="OUTPUT_FOLDER/reads/"


echo "sample TRIM: ${lineArray}"

java -jar TRIM_LOCATION/trimmomatic-0.39.jar SE -threads 14 \
    ${inputFolder}${lineArray}*.fastq.gz \
    ${outputFolder1}${lineArray}.fastq.gz \
	ILLUMINACLIP:/TRIM_ADAPTERS_LOC/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36







