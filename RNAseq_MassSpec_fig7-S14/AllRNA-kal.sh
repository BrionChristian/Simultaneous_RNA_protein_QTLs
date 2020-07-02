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


echo "sample kal: ${lineArray}"

	# using info from bioanalyzer chip on average library size; guessing at SD
	# run one file at a time

	#ls ${outputFolder1}${lineArray[0]}.fastq.gz

	# note that while bioanalyzer library has mean length of 336, this includes adapters – go with the standard 200 instead
	# genomebam: the chromosomes file had a line ending issue!

mkdir ${outputFolder2}${lineArray}

~/src/kallisto/kallisto quant \
    -i /INDEX_LOCATION/Saccharomyces_cerevisiae.R64-1-1.cdna.all.kallistoIndex \
    -t 24 --single -l 200 -s 30 --rf-stranded -b 100 \
    -o ${outputFolder2}${lineArray} \
    --genomebam --gtf /MOD_LOCATION/Saccharomyces_cerevisiae.R64-1-1.93_mod.gtf.gz \
    --chromosomes /LOC_LIST_CHROMOSOMES/chrLengths_ensemblFormat.txt \
    ${outputFolder1}${lineArray}.fastq.gz


# worked well; worst were 2C & 1A with 5.9% of reads dropped





