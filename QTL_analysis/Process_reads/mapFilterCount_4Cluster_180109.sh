#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=4,mem=32gb
#PBS -o /LOG_LOC/logs/
#PBS -e /LOG_LOC/logs/


# load software
# note that this will only be available at runtime; it will not be available after the script finished
bwaVersion='bwa/0.7.15'
module load ${bwaVersion}
samtoolsVersion='samtools/1.5' #was 1.4
module load ${samtoolsVersion}
module load bamtools

# careful here: if run on MSI, need to request enough threads if want to speed up
# might be easier to just run single thread, since this happens on highly multiplexed samples
#nThreads=24

#fore testing:
#readDat='/LOCATION/'
#BWADat='/LOCATION/'
#genome='/LOCATION/sacCer3.fa'
#SNPs='/LOCATION/SNPs_Maggie_170809_BY_positions.txt'
#fileRoot='DPFA0001_A01_S1_'
#fileType='fastq'
#pairID1="R1"
#pairID2="R2"
#postfix="_001"
#PE='PE'

# the definitions above were for setting up. Now do this right, with getopt:
# see
# http://frodo.looijaard.name/project/getopt/misc/examples
# http://www.bahmanm.com/blogs/command-line-options-how-to-parse-in-bash-using-getopt


fileR1=${readDat}${fileRoot}${pairID1}${postfix}.${fileType}
fileR2=${readDat}${fileRoot}${pairID2}${postfix}.${fileType}

if [ ! -e ${fileR1} ]
then
    echo "ERROR: File ${fileR1} not found. Exiting"; exit 1
fi

if [ "${PE}" = "PE" ]
then
    echo "Runing in PE mode"
    if [ ! -e ${fileR2} ]
    then
        echo "ERROR: File ${fileR2} not found in PE mode. Exiting"; exit 1
    fi
    bwa mem -t ${nThreads} ${genome} ${fileR1} ${fileR2} | samtools sort -@${nThreads} -O BAM -o ${BWADat}"/"${fileRoot}sort.bam -
else
if [ "${PE}" = "SE" ]
then
    echo "Running in SE mode"
    bwa mem -t ${nThreads} ${genome} ${fileR1} | samtools sort -@${nThreads} -O BAM -o ${BWADat}"/"${fileRoot}sort.bam -

else
    echo "ERROR: PE/SE mode not set correctly. Don't know what to do; Exiting"; exit 1
fi
fi

# the above will have either produced a sam or exited with a message

samtools view -q 30 ${BWADat}"/"${fileRoot}sort.bam | grep 'XS:i:0' | samtools view -b -T ${genome} - > ${BWADat}"/"${fileRoot}sort_filtered.bam #get ride of the missmatch
samtools rmdup -S ${BWADat}"/"${fileRoot}sort_filtered.bam ${BWADat}"/"${fileRoot}sort_filtered_rmdup.bam #get ride of the pcr duplicate

samtools mpileup -vu -t INFO/AD -l ${SNPs} -f ${genome} ${BWADat}"/"${fileRoot}sort_filtered_rmdup.bam > ${BWADat}"/"${fileRoot}sort_filtered_rmdup.vcf #counting the covering per SNPs

echo "done"
