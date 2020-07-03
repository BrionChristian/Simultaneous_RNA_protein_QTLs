#!/bin/sh

# once this finishes, make sure all files succeeded
# e.g. by counting the number of vcfs in the alignments folder: ls -alh *vcf | wc -l
# might need to rerun a few by hand, by switching to the non-cluster script and specifycing individual samples
# can all be done by commenting in/out lines below

# these should usually not need to change:
mapScript='/SCRIPT_LOCATION/mapFilterCount_4Cluster_180109.sh'


genome='/GENOME_LOCATION/sacCer3.fa'
# be sure the SNP set is two columns for the reference genome; and that it has the correct line breaks (safe in xcode once, if necessary)
SNPs='/SNP_LIST_LOCATION/SNPs_Maggie_170809_BY_positions.txt'
nThreads=24

# these are probably fine, but double check:
fileType='fastq'
pairID1='R1'
pairID2='R2'
postfix='_001'

# NEED TO adjust these as needed:
readDat='/READ_LOCATION/'
BWADat='/OUTPUR_LOCATION/'
fileType='fastq'
PE='PE'

# need to be careful that the underscores in front of ${pairID1} match
fileRoots=( $(ls ${readDat} | grep "${pairID1}" | sed "s/${pairID1}.*//") )

for thisRoot in ${fileRoots[@]}; do
#for thisRoot in ${fileRoots[@]:1:3}; do
    echo ${thisRoot}
    qsub -v fileRoot=${thisRoot},PE=${PE},nThreads=${nThreads},readDat=${readDat},BWADat=${BWADat},genome=${genome},SNPs=${SNPs},pairID1=${pairID1},pairID2=${pairID2},postfix=${postfix},fileType=${fileType} ${mapScript}
done

#sh ${mapScript} -f ${fileRoots[0]} -p ${PE} -t ${nThreads} --readsFolder ${readDat} --outputFolder ${BWADat} --genome ${genome} --SNPs ${SNPs} --r1 ${pairID1} --r2 ${pairID2} --postFix ${postfix} --fileType ${fileType}
#sh ${mapScript} -f "DPFA0001_C10_S30_" -p ${PE} -t ${nThreads} --readsFolder ${readDat} --outputFolder ${BWADat} --genome ${genome} --SNPs ${SNPs} --r1 ${pairID1} --r2 ${pairID2} --postFix ${postfix} --fileType ${fileType}
