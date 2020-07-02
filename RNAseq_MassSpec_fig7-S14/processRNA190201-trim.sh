#!/bin/sh

# once this finishes, make sure all files succeeded
# e.g. by counting the number of vcfs in the alignments folder: ls -alh *vcf | wc -l
# might need to rerun a few by hand, by switching to the non-cluster script and specifycing individual samples
# can all be done by commenting in/out lines below

# bash processRNA190201.sh


# these should usually not need to change:
processScript='AllRNA-trim.sh'
inputFiles="keyRNA.txt" # provided in this code directory

# these are probably fine, but double check:


# need to be careful that the underscores in front of ${pairID1} match
while IFS=$'\t' read -r -a lineArray
do
	echo ${lineArray[0]}
    qsub -v lineArray=${lineArray[0]} ${processScript}
done < "$inputFiles"