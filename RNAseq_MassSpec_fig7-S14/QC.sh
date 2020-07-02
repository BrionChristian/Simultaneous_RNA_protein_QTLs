#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=4,mem=32gb
#PBS -o /LOG_LOC/logs/
#PBS -e /LOG_LOC/logs/

#module load fastqc
module load python2


#fastqc /TRIM_READS_LOC/*.fastq.gz --/QC_LOC/QC 
/PYTHON_SCRIPT/geneBody_coverage.py -i `ls -dm /QC_LOC/QC/pseudoalignments*.bam | tr -d ' \n'` -r /MOD_BED_LOC/ensemblGenes_fromUCSC_mod.bed -o RSeQC_geneBodyCoverage
/PYTHON_SCRIPT/tin.py -i `ls -dm /QC_LOC/QC/pseudoalignments*.bam | tr -d ' \n'` -r /MOD_BED_LOC/ensemblGenes_fromUCSC_mod.bed > RSeQC_TIN