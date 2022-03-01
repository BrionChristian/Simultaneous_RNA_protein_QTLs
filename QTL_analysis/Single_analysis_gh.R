# Author: Christian Brion - 2020 - UMN
#
#-convert allele frequency vcf files into R data
#-prepare meta files with the list of populations
#-provide allele frequency figures for QC
#


# can then call this from external
i = 1
alignmentDir <- "YOUR/WORK/DIRECTORY/"
setwd(alignmentDir)

library("VariantAnnotation")
source(paste(alignmentDir,"/Scripts/x_qtl_seq_functions_170831.R",sep=""))



SNPs <- read.table(paste(alignmentDir,"/Scripts/SNPs_Maggie_170809_BY_positions.txt",sep=""), stringsAsFactors=FALSE, head=FALSE)
# see comments above. As of 8/31/17, the SNPs seem not to be fully filtered, and are out of sorting order
for (thisChr in unique(SNPs[,1])){SNPs[SNPs[,1] == thisChr, 2] <- sort(SNPs[SNPs[,1] == thisChr, 2])}
SNPs <- rbind(SNPs[SNPs[,1] == "chrI",], SNPs[SNPs[,1] == "chrII",], SNPs[SNPs[,1] == "chrIII",], SNPs[SNPs[,1] == "chrIV",], SNPs[SNPs[,1] == "chrV",], SNPs[SNPs[,1] == "chrVI",], SNPs[SNPs[,1] == "chrVII",], SNPs[SNPs[,1] == "chrVIII",], SNPs[SNPs[,1] == "chrIX",], SNPs[SNPs[,1] == "chrX",], SNPs[SNPs[,1] == "chrXI",], SNPs[SNPs[,1] == "chrXII",], SNPs[SNPs[,1] == "chrXIII",], SNPs[SNPs[,1] == "chrXIV",], SNPs[SNPs[,1] == "chrXV",], SNPs[SNPs[,1] == "chrXVI",])

alignmentDir2 <- paste(alignmentDir,"",sep="")


resultsFolder <- paste(alignmentDir2,"AlleleFreqresults/",sep="")
sgd_table <- paste(alignmentDir,"/Scripts/sacCer3ChromLenghts.txt",sep="")


# common annotations, functions, etc
geneInfo = read.table(paste(alignmentDir,"/Scripts/ensemblGenes_ensembl83_160307_MOD.txt",sep=""), stringsAsFactors=FALSE, sep="\t", header=TRUE)
rownames(geneInfo) <- geneInfo[,"geneID"]
allNames <- geneInfo[, "geneName"]
names(allNames) <- geneInfo[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

sepBetweenChr <- 1e5
LoessSpan = 0.1

# done reading annotation files, now get to work
FileList <- dir(paste(alignmentDir2,"alignments",sep=""), pattern="vcf", full.names=TRUE)

finalINFO<-as.data.frame(matrix("",nrow = length(FileList),ncol = 8),stringsAsFactors=F)
colnames(finalINFO)<-c("orf","gene","crispr","rep","tech","gate","seqnum","cover")

for (i in 1:length(FileList)) {
  File<-FileList[i]
  fileNAME<-basename(File)
  fileinfo<-strsplit(fileNAME, split = "[-_.]")[[1]]
  
  if (fileinfo[2]=="NAi") {
    finalINFO$gene[i]<-"NA"
    finalINFO$orf[i]<-"YFL039C"
    finalINFO$crispr[i]<-"i"
  } else {
    finalINFO$gene[i]<-substr(fileinfo[2],1,4)
    finalINFO$orf[i]<-geneInfo$geneID[geneInfo$geneName==finalINFO$gene[i]]
    finalINFO$crispr[i]<-substr(fileinfo[2],5,7)
  }
  finalINFO$rep[i]<-"0"
  finalINFO$tech[i]<-"1"
  finalINFO$gate[i]<-fileinfo[3]
  finalINFO$seqnum[i]<-fileinfo[4]
  
  # finalINFO$gene[i]<-fileinfo[2]
  # finalINFO$orf[i]<-geneInfo$geneID[geneInfo$geneName==gene]
  # finalINFO$crispr[i]<-"a"
  # finalINFO$rep[i]<-fileinfo[3]
  # finalINFO$tech[i]<-fileinfo[4]
  # finalINFO$gate[i]<-fileinfo[5]
  # finalINFO$seqnum[i]<-fileinfo[6]
  
  resultFileName<-paste(finalINFO$orf[i],finalINFO$gene[i],finalINFO$crispr[i],finalINFO$rep[i],finalINFO$tech[i],finalINFO$gate[i],finalINFO$seqnum[i],sep="_")
  
  FileVCF <- readVcf(File,"sacCer3")
  FileVCFCounts <- t(sapply(info(FileVCF)$AD, function(x){c(x[[1]], x[[2]])}))
  rownames(FileVCFCounts) <- sapply(names(rowRanges(FileVCF)), function(x){strsplit(x, "_")[[1]][1]})
  # these counts were read from the vcf, so they don't contain counts for sites without reads
  # need to build a common SNP table
  theseCounts <- data.frame(SNPs, matrix(0, nrow=nrow(SNPs), ncol=2))
  rownames(theseCounts) <- paste(theseCounts[,1], theseCounts[,2], sep=":")
  colnames(theseCounts) <- c("chr", "pos", "ref", "alt")
  theseCounts[rownames(FileVCFCounts), c("ref", "alt")] <- FileVCFCounts
  
  
  gcoords= getGcoords(theseCounts$chr, theseCounts$pos, sepBetweenChr, sgd.table=sgd_table)
  names(gcoords) = rownames(theseCounts)
  
  # get plot coordinates for the chromosome line dividers
  chrCutoffs <- sapply(unique(theseCounts$chr), function(x){gcoords[theseCounts$chr == x][1] - sepBetweenChr/2})
  names(chrCutoffs) <- unique(theseCounts$chr)
  chrLabels = sapply(1:(length(chrCutoffs)-1), function(i)(chrCutoffs[i] + chrCutoffs[i+1])/2)
  # add half the length of chrXVI
  chrLabels = c(chrLabels, chrCutoffs[16] + sepBetweenChr + 948066/2)
  names(chrLabels)[16] = "chrXVI"
  
  
  # compute BY allele frequencies
  theseCounts$BYAF <- theseCounts$ref / (theseCounts$ref + theseCounts$alt)
  theseCounts$Coverage <- theseCounts$ref + theseCounts$alt
  theseCounts$roll = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"], theseCounts[,"BYAF"], gcoords, median(theseCounts[,"Coverage"]), stringsAsFactors=FALSE), LoessSpan)
  finalINFO$cover[i]<-median(theseCounts[,"Coverage"])
  
  # write out the counts to be analyzed
  save(theseCounts, file=paste(resultsFolder,resultFileName, ".RData", sep=""))
  
  ###just plot for control
  png(filename = paste(resultsFolder,resultFileName,".png",sep=""),width = 500,height = 340)
  #pdf(file = "test.pdf",width = 10,height = 6.8)
  plot(gcoords, rep(0.5, length(gcoords)), ylim=c(0,1), xaxt='n', xlab="chromosome", ylab="BY allele frequency", type="n")
  
  # CAN1 is on chr05, pos 33466
  # MATALPHA is on chr03, 198671
  abline(v = getGcoords("chrV", 33466, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="green")
  abline(v = getGcoords("chrIII", 198671, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="green")
  abline(v = getGcoords(geneInfo[finalINFO$orf[i], "chr"], mean(as.numeric(geneInfo[finalINFO$orf[i], c("start", "end")])), sepBetweenChr, sgd.table=sgd_table), lwd=2, col="purple")
  
  points(gcoords, theseCounts[,"BYAF"], cex=.2, col="#00000022")
  
  
  for (j in unique(theseCounts[,"chr"])){
    points(gcoords[theseCounts$chr == j], theseCounts$roll[theseCounts$chr == j], type="l", lwd=2)
  }
  
  for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
  }
  abline(h = 0.5, lty=2, col="light blue")
  
  legend("topright", legend = c(paste0("median coverage: ", median(theseCounts[,"Coverage"]))), box.lty=0, bg = rgb(0,0,0,0))
  axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)
  dev.off()
}

write.table(x = finalINFO,file = paste(resultsFolder, "FinalINFO.txt", sep=""),quote = F,sep = "\t",col.names = T,row.names = F)


