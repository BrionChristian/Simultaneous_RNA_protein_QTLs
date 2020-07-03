alignmentDir <- "C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/190208-AllQTLprocess/"
setwd(alignmentDir)

#library("VariantAnnotation")
source(paste(alignmentDir,"/Scripts/gTest.R",sep=""))
source(paste(alignmentDir,"/Scripts/x_qtl_seq_functions_170831.R",sep=""))
source(paste(alignmentDir,"/Scripts/mp_JB_170901.R",sep=""))
source(paste(alignmentDir,"/Scripts/peaksFromVector.R",sep=""))



# can then call this from external
i = 1
sepBetweenChr <- 1e5
trimFromEnd = 15e3
obsMin <- 10
LoessSpan = 0.1
AFThres = 0.09653124 # same as in Albert 2014
multiThres = 2.8 # LOD threshold for multipool, if run with N=1000: 4.5
withMultipool <- TRUE

SNPs <- read.table(paste(alignmentDir,"/Scripts/SNPs_Maggie_170809_BY_positions.txt",sep=""), stringsAsFactors=FALSE, head=FALSE)
# see comments above. As of 8/31/17, the SNPs seem not to be fully filtered, and are out of sorting order
for (thisChr in unique(SNPs[,1])){SNPs[SNPs[,1] == thisChr, 2] <- sort(SNPs[SNPs[,1] == thisChr, 2])}
SNPs <- rbind(SNPs[SNPs[,1] == "chrI",], SNPs[SNPs[,1] == "chrII",], SNPs[SNPs[,1] == "chrIII",], SNPs[SNPs[,1] == "chrIV",], SNPs[SNPs[,1] == "chrV",], SNPs[SNPs[,1] == "chrVI",], SNPs[SNPs[,1] == "chrVII",], SNPs[SNPs[,1] == "chrVIII",], SNPs[SNPs[,1] == "chrIX",], SNPs[SNPs[,1] == "chrX",], SNPs[SNPs[,1] == "chrXI",], SNPs[SNPs[,1] == "chrXII",], SNPs[SNPs[,1] == "chrXIII",], SNPs[SNPs[,1] == "chrXIV",], SNPs[SNPs[,1] == "chrXV",], SNPs[SNPs[,1] == "chrXVI",])
sgd_table <- paste(alignmentDir,"/Scripts/sacCer3ChromLenghts.txt",sep="")
geneInfo = read.table(paste(alignmentDir,"/Scripts/ensemblGenes_ensembl83_160307_MOD.txt",sep=""), stringsAsFactors=FALSE, sep="\t", header=TRUE)
rownames(geneInfo) <- geneInfo[,"geneID"]
allNames <- geneInfo[, "geneName"]
names(allNames) <- geneInfo[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames
chromosome<-unique(SNPs[,1])


#frank data
xpQTL <- read.table("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/Biblio/Frank_eQTL_pQTL/xpQTLs.txt",header = T,sep = "\t", quote = "", na.strings = "NA",stringsAsFactors = F)
final_eQTL <- read.table("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/eQTLFrank/final_eQTL.txt",header = T,sep = "\t", quote = "", na.strings = "NA",stringsAsFactors = F,comment.char = "")
geneQTG <- read.table("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/eQTLFrank/genes_RMBY_variant.txt",header = T,sep = "\t", quote = "", na.strings = "",stringsAsFactors = F,comment.char = "")
geneQTG<-geneQTG[!(is.na(geneQTG$chr)),]

experimentFile <- read.table(paste(alignmentDir,"AllPoprecap.txt",sep=""), stringsAsFactors=FALSE, head=TRUE)

experimentFile$tube<-paste(experimentFile$gene,experimentFile$crispr,experimentFile$rep,experimentFile$tech,experimentFile$gate,sep="_")
experimentFile2<-experimentFile[experimentFile$tube!= "GPD1_a_3_2",] #experimentFile$crispr!='i' & 
levels(as.factor(experimentFile2$tube))

#resultsFolder <- paste(alignmentDir,"/SeqOutput/180515folder/results/",sep="")




experimentFile3<-experimentFile2[experimentFile2$crispr!="i",]
experimentFile3<-experimentFile3[paste(experimentFile3$gene,experimentFile3$rep)!="RPS10A 2",]
experimentFile3<-experimentFile3[experimentFile3$gene!="TDH3",]
experimentFile3<-experimentFile3[experimentFile3$gate<6,]


ref1<-c()
ref2<-c()
n=1
for (i in 1:nrow(experimentFile3)){
  tempref<-experimentFile3[i,]
  tempfile<-experimentFile3[experimentFile3$gene==tempref$gene & experimentFile3$gate == tempref$gate & experimentFile3$rep+experimentFile3$tech/10>=tempref$rep+tempref$tech/10,]
  if (nrow(tempfile)==0) {
    next
  } else {
    for (j in 1:nrow(tempfile)) {
      ref1[n]<-tempref$tube
      ref2[n]<-tempfile$tube[j]
      n=n+1
    }
  }

}
refinfo<-data.frame(ref1,ref2)
refinfo<-refinfo[refinfo$ref1!=refinfo$ref2,]
refinfo$conca<-paste(refinfo$ref1,refinfo$ref2)
levels(as.factor(refinfo$conca))

initialtable <- read.table(paste(alignmentDir,"ReplicateResult/ARO8_a_0_1_1 ARO8_a_1_1_1_QTLs.txt",sep=""), stringsAsFactors=FALSE, head=TRUE,sep = "\t")
initialtable <- initialtable[-(1:nrow(initialtable)),]

i=1



for (i in 1:nrow(refinfo)) { #loop every 5min
  

experiment<-refinfo$conca[i]
print(paste(experiment,date(),paste(i,"/",nrow(refinfo),sep=""),sep=" - "))


#getdata
tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(refinfo$ref1[i], "_.*RData$"), full.names=TRUE)
load(file = tempFile)
G1<-theseCounts
tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(refinfo$ref2[i], "_.*RData$"), full.names=TRUE)
load(file = tempFile)
G2<-theseCounts


#do multipool
# date()
# multipoolOutput <- lapply(unique(G1[,"chr"]), function(j){
#   doMultiPoolFromWithinR(G1[G1$chr == j, c("pos", "ref", "alt")], G2[G2$chr == j, c("pos", "ref", "alt")])
# })
# date() #9min
# 
# # call peaks from multipool LODs
# multiPeaks <- lapply(multipoolOutput, function(x){
#   thisLODTrace = x[[2]]
#   theseChrPeaks = callPeaks(thisLODTrace[,2], multiThres, 2)
#   theseChrPeaks
# })
# 
# save(multiPeaks, multipoolOutput, file=paste(alignmentDir,"ReplicateResult/",paste(experiment, sep="_"), "_multipoolResults.RData", sep=""))
#####

#load multipooldata
tempFile <- dir(paste(alignmentDir,"ReplicateResult/",sep=""), pattern=paste0(experiment, "_multi.*RData$"), full.names=TRUE)
load(file = tempFile)

#getQTL

curentEXP<-refinfo[i,]
exportpeak<-cbind(curentEXP,c(1),c(1),c(1),c(1),c(1))
colnames(exportpeak) = c(colnames(curentEXP),"chromo","maxIndex", "maxValue", "leftIndex", "rightIndex")
exportpeak<-exportpeak[-1,]
for (chr in 1:length(multiPeaks)) {
  if (is.null(multiPeaks[[chr]])) {next}
  for (peak in 1:nrow(multiPeaks[[chr]])) {
    exportpeak<-rbind(exportpeak,c(curentEXP,"chromo"=chr,multiPeaks[[chr]][peak,]))
  }
}
colF <- sapply(exportpeak, is.factor)
exportpeak[colF] <- lapply(exportpeak[colF], as.character)


if (nrow(exportpeak)>0) {
  deltaHight_Low<-c()
  for (j in (1:nrow(exportpeak))) {
    HLcurve<-G1$roll-G2$roll
    HLcurve<-HLcurve[G1$chr==chromosome[exportpeak$chromo[j]]]
    HLcurve2<-HLcurve[!(is.na(HLcurve))]
    postemp<-G1$pos[G1$chr==chromosome[exportpeak$chromo[j]]]
    postemp<-postemp[!(is.na(HLcurve))]
    bestpos<-abs(postemp-exportpeak$maxIndex[j])==min(abs(postemp-exportpeak$maxIndex[j]))
    deltaHight_Low[j]<-HLcurve2[bestpos][1]
  }
  exportpeak$deltaHight_Low<-deltaHight_Low
} else {
  exportpeak$deltaHight_Low<-c()
}

write.table(exportpeak,paste(alignmentDir,"ReplicateResult/",paste(experiment, sep="_"), "_QTLs.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)

initialtable<-rbind(initialtable,exportpeak)

#}

#now ploting


pdf(file = paste(alignmentDir,"ReplicateResult/",paste(experiment, sep="_"), "_plot.pdf", sep=""), width=11, height=16)
par(mfrow=c(3,1))
# get plot coordinates for the chromosome line dividers
gcoords= getGcoords(G1$chr, G1$pos, sepBetweenChr, sgd.table=sgd_table)
names(gcoords) = rownames(G1)
chrCutoffs <- sapply(unique(G1$chr), function(x){gcoords[G1$chr == x][1] - sepBetweenChr/2})
names(chrCutoffs) <- unique(G1$chr)
chrLabels = sapply(1:(length(chrCutoffs)-1), function(i)(chrCutoffs[i] + chrCutoffs[i+1])/2)
# add half the length of chrXVI
chrLabels = c(chrLabels, chrCutoffs[16] + sepBetweenChr + 948066/2)
names(chrLabels)[16] = "chrXVI"

plot(gcoords, rep(0.5, length(gcoords)), ylim=c(0,1), main = paste(experiment, sep="_"), xaxt='n', xlab="chromosome", ylab="BY allele frequency", type="n")
abline(v = getGcoords("chrV", 33466, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # CAN1 is on chr05, pos 33466
abline(v = getGcoords("chrIII", 198671, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # MATALPHA is on chr03, 198671
for (j in unique(G1[,"chr"])){
  points(gcoords[G1$chr == j], G1$roll[G1$chr == j], type="l", lwd=2,col="black",lty=1)
  points(gcoords[G2$chr == j], G2$roll[G2$chr == j], type="l", lwd=2,col="red",lty=1)
}
for (j in chrCutoffs){
  abline(v = j, lty=2, col="light blue")
}
abline(h = 0.5, lty=2, col="light blue")
axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)

plot(gcoords, rep(0, length(gcoords)), ylim=c(-0.4,0.4), main = paste(experiment, sep="_"), xaxt='n', xlab="chromosome", ylab="BY allele frequency", type="n")
#abline(v = getGcoords("chrV", 33466, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # CAN1 is on chr05, pos 33466
#abline(v = getGcoords("chrIII", 198671, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # MATALPHA is on chr03, 198671
for (j in unique(G1[,"chr"])){
  points(gcoords[G1$chr == j], G1$roll[G1$chr == j]-G2$roll[G1$chr == j], type="l", lwd=2,col="red",lty=1)
}
for (j in chrCutoffs){
  abline(v = j, lty=2, col="light blue")
}
abline(h = 0, lty=2, col="light blue")
# these difference thresholds from the nullSorts:
abline(h = AFThres, lty = 2, lwd=2, col = c(1,1,1,0.8))
abline(h = -AFThres, lty = 2, lwd=2, col = c(1,1,1,0.8))
axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)


ylimMax = max(c(multiThres, sapply(multipoolOutput, function(x){max(x[[2]][,2])}))) + 1
ylimMax <- min(ylimMax,70)
plot(gcoords, rep(0, length(gcoords)), main = paste(experiment, sep="_"), type="n", xaxt='n', xlab="chromosome", ylab="Multipool LOD", ylim=c(0, ylimMax))
for (j in 1:16){
  points(getGcoords(paste0("chr", as.roman(j)), multipoolOutput[[j]][[2]][,1], sepBetweenChr, sgd.table=sgd_table), multipoolOutput[[j]][[2]][,2], type="l", lwd=2, col=rgb(0,0.8,0.5))
  # add stars at peaks
  if(!is.null(multiPeaks[[j]])){
    for (thisPeak in 1:nrow(multiPeaks[[j]])){
      thisPeakPlotPos <- getGcoords(j, multiPeaks[[j]][thisPeak, "maxIndex"], sepBetweenChr, sgd.table=sgd_table)
      text(thisPeakPlotPos, multiPeaks[[j]][thisPeak, "maxValue"] + 0.5, labels="*", col=rgb(0,0.8,0.5), cex=3)
    }
  }
}
abline(h = 0, lty=2, col="light blue")
for (j in chrCutoffs){
  abline(v = j, lty=2, col="light blue")
}
axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)
abline(h = multiThres, lty=2, col="red", lwd=2)


dev.off()

}

write.table(initialtable,paste(alignmentDir,"All_rep.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)


#===========================================================



trueQTL<-read.table("All_QTLs.txt",header = T,sep = "\t",stringsAsFactors = F)

refinfo2<-refinfo
refinfo2<-refinfo[-grep(pattern = "_0_",x = refinfo$conca),]
#refinfo2<-refinfo2[-grep(pattern = "TPO1",x = refinfo2$conca),] #just to test
refinfo3<-refinfo[grep(pattern = "aF_0_",x = refinfo$ref1),]
refinfo3<-refinfo3[-grep(pattern = "a_1_",x = refinfo3$conca),]
refinfo2<-rbind(refinfo2,refinfo3)
nf<-length(levels(as.factor(refinfo2$conca)))

length(levels(as.factor(initialtable$conca)))
levels(as.factor(initialtable$conca))
initialtable2<-initialtable
initialtable2<-initialtable[-grep(pattern = "_0_",x = initialtable$conca),]
#initialtable2<-initialtable2[-grep(pattern = "TPO1",x = initialtable2$conca),]#just to test
initialtable3<-initialtable[grep(pattern = "aF_0_",x = initialtable$ref1),]
initialtable3<-initialtable3[-grep(pattern = "a_1_",x = initialtable3$conca),]
initialtable2<-rbind(initialtable2,initialtable3)
nrow(initialtable2[initialtable2$maxValue>3,])

trueQTL2<-trueQTL[trueQTL$crispr %in% c("a","aF") & trueQTL$gene != "TDH3" & paste(trueQTL$gene,trueQTL$rep)!="RPS10A 2",]
nt<-(length(levels(as.factor(paste(trueQTL2$gene,trueQTL2$crispr,trueQTL2$rep,trueQTL2$tech)))))*2

x<-0:400/20
rank<-1:401
yt<-0:400
yf<-0:400
for (i in rank) {
  yf[i]<-nrow(initialtable2[initialtable2$maxValue>x[i],])
  yt[i]<-nrow(trueQTL2[trueQTL2$maxValue>x[i],])
}

Hi<-5
Wi<-6

svg(filename = "FDRcount.svg",width = Wi,height = Hi)
plot(x,yt/nt,ylim=c(0,max(yt/nt)),type="l",xlab="lod threshold",ylab = "number of QTL")
points(x,yf/nf,type="l",col="red")
dev.off()

svg(filename = "FDRratio.svg",width = Wi,height = Hi)
plot(x,(yf/nf)/(yt/nt),ylim=c(0,0.2),type="l",xlab="lod threshold",ylab = "FDR")
dev.off()

#plot(x,yt/nt,ylim=c(0,max(yt/nt)),xlim=c(3,6),type="l",xlab="lod threshold",ylab = "number of QTL")
#points(x,yf/nf,ylim=c(0,max(yt/nt)),type="l",col="red")

FDR<-(yf/nf)/(yt/nt)

thr1=3.5
thr2=4.5
FDR1<-FDR[x==thr1]
FDR2<-FDR[x==thr2]

svg(filename = "FDRratiozoom.svg",width = Wi,height = Hi)
plot(x,(yf/nf)/(yt/nt),xlim=c(3,6),ylim=c(0.04,0.14),type="l",xlab="lod threshold",ylab = "FDR")
lines(c(0,thr1),c(FDR1,FDR1),col="grey")
lines(c(thr1,thr1),c(0,FDR1),col="grey")
lines(c(0,thr2),c(FDR2,FDR2),col="red",lwd=1.5)
lines(c(thr2,thr2),c(0,FDR2),col="red",lwd=1.5)
dev.off()

initialtable2tech<-initialtable[c(grep(pattern = "UGP1_aF_0_1_[0-9] UGP1_a_0_1_[0-9]",x = initialtable$conca),
                                  grep(pattern = "GPD1_a_1_1_[0-9] GPD1_a_1_2_[0-9]",x = initialtable$conca),
                                  grep(pattern = "GPD1_a_3_1_[0-9] GPD1_a_3_2_[0-9]",x = initialtable$conca)),]


