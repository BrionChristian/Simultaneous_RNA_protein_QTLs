# Author: Christian Brion - 2020 - UMN
#
#-perform the multipool bulk QTL analysis experiment per experiment
#-save raw QTL data in txt file
#-save processed comparative data in R files
#-provide figures for each QTL analysis
#


alignmentDir <- "/home/christian/Dropbox/Gdrive/MN_postdoc/diverR/190208-AllQTLprocess/"
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
multiThres = 3 # LOD threshold for multipool, if run with N=1000: 4.5
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
xpQTL <- read.table("/home/christian/Dropbox/Gdrive/MN_postdoc/Biblio/Frank_eQTL_pQTL/xpQTLs.txt",header = T,sep = "\t", quote = "", na.strings = "NA",stringsAsFactors = F)
final_eQTL <- read.table("/home/christian/Dropbox/Gdrive/MN_postdoc/diverR/eQTLFrank/final_eQTL.txt",header = T,sep = "\t", quote = "", na.strings = "NA",stringsAsFactors = F,comment.char = "")
geneQTG <- read.table("/home/christian/Dropbox/Gdrive/MN_postdoc/diverR/eQTLFrank/genes_RMBY_variant.txt",header = T,sep = "\t", quote = "", na.strings = "",stringsAsFactors = F,comment.char = "")
geneQTG<-geneQTG[!(is.na(geneQTG$chr)),]

experimentFile <- read.table(paste(alignmentDir,"/Scripts/AllPoprecap.txt",sep=""), stringsAsFactors=FALSE, head=TRUE)

experimentFile$tube<-paste(experimentFile$gene,experimentFile$crispr,experimentFile$rep,experimentFile$tech,sep="_")
experimentFile2<-experimentFile[experimentFile$tube!= "GPD1_a_3_2",] #experimentFile$crispr!='i' & 
levels(as.factor(experimentFile2$tube))

#resultsFolder <- paste(alignmentDir,"/SeqOutput/180515folder/results/",sep="")

initialtable <- read.table(paste(alignmentDir,"QTLresults/YDL022W_GPD1_a_0_1_QTLs.txt",sep=""), stringsAsFactors=FALSE, head=TRUE)
initialtable <- initialtable[-(1:nrow(initialtable)),]


for (i in 1:length(levels(as.factor(experimentFile2$tube)))) { #loop every 5min
  

experiment<-levels(as.factor(experimentFile2$tube))[i]
orf<-experimentFile2$orf[experimentFile2$tube==experiment][1]
print(paste(orf,experiment,date(),paste(i,"/",length(levels(as.factor(experimentFile2$tube))),sep=""),sep=" - "))


#getdata
tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment, "_1_.*RData$"), full.names=TRUE)
load(file = tempFile)
G1<-theseCounts
tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment, "_2_.*RData$"), full.names=TRUE)
load(file = tempFile)
G2<-theseCounts
tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment, "_3_.*RData$"), full.names=TRUE)
load(file = tempFile)
G3<-theseCounts
tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment, "_4_.*RData$"), full.names=TRUE)
load(file = tempFile)
G4<-theseCounts
tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment, "_5_.*RData$"), full.names=TRUE)
load(file = tempFile)
G5<-theseCounts


#do multipool
# date()
# multipoolOutputGFP <- lapply(unique(G1[,"chr"]), function(j){
#   doMultiPoolFromWithinR(G5[G5$chr == j, c("pos", "ref", "alt")], G4[G4$chr == j, c("pos", "ref", "alt")])
# })
# date() #4min
# 
# # call peaks from multipool LODs
# multiPeaksGFP <- lapply(multipoolOutputGFP, function(x){
#   thisLODTrace = x[[2]]
#   theseChrPeaks = callPeaks(thisLODTrace[,2], multiThres, 2)
#   theseChrPeaks
# })
# 
# date()
# multipoolOutputmCH <- lapply(unique(G1[,"chr"]), function(j){
#   doMultiPoolFromWithinR(G3[G3$chr == j, c("pos", "ref", "alt")], G2[G2$chr == j, c("pos", "ref", "alt")])
# })
# date() #4min
# 
# # call peaks from multipool LODs
# multiPeaksmCH <- lapply(multipoolOutputmCH, function(x){
#   thisLODTrace = x[[2]]
#   theseChrPeaks = callPeaks(thisLODTrace[,2], multiThres, 2)
#   theseChrPeaks
# })
# 
# save(multiPeaksGFP, multipoolOutputGFP, multiPeaksmCH,multipoolOutputmCH, file=paste(alignmentDir,"QTLresults/",paste(orf,experiment, sep="_"), "_multipoolResults.RData", sep=""))
######

#load multipooldata
tempFile <- dir(paste(alignmentDir,"QTLresults/",sep=""), pattern=paste0(experiment, "_multi.*RData$"), full.names=TRUE)
load(file = tempFile)

#getQTL
exportpeak<-experimentFile[1,-c(6,7,8)]
curentEXP<-experimentFile[experimentFile$tube==experiment,-c(6,7,8)][1,]
exportpeak<-cbind(exportpeak,c("gfp"),c(1),c(1),c(1),c(1),c(1))
colnames(exportpeak) = c(colnames(exportpeak)[1:7],"fluorecence","chromo","maxIndex", "maxValue", "leftIndex", "rightIndex")
exportpeak$fluorecence<-as.character(exportpeak$fluorecence)
exportpeak<-exportpeak[-1,]
for (chr in 1:length(multiPeaksGFP)) {
  if (is.null(multiPeaksGFP[[chr]])) {next}
  for (peak in 1:nrow(multiPeaksGFP[[chr]])) {
    exportpeak<-rbind(exportpeak,c(curentEXP,"fluorecence"="gfp","chromo"=chr,multiPeaksGFP[[chr]][peak,]))
  }
}
colF <- sapply(exportpeak, is.factor)
exportpeak[colF] <- lapply(exportpeak[colF], as.character)
for (chr in 1:length(multiPeaksmCH)) {
  if (is.null(multiPeaksmCH[[chr]])) {next}
  for (peak in 1:nrow(multiPeaksmCH[[chr]])) {
    exportpeak<-rbind(exportpeak,c(curentEXP,"fluorecence"="mch","chromo"=chr,multiPeaksmCH[[chr]][peak,]))
  }
}

deltaHight_Low<-c()
deltaOther<-c()
LOD_Other<-c()
effectRNAPROT<-c()
for (j in (1:nrow(exportpeak))) {
  if (exportpeak$fluorecence[j] =="gfp") {
    HLcurve<-G5$roll-G4$roll
    HLcurveOther<-G3$roll-G2$roll
    lodOther<-multipoolOutputmCH[[exportpeak$chromo[j]]][[2]]
  } else {
    HLcurve<-G3$roll-G2$roll
    HLcurveOther<-G5$roll-G4$roll
    lodOther<-multipoolOutputGFP[[exportpeak$chromo[j]]][[2]]
  }
  HLcurve<-HLcurve[G1$chr==chromosome[exportpeak$chromo[j]]]
  HLcurve2<-HLcurve[!(is.na(HLcurve))]
  postemp<-G1$pos[G1$chr==chromosome[exportpeak$chromo[j]]]
  postemp<-postemp[!(is.na(HLcurve))]
  HLcurveOther<-HLcurveOther[G1$chr==chromosome[exportpeak$chromo[j]]]
  HLcurveOther2<-HLcurveOther[!(is.na(HLcurveOther))]
  postempOther<-G1$pos[G1$chr==chromosome[exportpeak$chromo[j]]]
  postempOther<-postempOther[!(is.na(HLcurveOther))]
  bestpos<-abs(postemp-exportpeak$maxIndex[j])==min(abs(postemp-exportpeak$maxIndex[j]))
  bestposOther<-abs(postempOther-exportpeak$maxIndex[j])==min(abs(postempOther-exportpeak$maxIndex[j]))
  deltaHight_Low[j]<-HLcurve2[bestpos][1]
  deltaOther[j]<-HLcurveOther2[bestposOther][1]
  LOD_Other[j]<-lodOther[lodOther[,1]==exportpeak$maxIndex[j],2]
  if (LOD_Other[j]<3) {
    effectRNAPROT[j]<-"specific"
  } else {
    if (sign(deltaHight_Low[j])==sign(deltaOther[j])) {
      if (exportpeak$crispr[j]=="i") {
        effectRNAPROT[j]<-"oposite"
      } else {
        effectRNAPROT[j]<-"together"
      }
    } else {
      if (exportpeak$crispr[j]=="i") {
        effectRNAPROT[j]<-"together"
      } else {
        effectRNAPROT[j]<-"oposite"
      }
    }
  }
}
exportpeak$deltaHight_Low<-deltaHight_Low
exportpeak$deltaOther<-deltaOther
exportpeak$LOD_Other<-LOD_Other
exportpeak$effectRNAPROT<-effectRNAPROT

write.table(exportpeak,paste(alignmentDir,"QTLresults/",paste(orf,experiment, sep="_"), "_QTLs.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)

initialtable<-rbind(initialtable,exportpeak)

#}

#now ploting


pdf(file = paste(alignmentDir,"QTLresults/",paste(orf,experiment, sep="_"), "_plot.pdf", sep=""), width=11, height=16)
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

plot(gcoords, rep(0.5, length(gcoords)), ylim=c(0,1), main = paste(orf,experiment, sep="_"), xaxt='n', xlab="chromosome", ylab="BY allele frequency", type="n")
abline(v = getGcoords("chrV", 33466, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # CAN1 is on chr05, pos 33466
abline(v = getGcoords("chrIII", 198671, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # MATALPHA is on chr03, 198671
abline(v = getGcoords(geneInfo[orf, "chr"], mean(as.numeric(geneInfo[orf, c("start", "end")])), sepBetweenChr, sgd.table=sgd_table), lwd=2, col="purple")
for (j in unique(G1[,"chr"])){
  points(gcoords[G1$chr == j], G1$roll[G1$chr == j], type="l", lwd=2,col="black",lty=1)
  points(gcoords[G2$chr == j], G2$roll[G2$chr == j], type="l", lwd=1.3,col="red",lty=2)
  points(gcoords[G3$chr == j], G3$roll[G3$chr == j], type="l", lwd=2,col="red",lty=1)
  points(gcoords[G4$chr == j], G4$roll[G4$chr == j], type="l", lwd=1.3,col=rgb(0,0.8,0.5),lty=2)
  points(gcoords[G5$chr == j], G5$roll[G5$chr == j], type="l", lwd=2,col=rgb(0,0.8,0.5),lty=1)
}
for (j in chrCutoffs){
  abline(v = j, lty=2, col="light blue")
}
abline(h = 0.5, lty=2, col="light blue")
axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)

plot(gcoords, rep(0, length(gcoords)), ylim=c(-1,1), main = paste(orf,experiment, sep="_"), xaxt='n', xlab="chromosome", ylab="BY allele frequency", type="n")
#abline(v = getGcoords("chrV", 33466, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # CAN1 is on chr05, pos 33466
#abline(v = getGcoords("chrIII", 198671, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # MATALPHA is on chr03, 198671
abline(v = getGcoords(geneInfo[orf, "chr"], mean(as.numeric(geneInfo[orf, c("start", "end")])), sepBetweenChr, sgd.table=sgd_table), lwd=2, col="purple")
for (j in unique(G1[,"chr"])){
  points(gcoords[G3$chr == j], G3$roll[G3$chr == j]-G2$roll[G3$chr == j], type="l", lwd=2,col="red",lty=1)
  points(gcoords[G5$chr == j], G5$roll[G3$chr == j]-G4$roll[G5$chr == j], type="l", lwd=2,col=rgb(0,0.8,0.5),lty=1)
}
for (j in chrCutoffs){
  abline(v = j, lty=2, col="light blue")
}
abline(h = 0, lty=2, col="light blue")
# these difference thresholds from the nullSorts:
abline(h = AFThres, lty = 2, lwd=2, col = c(1,1,1,0.8))
abline(h = -AFThres, lty = 2, lwd=2, col = c(1,1,1,0.8))
this_eQTL<- final_eQTL[final_eQTL$gene==orf & final_eQTL$cis == F,]
if (nrow(this_eQTL)!=0) {
  for (j in 1:nrow(this_eQTL)) {
    abline(v = getGcoords(this_eQTL$chr.x[j], this_eQTL$maxQTLpos[j], sepBetweenChr, sgd.table=sgd_table), lwd=2, col="red")
    points(getGcoords(this_eQTL$chr.x[j], this_eQTL$maxQTLpos[j], sepBetweenChr, sgd.table=sgd_table),-this_eQTL$r[j],pch = 3, lwd=2, col="red")
  }
}
this_pQTL<- xpQTL[xpQTL$gene==orf,]
this_pQTL$chr.x<-as.numeric(gsub(pattern = "chr",replacement = "",this_pQTL$chromosome))
if (nrow(this_pQTL)!=0) {
  for (j in 1:nrow(this_pQTL)) {
    abline(v = getGcoords(chromosome[this_pQTL$chr.x[j]], this_pQTL$peakPosition[j], sepBetweenChr, sgd.table=sgd_table), lwd=2, col=rgb(0,0.8,0.5))
    points(getGcoords(chromosome[this_pQTL$chr.x[j]], this_pQTL$peakPosition[j], sepBetweenChr, sgd.table=sgd_table),this_pQTL$alleleFrequencyDifference[j],pch = 3, lwd=2, col=rgb(0,0.8,0.5))
  }
}
axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)


ylimMaxGFP = max(c(multiThres, sapply(multipoolOutputGFP, function(x){max(x[[2]][,2])}))) + 1
ylimMaxmch = max(c(multiThres, sapply(multipoolOutputmCH, function(x){max(x[[2]][,2])}))) + 1
ylimMax <- max(ylimMaxGFP,ylimMaxmch)
ylimMax <- min(ylimMax,70)
plot(gcoords, rep(0, length(gcoords)), main = paste(orf,experiment, sep="_"), type="n", xaxt='n', xlab="chromosome", ylab="Multipool LOD", ylim=c(0, ylimMax))
abline(v = getGcoords(geneInfo[orf, "chr"], mean(as.numeric(geneInfo[orf, c("start", "end")])), sepBetweenChr, sgd.table=sgd_table), lwd=2, col="purple")
for (j in 1:16){
  points(getGcoords(paste0("chr", as.roman(j)), multipoolOutputGFP[[j]][[2]][,1], sepBetweenChr, sgd.table=sgd_table), multipoolOutputGFP[[j]][[2]][,2], type="l", lwd=2, col=rgb(0,0.8,0.5))
  # add stars at peaks
  if(!is.null(multiPeaksGFP[[j]])){
    for (thisPeak in 1:nrow(multiPeaksGFP[[j]])){
      thisPeakPlotPos <- getGcoords(j, multiPeaksGFP[[j]][thisPeak, "maxIndex"], sepBetweenChr, sgd.table=sgd_table)
      text(thisPeakPlotPos, multiPeaksGFP[[j]][thisPeak, "maxValue"] + 0.5, labels="*", col=rgb(0,0.8,0.5), cex=3)
    }
  }
  points(getGcoords(paste0("chr", as.roman(j)), multipoolOutputmCH[[j]][[2]][,1], sepBetweenChr, sgd.table=sgd_table), multipoolOutputmCH[[j]][[2]][,2], type="l", lwd=2, col="red")
  # add stars at peaks
  if(!is.null(multiPeaksmCH[[j]])){
    for (thisPeak in 1:nrow(multiPeaksmCH[[j]])){
      thisPeakPlotPos <- getGcoords(j, multiPeaksmCH[[j]][thisPeak, "maxIndex"], sepBetweenChr, sgd.table=sgd_table)
      text(thisPeakPlotPos, multiPeaksmCH[[j]][thisPeak, "maxValue"] + 0.5, labels="*", col="red", cex=3)
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

write.table(initialtable,paste(alignmentDir,"All_QTLs.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)

summary(as.factor(initialtable$effectRNAPROT))
summary(as.factor(initialtable$fluorecence))
summary(as.factor(initialtable$effectRNAPROT[initialtable$crispr!="i"]))

finalTable<-initialtable[initialtable$crispr!="i" & !(initialtable$tube%in%c("RPS10A_a_2_1","TDH3_a_0_1")),]
#finalTable<-initialtable
summary(as.factor(finalTable$effectRNAPROT))
summary(as.factor(finalTable$fluorecence))
summary(as.factor(paste(finalTable$fluorecence,finalTable$effectRNAPROT)))


cexall<-sapply(finalTable$maxValue,FUN =  function(x) {max(0.4,min(3,x/20))})
finalTable2<-finalTable
finalTable2$cexall<-cexall
finalTable45<-finalTable2[finalTable2$maxValue>4.5 & abs(finalTable2$deltaHight_Low)>0.05,]
summary(as.factor(paste(finalTable45$fluorecence,finalTable45$effectRNAPROT)))
finalTable45rep1<-finalTable45[finalTable45$rep==1 & finalTable45$tech==1,]
summary(as.factor(paste(finalTable45rep1$fluorecence,finalTable45rep1$effectRNAPROT)))
finalTable2<-finalTable45
finalTable2<-finalTable45rep1
summary(as.factor(paste(finalTable2$fluorecence,finalTable2$effectRNAPROT)))
plot(finalTable2$deltaHight_Low[finalTable2$fluorecence=="gfp"],finalTable2$deltaOther[finalTable2$fluorecence=="gfp"],
     xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0.8,0.5),cex=finalTable2$cexall[finalTable2$fluorecence=="gfp"],
     xlab="AF(BY/all) difference (high-low) - GFP gate", ylab="AF(BY/all) difference (high-low) - mCherry gate",
     main="QTL effect: protein vs mRNA")
points(finalTable2$deltaOther[finalTable2$fluorecence=="mch"],finalTable2$deltaHight_Low[finalTable2$fluorecence=="mch"],xlim=c(-1,1),ylim=c(-1,1),col="red",cex=finalTable2$cexall[finalTable2$fluorecence=="mch"])
abline(v=0)
abline(h=0)
#subtable<-finalTable2[finalTable2$gene=="GPD1" & finalTable2$chromo == 10 & finalTable2$maxIndex > 100000  & finalTable2$maxIndex < 200000,]
subtable<-finalTable2[finalTable2$chromo == 10 & finalTable2$maxIndex > 100000  & finalTable2$maxIndex < 200000,] #TIF2
#subtable<-finalTable2[finalTable2$chromo == 12 & finalTable2$maxIndex > 620000  & finalTable2$maxIndex < 680000,] #HAP1
#subtable<-finalTable2[finalTable2$chromo == 14 & finalTable2$maxIndex > 440000  & finalTable2$maxIndex < 500000,] #MKT1
#subtable<-finalTable2[finalTable2$chromo == 4 & finalTable2$maxIndex > 1000000  & finalTable2$maxIndex < 1100000,] #GCN2
#subtable<-finalTable2[finalTable2$chromo == 15 & finalTable2$maxIndex > 140000  & finalTable2$maxIndex < 200000,] #IRA2
#subtable<-finalTable2[finalTable2$gene=="GPD1" & finalTable2$chromo == 5,]
points(subtable$deltaHight_Low[subtable$fluorecence=="gfp"],subtable$deltaOther[subtable$fluorecence=="gfp"],xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,1,0.5),cex=subtable$cexall[subtable$fluorecence=="gfp"])
points(subtable$deltaOther[subtable$fluorecence=="mch"],subtable$deltaHight_Low[subtable$fluorecence=="mch"],xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,1,0.5),cex=subtable$cexall[subtable$fluorecence=="mch"])

subtable<-finalTable2[finalTable2$fluorecence=="gfp" & finalTable2$maxValue > 30 & finalTable2$effectRNAPROT == "specific" ,] 

finalTable[finalTable$deltaHight_Low < -0.4,]

finalTable[finalTable$gene=="GPD1" & finalTable$chromo == 10 & finalTable$maxIndex < 200000,]

xpQTL2<-merge(xpQTL,geneInfo,by.x=1,by.y=1,all.x =T,all.y=F)

eQTL_OI<-final_eQTL[final_eQTL$chr.x == chromosome[10] & final_eQTL$maxQTLpos > 130000  & final_eQTL$maxQTLpos < 170000,] #TIF2
pQTL_OI<-xpQTL2[gsub(pattern = "chr", "", xpQTL2$chromosome) == "10" & xpQTL2$peakPosition > 130000  & xpQTL$peakPosition < 170000,] #TIF2
QTG_OI<-geneQTG[geneQTG$chr == 10 & geneQTG$posup > 130000  & geneQTG$posup < 170000,] #TIF2


