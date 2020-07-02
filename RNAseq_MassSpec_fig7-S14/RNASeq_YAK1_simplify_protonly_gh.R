

library(tidyverse)
library(DESeq2)
library(readxl)
library(readr)
library(qvalue)
library(FactoMineR)
library(sva)
library(lme4)
library(smatr)
library(reshape2)

setwd("YOUR/WORKING/DIRECTORY")
WD<-getwd()


geneAnnotation = read.table("ensemblGenes_ensembl83_160307_MOD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE) # provided
OD = read.table("OD.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE) # provided
sampleInfo = read.table("sampleInfo.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE) # provided
samples<-paste(sampleInfo$const,sampleInfo$rep,sep="")
genesINFO<-read.table(paste(WD,"/","geneINFOsc.txt", sep=""), sep="\t", head=TRUE, stringsAsFactors=FALSE, quote = "")
genesINFO2<-genesINFO[genesINFO$type=="ORF",] #gene annotation

GOI <- read.table("GOI.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)

ODs<-paste(OD$const,OD$rep,sep="")
OD600<-c()
i<-1
for (i in 1:10) {
  strain<-samples[i]
  ODt<-sampleInfo$timeSample[i]
  ODtx<-OD$time[ODs==strain]
  ODox<-OD$OD[ODs==strain]
  index<-length(ODtx[ODtx<ODt])
  time1<-ODtx[index]
  OD1<-ODox[index]
  time2<-ODtx[index+1]
  OD2<-ODox[index+1]
  OD600[i]<-OD1+((OD2-OD1)/(time2-time1)*(ODt-time1))
}
sampleInfo$OD600<-OD600

# CAREFUL: sampleInfo is initially in a different order than the samples vector; be sure to always address via the vector
# or better, reorder sampleInfo as we do below
# reorder to match "samples"

# make some batches factors (not numerics)
sampleInfo$batch <- as.factor(sampleInfo$batch)
#OD <- sampleInfo[samples,"OD600"]

GPD1<-"YDL022W"

geneCounts <- sapply(samples, function(x){
    thisDat <- read.table(paste(WD,"/","abundance/abundance",x,".tsv", sep=""), sep="\t", head=TRUE)
    thisDat[,"est_counts"]
})
rownames(geneCounts) <- read.table(paste(WD,"/","abundance/abundance","B1",".tsv", sep=""), sep="\t", head=TRUE, stringsAsFactors=FALSE)$target_id
# write out for paper/GEO
# write.table(geneCounts, file="geneCounts4paper.txt", quote=FALSE, sep="\t")

#Visual control
boxplot(log(geneCounts))
pcaRAW<-PCA(t(geneCounts),graph=F)
plot(pcaRAW,choix="ind")
GPD1rawexp<-geneCounts[rownames(geneCounts)==GPD1,]
plot(GPD1rawexp)



TINPerGene <- sapply(samples, function(x){
    theseTINs <- read.table(paste(WD,"/","abundance/pseudoalignments", x, ".tin.xls", sep=""), head=TRUE, stringsAsFactors=FALSE)
    rownames(theseTINs) <- theseTINs[,1]
    theseTINs <- theseTINs[rownames(geneCounts),]
    theseTINs[,"TIN"]
})
rownames(TINPerGene) <- rownames(geneCounts)


meanTINPerGene <- rowMeans(TINPerGene)
# effective length is annotation based, and therefore can be pulled from any one sample
effectiveLengthPerGene <- read.table(paste(WD,"/","abundance/abundance","B1",".tsv", sep=""), sep="\t", head=TRUE)[,"eff_length"]
names(effectiveLengthPerGene) <- read.table(paste(WD,"/","abundance/abundance","B1",".tsv", sep=""), sep="\t", head=TRUE)[,"target_id"]


geneCountsFiltered <- geneCounts[rowMin(round(geneCounts)) > 0 & apply(TINPerGene, 1, min) > 0 & effectiveLengthPerGene > 0,]
# down to 5755 from 6713
# the "round" is important, otherwise get Inf norm factors for three genes with "counts" very close to, but not exactly, zero

write.table(geneCountsFiltered, file="geneCountsFiltered4paper.txt", quote=FALSE, sep="\t")




############################
# SVA

# need a DESeq model first to cull expression values from
# can be any in terms of factors, just be sure to below include the ones we want
thisColDat <- sampleInfo[,colnames(sampleInfo) %in% c("const", "batch", "RIN", "consE","OD600")]
dds <- DESeqDataSetFromMatrix(countData = round(geneCountsFiltered), colData = thisColDat, design = ~ batch + const)
dds <- DESeq(dds, betaPrior=TRUE)

# SVA
dat4SVA  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat4SVA) > 1
dat4SVA  <- dat4SVA[idx, ]
mod  <- model.matrix(~ batch + const, colData(dds))
mod0 <- model.matrix(~ batch, colData(dds))
svseq <- svaseq(dat4SVA, mod, mod0, n.sv = 2)

ddssva <- dds
thisColDat_SVA <- thisColDat
thisColDat_SVA$SV1 <- svseq$sv[,1]
thisColDat_SVA$SV2 <- svseq$sv[,2]

ddssva <- DESeqDataSetFromMatrix(countData = round(geneCountsFiltered), colData = thisColDat_SVA, design = ~ batch + const) # + SV1
ddssva <- DESeq(ddssva, betaPrior=TRUE)

# what does this look like in PCA?
pdf("PCA_SVA.pdf")
vsd <- vst(ddssva, blind=TRUE)
plotPCA(vsd, intgroup="SV1")
plotPCA(vsd, intgroup="SV2")
plotPCA(vsd, intgroup="const")
plotPCA(vsd, intgroup="OD600")
plotPCA(vsd, intgroup="batch")
plotPCA(vsd, intgroup="RIN")
plotPCA(vsd, intgroup="consE")
dev.off()

ddssvaN <- DESeqDataSetFromMatrix(countData = round(geneCountsFiltered), colData = thisColDat_SVA, design = ~ batch + const)
ddssvaN <- DESeq(ddssvaN, betaPrior=TRUE)
vstC<-vst(ddssvaN,blind = F)
plotPCA(vstC, intgroup="const")
plotPCA(vstC, intgroup="batch")

expF<-assay(vstC)
expT<-assay(vsd)

plot(expT[rownames(expT)=="YDL022W",],expF[rownames(expF)=="YDL022W",],col=rep(c("black","blue"),5))
abline(a=0,b=1)

expTnromBatch<-expT
batch<-c(rep("1",4),rep("2",6))
for (i in 1:nrow(expT)) {
  M<-mean(expT[i,])
  mod<-lm(expT[i,]~as.factor(batch))
  expTnromBatch[i,]<-M+mod$residuals
}

plot(expT[rownames(expT)=="YDL022W",],expTnromBatch[rownames(expTnromBatch)=="YDL022W",],col=rep(c("black","blue"),5))
abline(a=0,b=1)

##################
# DESeq2

DESeq_results <- results(ddssva, contrast=c("const", "V", "B"))
plotMA(DESeq_results, ylim=c(-2,2))

DEresulttable <- as.data.frame(DESeq_results@rownames)
colnames(DEresulttable)<-"Gene"
DEresulttable$meanexp <- DESeq_results$baseMean
DEresulttable$logFC <- DESeq_results$log2FoldChange
DEresulttable$pval <- DESeq_results$pvalue
DEresulttable$padj <- DESeq_results$padj

DEresulttable[DEresulttable$Gene==GPD1,]
sum(DEresulttable$padj<0.05)


plot(DEresulttable$logFC,-log10(DEresulttable$padj))
points(DEresulttable$logFC[DEresulttable$Gene=="YDL022W"],-log10(DEresulttable$padj[DEresulttable$Gene=="YDL022W"]),pch=16,col="red")

plot(log2(DEresulttable$meanexp),DEresulttable$logFC,ylim=c(-1,1))
abline(h=0)
points(log2(DEresulttable$meanexp[DEresulttable$Gene=="YDL022W"]),DEresulttable$logFC[DEresulttable$Gene=="YDL022W"],pch=16,col="red")

hist(DEresulttable$logFC,breaks=100)
abline(v=0)

DEresulttable2<-merge(genesINFO2,DEresulttable,by.x = 1,by.y=1,all.x=F,all.y=T)

sign<-DEresulttable2[DEresulttable2$padj<0.05,]

write.table(DEresulttable2, file="difExp.txt", quote=FALSE, sep="\t",row.names = F)
write.table(sign, file="difExpSign.txt", quote=FALSE, sep="\t",row.names = F)

write.table(DEresulttable[DEresulttable$padj<0.05 & DEresulttable$logFC<0,1], file="RtL.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(DEresulttable[DEresulttable$padj<0.05 & DEresulttable$logFC>0,1], file="RtH.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(DEresulttable[,1], file="Rtot.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)

dplus<-density(DEresulttable$logFC[DEresulttable$padj<0.05 & DEresulttable$logFC>0])
dminus<-density(-DEresulttable$logFC[DEresulttable$padj<0.05 & DEresulttable$logFC<0])
plot(dminus,col="#34ab83AA",lty=2,lwd=2,xlim=c(0,1))
lines(dplus,col="#34ab83AA",lty=1,lwd=2)

###############################

#data from Kemmeren deleteome 2014

# Kemmerdata<-read.table("Kemmerman_deleteome_all_mutants_ex_wt_var_controls.txt",sep="\t",head = T, stringsAsFactors=FALSE, quote = "",na.strings = "NA")
# KemmerdataGenelis<-Kemmerdata[-1,1:3]
# KemmerdataYAK1<-Kemmerdata[-1,grep("yak1",x = colnames(Kemmerdata))]
# KemmerdataYAK1tot<-cbind(KemmerdataGenelis,KemmerdataYAK1)
# colnames(KemmerdataYAK1tot)<-c("reporterId","systematicName","geneSymbol","M","A","pvalue")
# write.table(KemmerdataYAK1tot,"KemmerdataYAK1.txt",sep="\t",col.names = T,quote = F,row.names = F)
KemmerdataYAK1tot<-read.table("KemmerdataYAK1.txt",sep="\t",head = T, stringsAsFactors=FALSE, quote = "",na.strings = "NA")

KemvsMe<-merge(DEresulttable,KemmerdataYAK1tot,all.x=F,all.y=F,by.x=1,by.y=2)
plot(log2(KemvsMe$meanexp),KemvsMe$A,pch=16,col=rgb(0,0,0,0.2),cex=0.5)
points(log2(KemvsMe$meanexp)[KemvsMe$geneSymbol=="GPD1"],KemvsMe$A[KemvsMe$geneSymbol=="GPD1"],pch=16,col="blue")
abline(a=0,b=1)
reg=lm(KemvsMe$A~log2(KemvsMe$meanexp))
abline(reg,col="red")
summary(reg)

plot(KemvsMe$logFC,KemvsMe$M,pch=16,col=rgb(0,0,0,0.2),cex=0.5,xlim=c(-0.8,0.8),ylim=c(-0.8,0.8))
reg=lm(KemvsMe$M~KemvsMe$logFC)
abline(a=0,b=1)
abline(h=0)
abline(v=0)
abline(reg,col="red")
summary(reg)
points(KemvsMe$logFC[KemvsMe$padj<0.05 & KemvsMe$pvalue<0.05],KemvsMe$M[KemvsMe$padj<0.05 & KemvsMe$pvalue<0.05],pch=16,col="red")
points(KemvsMe$logFC[KemvsMe$geneSymbol=="GPD1"],KemvsMe$M[KemvsMe$geneSymbol=="GPD1"],pch=16,col="blue")
c(length(KemvsMe$logFC[KemvsMe$padj<0.05 & KemvsMe$pvalue<0.05]),length(KemvsMe$logFC[KemvsMe$padj<0.05]),length(KemvsMe$logFC[KemvsMe$pvalue<0.05]))
plot(-log10(KemvsMe$padj),-log10(KemvsMe$pvalue),pch=16,col=rgb(0,0,0,0.2),cex=0.5)

KemvsMe$cat<-rep("no",nrow(KemvsMe))
KemvsMe$cat[KemvsMe$padj<0.05 & KemvsMe$pvalue<0.05 & KemvsMe$logFC<0 & KemvsMe$M<0]<-"low"
KemvsMe$cat[KemvsMe$padj<0.05 & KemvsMe$pvalue<0.05 & KemvsMe$logFC>0 & KemvsMe$M>0]<-"high"

KemvsMe2<-merge(genesINFO2,KemvsMe,by.x=1,by.y=1,all.x=F,all.y=T)

write.table(KemvsMe2, file="Kem_vs_Brion.txt", quote=FALSE, sep="\t",row.names = F)
write.table(KemvsMe2[KemvsMe2$cat=="low",1], file="kemvsL.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(KemvsMe2[KemvsMe2$cat=="high",1], file="kemvsH.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)

#=================
#protdata

protdata<-read.table(paste(WD,"/","allProt.txt", sep=""), sep="\t", head=TRUE, stringsAsFactors=FALSE, quote = "",na.strings = "NA")
grep("std",colnames(protdata))
protdatanoSTD<-protdata[,-grep("std",colnames(protdata))]
protdatanoSTD$name[protdatanoSTD$name=="ARG5"]<-"ARG5,6"
protdatanoSTD$name[protdatanoSTD$name=="ADE5"]<-"ADE5,7"
protdatanoSTD$name[protdatanoSTD$name=="AIM1"]<-"BOL3"
protdatanoSTD$name[protdatanoSTD$name=="NAT10"]<-"KRE33"
protdatanoSTD$name[protdatanoSTD$IDyeast=="SYV_YEAST"]<-"VAS1"
protdatanoSTD$name[protdatanoSTD$IDyeast=="INV2_YEAST"]<-"SUC2"
protdatanoSTD$name[protdatanoSTD$IDyeast=="SYG_YEAST [2]"]<-"GRS1"
protdatanoSTD$name[protdatanoSTD$IDyeast=="DNLI1_YEAST"]<-"CDC9"
protdatanoSTD$name[protdatanoSTD$IDyeast=="SYH_YEAST"]<-"HTS1"
protdatanoSTD$name[protdatanoSTD$IDyeast=="BLH1_YEAST"]<-"LAP3"
protdatanoSTD$name[protdatanoSTD$IDyeast=="IF5_YEAST"]<-"TIF5"
protdatanoSTD$name[protdatanoSTD$IDyeast=="LEU1_YEAST [2]"]<-"LEU4"
protdatanoSTD$name[protdatanoSTD$IDyeast=="SYA_YEAST"]<-"ALA1"
pt<-protdatanoSTD[protdatanoSTD$name%in%genesINFO2$name,]
pF<-protdatanoSTD[!(protdatanoSTD$name %in% genesINFO2$name),]
pFt<-pF[(pF$name %in% genesINFO2$sysname),]
pFF<-pF[!(pF$name %in% genesINFO2$sysname),]
testgene<-c()
keep<-c()
for (i in 1:nrow(pFF)){
  tempsinfo<-strsplit(pFF$IDyeast[i],split = "_")[[1]]
  if (tempsinfo[1]%in%genesINFO2$name & !(tempsinfo[1]%in%pt$name) & tempsinfo[2]!="YEAST-DECOY") {
    keep[length(keep)+1]<-i
    pFF$name[i]<-tempsinfo[1]
  } else {
    testgene[length(testgene)+1]<-pFF$IDyeast[i]
  }
}
pFFt<-pFF[keep,]
pFFF<-pFF[-keep,]
protdataMt<-merge(genesINFO2,pt,by.x="name",by.y="name",all.x=F,all.y=F)
protdataMt<-protdataMt[,c(2,1,3:ncol(protdataMt))]
protdataMf<-merge(genesINFO2,pFt,by.x="sysname",by.y="name",all.x=F,all.y=F)
protdataMft<-merge(genesINFO2,pFFt,by.x="name",by.y="name",all.x=F,all.y=F)
protdataMft<-protdataMft[,c(2,1,3:ncol(protdataMft))]
protdataM<-rbind(protdataMt,protdataMf,protdataMft)
protdataM$adjpV<-as.numeric(sub(pattern = "< ",replacement = "",x = protdataM$adjpV))
protdataN<-as.matrix(protdataM[,c((ncol(protdataM)-9):ncol(protdataM))])
rownames(protdataN)<-protdataM$sysname
protdataM$Average_prot<-rowMeans(protdataN,na.rm = T)

protdataM$BHadj<-p.adjust(protdataM$adjpV, method = "BH")
hist(-log10(protdataM$BHadj))
sum(protdataM$BHadj<0.05)
sum(protdataM$BHadj<0.01)
protdataM[protdataM$sysname==GPD1,]

write.table(protdataM, file="difProt.txt", quote=FALSE, sep="\t",row.names = F)
write.table(protdataM[protdataM$BHadj<0.05,], file="difProtSign.txt", quote=FALSE, sep="\t",row.names = F)

write.table(protdataM[protdataM$BHadj<0.05 & protdataM$effectV<0,1], file="PtL.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(protdataM[protdataM$BHadj<0.05 & protdataM$effectV>0,1], file="PtH.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(protdataM[,1], file="Ptot.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)

dplus<-density(protdataM$effectV[protdataM$BHadj<0.05 & protdataM$effectV>0])
dminus<-density(-protdataM$effectV[protdataM$BHadj<0.05 & protdataM$effectV<0])
plot(dplus,col="red",lty=1,lwd=2,xlim=c(0,1))
lines(dminus,col="red",lty=2,lwd=2)
dall<-density(protdataM$effectV)
plot(dall,col="red",lty=1,lwd=2,xlim=c(-1,1))
abline(v=0)
mean(protdataM$effectV)
t.test(protdataM$effectV,alternative = "greater")
t.test(protdataM$effectV,alternative = "less")
wilcox.test(protdataM$effectV,alternative="less")

#test of degradation (n-en dregon effect library (Kats et al 2018))
# deg<-read.table("PSIdegron.txt",header = T,sep="\t",quote = "",stringsAsFactors = F,na.strings = "")
# deg2<-deg[,c(1,48)]
# deg2<-deg2[!(is.na(deg2$wild.type)),]
# degvsprot<-merge(protdataM,deg2,by.x=1,by.y=1,all.x=F,all.y=T)
# degvsprot$cat<-rep("Q",nrow(degvsprot))
# degvsprot$cat[is.na(degvsprot$name)]<-"nQ"
# t.test(degvsprot$wild.type[degvsprot$cat=="nQ"],degvsprot$wild.type[degvsprot$cat=="Q"])
# boxplot(degvsprot$wild.type[degvsprot$cat=="nQ"],degvsprot$wild.type[degvsprot$cat=="Q"])
# plot(degvsprot$wild.type,degvsprot$effectV)
# reg<-lm(degvsprot$effectV~degvsprot$wild.type)
# t.test(degvsprot$wild.type[degvsprot$effectV>0 & degvsprot$BHadj<0.05],degvsprot$wild.type[degvsprot$effectV<0 & degvsprot$BHadj<0.05])

valpos<-c()
valneg<-c()
for (i in 0:100) {
  valpos[i+1]<-sum(protdataM$effectV > i/100 & protdataM$effectV < (i+4)/100)
  valneg[i+1]<-sum(protdataM$effectV < -i/100 & protdataM$effectV > -(i+4)/100)
}
plot((0:100)/100,valpos,typ="l")
points((0:100)/100,valneg,typ="l",col="red")
plot((0:100)/100,valpos-valneg,typ="l")
abline(h=0)


#===========================================
#compare


dif_redVP<-merge(protdataM,DEresulttable,by.x=1,by.y=1,all.x=F,all.y=F)
colnames(dif_redVP)
colnames(dif_redVP)[c(1:10,33:36,31,19,18,32)]
dif_redVP<-dif_redVP[,c(1:10,33:36,31,19,18,32)]
colnames(dif_redVP)<-c(colnames(genesINFO),"R_A","R_FC","R_pv","R_apv","P_A","P_FC","P_pv","P_apv")
dif_redVP$R_A<-log2(dif_redVP$R_A)
dif_redVP[dif_redVP$sysname==GPD1,]

plot(dif_redVP$R_FC,dif_redVP$P_FC,cex=0.5,xlim=c(-1,1),ylim=c(-1,1))
abline(a=0,b=1)
abline(h=0)
abline(v=0)
reg<-lm(dif_redVP$P_FC~dif_redVP$R_FC)
abline(reg,col="red")
summary(reg)

dif_redVP$P_aC<- -log10(dif_redVP$P_apv)*sign(dif_redVP$P_FC)
dif_redVP$R_aC<- -log10(dif_redVP$R_apv)*sign(dif_redVP$R_FC)
dif_redVP$P_pC<- -log10(dif_redVP$P_pv)*sign(dif_redVP$P_FC)
dif_redVP$R_pC<- -log10(dif_redVP$R_pv)*sign(dif_redVP$R_FC)
thre1 <- -log10(0.05)
thre2 <- -log10(0.05)

dif_redVP$cat<-rep("no",nrow(dif_redVP))
dif_redVP$cat[dif_redVP$R_aC < -thre1 & dif_redVP$P_pC > -thre2]<-"RL"
dif_redVP$cat[dif_redVP$R_aC > thre1 & dif_redVP$P_pC < thre2]<-"RH"
dif_redVP$cat[dif_redVP$R_pC > -thre2 & dif_redVP$P_aC < -thre1]<-"PL"
dif_redVP$cat[dif_redVP$R_pC < thre2 & dif_redVP$P_aC > thre1]<-"PH"
dif_redVP$cat[dif_redVP$R_aC < -thre1 & dif_redVP$P_aC < -thre1]<-"RPL"
dif_redVP$cat[dif_redVP$R_aC > thre1 & dif_redVP$P_aC > thre1]<-"RPH"
summary(as.factor(dif_redVP$cat))

cate<-c("no","RL","RH","PL","PH","RPL","RPH")
colcat<-c("#AAAAAA22","#ff0000AA","#ff0000AA","#34ab83AA","#34ab83AA","#000000AA","#000000AA")
colcat2<-c("#AAAAAAAA","#df0000ff","#df0000ff","#249b73ff","#249b73ff","#000000ff","#000000ff")
cexcat2<-c(0.3,0.7,0.7,0.7,0.7,0.7,0.7)

write.table(dif_redVP[dif_redVP$cat=="PL",1], file="PL.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(dif_redVP[dif_redVP$cat=="PH",1], file="PH.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(dif_redVP[dif_redVP$cat=="RL",1], file="RL.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(dif_redVP[dif_redVP$cat=="RH",1], file="RH.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(dif_redVP[dif_redVP$cat=="RPL",1], file="RPL.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(dif_redVP[dif_redVP$cat=="RPH",1], file="RPH.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)
write.table(dif_redVP[,1], file="allSynt.txt", quote=FALSE, sep="\t",row.names = F,col.names = F)

svg("RNAseq_VS_prot_Yak1.svg",width = 5.7,height = 6)
plot(dif_redVP$R_FC,dif_redVP$P_FC,cex=0.5,typ="n",xlim=c(-1,1),ylim=c(-1,1),pch=16,col=rgb(0,0,0,0.2))
for (i in 1:7) {
  points(dif_redVP$R_FC[dif_redVP$cat==cate[i]],dif_redVP$P_FC[dif_redVP$cat==cate[i]],pch=16,col=colcat[i],cex=cexcat2[i])
  points(dif_redVP$R_FC[dif_redVP$cat==cate[i]],dif_redVP$P_FC[dif_redVP$cat==cate[i]],pch=1,col=colcat2[i],cex=cexcat2[i])
}
abline(a=0,b=1)
abline(h=0)
abline(v=0)
points(dif_redVP$R_FC[dif_redVP$sysname==GPD1],dif_redVP$P_FC[dif_redVP$sysname==GPD1],pch=16,col="#1177FFFF",cex=1.2)
points(dif_redVP$R_FC[dif_redVP$sysname==GPD1],dif_redVP$P_FC[dif_redVP$sysname==GPD1],pch=1,col="blue",cex=1.2)
reg<-lm(dif_redVP$P_FC~dif_redVP$R_FC)
#abline(reg,col="red")
summary(reg)
dev.off()

write.table(dif_redVP, file="difCOMPA.txt", quote=FALSE, sep="\t",row.names = F)

#GPD1 data 
#(prot/rna)

GOItable<-read.table(file = "genetoplot.txt",header = T,sep="\t",quote = "",stringsAsFactors = F)

plot(dif_redVP$R_FC,dif_redVP$P_FC,cex=0.5,typ="n",xlim=c(-1,1),ylim=c(-1,1),pch=16,col=rgb(0,0,0,0.2))
for (i in 1:7) {
  points(dif_redVP$R_FC[dif_redVP$cat==cate[i]],dif_redVP$P_FC[dif_redVP$cat==cate[i]],pch=16,col=colcat[i],cex=cexcat2[i])
  points(dif_redVP$R_FC[dif_redVP$cat==cate[i]],dif_redVP$P_FC[dif_redVP$cat==cate[i]],pch=1,col=colcat2[i],cex=cexcat2[i])
}
abline(a=0,b=1)
abline(h=0)
abline(v=0)
points(dif_redVP$R_FC[dif_redVP$sysname%in%GOItable$gene],dif_redVP$P_FC[dif_redVP$sysname%in%GOItable$gene],pch=16,col="#1177FFFF",cex=1.2)
points(dif_redVP$R_FC[dif_redVP$sysname%in%GOItable$gene],dif_redVP$P_FC[dif_redVP$sysname%in%GOItable$gene],pch=1,col="blue",cex=1.2)
text(dif_redVP$R_FC[dif_redVP$sysname%in%GOItable$gene],dif_redVP$P_FC[dif_redVP$sysname%in%GOItable$gene],labels = dif_redVP$name[dif_redVP$sysname%in%GOItable$gene],col="black",cex=0.7)


goi_info<-dif_redVP[dif_redVP$sysname%in%GOItable$gene,]
goi_info2<-goi_info[,c(1,2,9,10,12,14,16,18,23)]

GOI<-GOItable[1,1]
const<-c(rep("black",5),rep("blue",5))
rnaGOI<-expTnromBatch[rownames(expTnromBatch)==GOI,c(1,3,5,7,9,2,4,6,8,10)]
protGOI<-protdataN[rownames(protdataN)==GOI,]
plot(rnaGOI,protGOI,col=c(rep("black",5),rep("blue",5)))
abline(a=0,b=1)
expr<-c(rnaGOI-mean(rnaGOI),protGOI-mean(protGOI))
const<-as.factor(rep(c(rep("B",5),rep("V",5)),2))
type<-as.factor(c(rep("RNA",10),rep("PROT",10)))
sample<-as.factor(names(expr))
boxplot(expr~paste(type,const),col=c("#34ab83FF","#34ab83FF","red","red"))

finEffect<-lm(expr~const*type)
summary(finEffect)

finEffect2<-lmer(expr~const*type +(1|sample) )
finEffect3<-lmer(expr~const+type +(1|sample) )
anova(finEffect3,finEffect2)

#=========
#cito-translation

citotrans<-read.table(file = "citotrans.txt",sep="\t",quote = "",stringsAsFactors = F)
citotransexp<-dif_redVP[dif_redVP$sysname %in% citotrans[,3],]
sum(citotransexp$R_aC<log10(0.05))
sum(dif_redVP$R_aC<log10(0.05))
citotransexpR<-DEresulttable2[DEresulttable2$sysname %in% citotrans[,3],]
sum(citotransexpR$padj<0.05 & citotransexpR$logFC<0)

svg("RNAseq_VS_prot_Yak1_translation.svg",width = 5.7,height = 6)
plot(dif_redVP$R_FC,dif_redVP$P_FC,cex=0.5,typ="n",xlim=c(-1,1),ylim=c(-1,1),pch=16,col=rgb(0,0,0,0.2))
points(dif_redVP$R_FC,dif_redVP$P_FC,pch=16,col=colcat[1],cex=cexcat2[1])
points(dif_redVP$R_FC,dif_redVP$P_FC,pch=1,col=colcat2[1],cex=cexcat2[1])
points(citotransexp$R_FC[citotransexp$R_aC>log10(0.05)],citotransexp$P_FC[citotransexp$R_aC>log10(0.05)],pch=16,col="#d28c8cff",cex=0.85)
points(citotransexp$R_FC[citotransexp$R_aC<log10(0.05)],citotransexp$P_FC[citotransexp$R_aC<log10(0.05)],pch=16,col="#e23f3fff",cex=0.85)
points(citotransexp$R_FC[citotransexp$R_aC>log10(0.05)],citotransexp$P_FC[citotransexp$R_aC>log10(0.05)],pch=1,col="#e23f3fff",cex=0.85)
points(citotransexp$R_FC[citotransexp$R_aC<log10(0.05)],citotransexp$P_FC[citotransexp$R_aC<log10(0.05)],pch=1,col="#ae1a1aff",cex=0.85)
abline(a=0,b=1)
abline(h=0)
abline(v=0)
dev.off()

phyper(q=32, m=113, n=2577-113, k=170, lower.tail = F)
