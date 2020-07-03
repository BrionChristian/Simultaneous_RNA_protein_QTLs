library(plyr)
library(reshape2)


alignmentDir <- "C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/190208-AllQTLprocess/"
setwd(alignmentDir)

#library("VariantAnnotation")
source(paste(alignmentDir,"/Scripts/gTest.R",sep=""))
source(paste(alignmentDir,"/Scripts/x_qtl_seq_functions_170831.R",sep=""))
source(paste(alignmentDir,"/Scripts/mp_JB_170901.R",sep=""))
source(paste(alignmentDir,"/Scripts/peaksFromVector.R",sep=""))


#set all parameter for politing and analysis
i = 1
sepBetweenChr <- 1e5
trimFromEnd = 15e3
obsMin <- 10
LoessSpan = 0.1
AFThres = 0.09653124 # same as in Albert 2014
multiThres = 3 # LOD threshold for multipool, if run with N=1000: 4.5
withMultipool <- TRUE

#get external table: variants, chromosomes, genes
SNPs <- read.table(paste(alignmentDir,"/Scripts/SNPs_Maggie_170809_BY_positions.txt",sep=""), stringsAsFactors=FALSE, head=FALSE)
# see comments above. As of 8/31/17, the SNPs seem not to be fully filtered, and are out of sorting order
for (thisChr in unique(SNPs[,1])){SNPs[SNPs[,1] == thisChr, 2] <- sort(SNPs[SNPs[,1] == thisChr, 2])}
SNPs <- rbind(SNPs[SNPs[,1] == "chrI",], SNPs[SNPs[,1] == "chrII",], SNPs[SNPs[,1] == "chrIII",], SNPs[SNPs[,1] == "chrIV",], SNPs[SNPs[,1] == "chrV",], SNPs[SNPs[,1] == "chrVI",], SNPs[SNPs[,1] == "chrVII",], SNPs[SNPs[,1] == "chrVIII",], SNPs[SNPs[,1] == "chrIX",], SNPs[SNPs[,1] == "chrX",], SNPs[SNPs[,1] == "chrXI",], SNPs[SNPs[,1] == "chrXII",], SNPs[SNPs[,1] == "chrXIII",], SNPs[SNPs[,1] == "chrXIV",], SNPs[SNPs[,1] == "chrXV",], SNPs[SNPs[,1] == "chrXVI",])
sgd_table <- paste(alignmentDir,"/Scripts/sacCer3ChromLenghts.txt",sep="")
chrInfo = read.table(sgd_table, stringsAsFactors=FALSE, sep="\t", header=F)
geneInfo = read.table(paste(alignmentDir,"/Scripts/ensemblGenes_ensembl83_160307_MOD.txt",sep=""), stringsAsFactors=FALSE, sep="\t", header=TRUE)
rownames(geneInfo) <- geneInfo[,"geneID"]
allNames <- geneInfo[, "geneName"]
names(allNames) <- geneInfo[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames
chromosome<-unique(SNPs[,1])

#Albert et al data: pQTL and eQTL
xpQTL <- read.table("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/Biblio/Frank_eQTL_pQTL/xpQTLs.txt",header = T,sep = "\t", quote = "", na.strings = "NA",stringsAsFactors = F)
final_eQTL <- read.table("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/eQTLFrank/final_eQTL.txt",header = T,sep = "\t", quote = "", na.strings = "NA",stringsAsFactors = F,comment.char = "")
geneQTG <- read.table("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/eQTLFrank/genes_RMBY_variant.txt",header = T,sep = "\t", quote = "", na.strings = "",stringsAsFactors = F,comment.char = "")
geneQTG<-geneQTG[!(is.na(geneQTG$chr)),]

#table with summury of all the QTL experiment
experimentFile <- read.table(paste(alignmentDir,"AllPoprecap.txt",sep=""), stringsAsFactors=FALSE, head=TRUE, na.strings = "")
experimentFile$tube<-paste(experimentFile$gene,experimentFile$crispr,experimentFile$rep,experimentFile$tech,sep="_")
experimentFile$exp<-paste(experimentFile$gene,substr(experimentFile$crispr,1,1),sep="_")
experimentFile$globrep<-paste(substr(experimentFile$crispr,2,2),experimentFile$rep,experimentFile$tech,sep="_")
experimentFile2<-experimentFile[experimentFile$crispr!='i' & experimentFile$tube!= "GPD1_a_3_2" & experimentFile$tube!="RPS10A_a_2_1",]
experimentFile3<-experimentFile[experimentFile$gene=='NA',] #experimentFile$crispr!='i' &
experimentFile2<-rbind(experimentFile2,experimentFile3)
levels(as.factor(experimentFile2$tube))
levels(as.factor(experimentFile2$exp))
levels(as.factor(experimentFile2$globrep))

initialtable<-read.table(paste(alignmentDir,"All_QTL_fused.txt",sep=""), stringsAsFactors=FALSE, head=TRUE)
#====================================================
# separate analysis of QTL

summary(as.factor(initialtable$effectRNAPROT))
summary(as.factor(initialtable$fluorecence))
summary(as.factor(initialtable$effectRNAPROT[initialtable$crispr!="i"]))

finalTable<-initialtable[initialtable$crispr!="i" & !(initialtable$tube%in%c("RPS10A_a_2_1","TDH3_a_0_1")),]
#finalTable<-initialtable
summary(as.factor(finalTable$effectRNAPROT))
summary(as.factor(finalTable$fluorecence))
summary(as.factor(paste(finalTable$fluorecence,finalTable$effectRNAPROT)))
summary(as.factor(paste(finalTable$fluorecence,finalTable$gene)))

cexall<-sapply(finalTable$maxValue,FUN =  function(x) {max(0.4,min(3,x/20))})
finalTable2<-finalTable
finalTable2$cexall<-cexall

finalTable45<-finalTable2[finalTable2$maxValue>4.5 & abs(finalTable2$deltaHight_Low)>0.05,]
summary(as.factor(paste(finalTable45$fluorecence,finalTable45$effectRNAPROT)))
finalTable2<-finalTable45

finalTable2<-finalTable
finalTable2$cexall<-cexall
plot(finalTable2$deltaHight_Low[finalTable2$fluorecence=="gfp"],finalTable2$deltaOther[finalTable2$fluorecence=="gfp"],
     xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0.8,0.5),cex=finalTable2$cexall[finalTable2$fluorecence=="gfp"],
     xlab="AF(BY/all) difference (high-low) - GFP gate", ylab="AF(BY/all) difference (high-low) - mCherry gate",
     main="QTL effect: protein vs mRNA")
points(finalTable2$deltaOther[finalTable2$fluorecence=="mch"],finalTable2$deltaHight_Low[finalTable2$fluorecence=="mch"],xlim=c(-1,1),ylim=c(-1,1),col="red",cex=finalTable2$cexall[finalTable2$fluorecence=="mch"])
abline(v=0)
abline(h=0)
#subtable<-finalTable2[finalTable2$gene=="GPD1" & finalTable2$chromo == 10 & finalTable2$maxIndex > 100000  & finalTable2$maxIndex < 200000,]
#subtable<-finalTable2[finalTable2$chromo == 10 & finalTable2$maxIndex > 100000  & finalTable2$maxIndex < 200000,] #TIF2
#subtable<-finalTable2[finalTable2$chromo == 12 & finalTable2$maxIndex > 620000  & finalTable2$maxIndex < 680000,] #HAP1
#subtable<-finalTable2[finalTable2$chromo == 14 & finalTable2$maxIndex > 440000  & finalTable2$maxIndex < 500000,] #MKT1
#subtable<-finalTable2[finalTable2$chromo == 4 & finalTable2$maxIndex > 1000000  & finalTable2$maxIndex < 1100000,] #GCN2
#subtable<-finalTable2[finalTable2$chromo == 15 & finalTable2$maxIndex > 140000  & finalTable2$maxIndex < 200000,] #IRA2
#subtable<-finalTable2[finalTable2$gene=="GPD1" & finalTable2$chromo == 5,]
#points(subtable$deltaHight_Low[subtable$fluorecence=="gfp"],subtable$deltaOther[subtable$fluorecence=="gfp"],xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,1,0.5),cex=subtable$cexall[subtable$fluorecence=="gfp"])
#points(subtable$deltaOther[subtable$fluorecence=="mch"],subtable$deltaHight_Low[subtable$fluorecence=="mch"],xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,1,0.5),cex=subtable$cexall[subtable$fluorecence=="mch"])

summary(as.factor(finalTable45$numbAcross[finalTable45$gene!="RPS10A"]))

subtable<-finalTable2[finalTable2$fluorecence=="gfp" & finalTable2$maxValue > 30 & finalTable2$effectRNAPROT == "specific" ,] 

finalTable[finalTable$deltaHight_Low < -0.4,]
finalTable[abs(finalTable$deltaHight_Low) < 0.04,]
finalTable<-finalTable[abs(finalTable$deltaHight_Low) > 0.04,]

trsh<-4.5
sum(sign(finalTable$deltaHight_Low[finalTable$maxValue>trsh])==sign(finalTable$deltaOther[finalTable$maxValue>trsh]))/nrow(finalTable[finalTable$maxValue>trsh,])
x<-1:20
y<-1:20
for (i in x) {
  y[i]<-sum(sign(finalTable$deltaHight_Low[finalTable$maxValue>i])==sign(finalTable$deltaOther[finalTable$maxValue>i]))/nrow(finalTable[finalTable$maxValue>i,])
}
plot(x,y,ylim=c(0.5,1),type="l",xlab="lod threshold",ylab = "percent of loci with similare AF direction")

finalTable[finalTable$chromo == 10 & finalTable$maxIndex < 200000,]

#================================================================#
#fuse mCherry and GFP QTLs

finalcumul<-finalTable[0,c(1,2,9,18)]
type<-c()
GFP_LOD<-c()
GFP_effect<-c()
GFP_maxPOS<-c()
GFP_leftPOS<-c()
GFP_rightPOS<-c()
mCH_LOD<-c()
mCH_effect<-c()
mCH_maxPOS<-c()
mCH_leftPOS<-c()
mCH_rightPOS<-c()
temp_mCH<-finalTable[finalTable$fluorecence=="mch",]
temp_gfp<-finalTable[finalTable$fluorecence=="gfp",]
#first fuse the GFP
for (i in 1:nrow(temp_gfp)) {
  finalcumul<-rbind(finalcumul,temp_gfp[i,c(1,2,9,18)])
  GFP_LOD[i]<-temp_gfp$maxValue[i]
  GFP_effect[i]<-temp_gfp$deltaHight_Low[i]
  GFP_maxPOS[i]<-temp_gfp$maxIndex[i]
  GFP_leftPOS[i]<-temp_gfp$leftIndex[i]
  GFP_rightPOS[i]<-temp_gfp$rightIndex[i]
  I_mCH<-(1:nrow(temp_mCH))[temp_mCH$orf==temp_gfp$orf[i] & temp_mCH$chromo==temp_gfp$chromo[i] & temp_mCH$rightIndex>min(temp_gfp$leftIndex[i],temp_gfp$maxIndex[i]) & temp_mCH$leftIndex<max(temp_gfp$rightIndex[i],temp_gfp$maxIndex[i])]
  if (length(I_mCH)==1) {
    mCH_LOD[i]<-temp_mCH$maxValue[I_mCH]
    mCH_effect[i]<-temp_mCH$deltaHight_Low[I_mCH]
    mCH_maxPOS[i]<-temp_mCH$maxIndex[I_mCH]
    mCH_leftPOS[i]<-temp_mCH$leftIndex[I_mCH]
    mCH_rightPOS[i]<-temp_mCH$rightIndex[I_mCH]
    if (sign(mCH_effect[i])==sign(GFP_effect[i])) {
      type[i]<-"together"
    } else {
      type[i]<-"oposite"
    }
    temp_mCH<-temp_mCH[-I_mCH,]
  } else if (length(I_mCH)>1) {
    print(temp_mCH[I_mCH,])
    break()
  } else {
    mCH_LOD[i]<-temp_gfp$LOD_Other[i]
    mCH_effect[i]<-temp_gfp$deltaOther[i]
    mCH_maxPOS[i]<-temp_gfp$maxIndex[i]
    mCH_leftPOS[i]<-temp_gfp$leftIndex[i]
    mCH_rightPOS[i]<-temp_gfp$rightIndex[i]
    type[i]<-paste("gfp",sep="")
  }
}
tempdata<-data.frame(type,GFP_LOD,GFP_effect,GFP_maxPOS,GFP_leftPOS,GFP_rightPOS
                     ,mCH_LOD,mCH_effect,mCH_maxPOS,mCH_leftPOS,mCH_rightPOS)
finalcumul2<-cbind(finalcumul,tempdata)
n<-i
#then fuse the mCherry
for (i in 1:nrow(temp_mCH)) {
  finalcumul<-rbind(finalcumul,temp_mCH[i,c(1,2,9,18)])
  mCH_LOD[n+i]<-temp_mCH$maxValue[i]
  mCH_effect[n+i]<-temp_mCH$deltaHight_Low[i]
  mCH_maxPOS[n+i]<-temp_mCH$maxIndex[i]
  mCH_leftPOS[n+i]<-temp_mCH$leftIndex[i]
  mCH_rightPOS[n+i]<-temp_mCH$rightIndex[i]
  GFP_LOD[n+i]<-temp_mCH$LOD_Other[i]
  GFP_effect[n+i]<-temp_mCH$deltaOther[i]
  GFP_maxPOS[n+i]<-temp_mCH$maxIndex[i]
  GFP_leftPOS[n+i]<-temp_mCH$leftIndex[i]
  GFP_rightPOS[n+i]<-temp_mCH$rightIndex[i]
  if (temp_mCH$effectRNAPROT[i]!="specific") {
    type[n+i]<-paste(temp_mCH$effectRNAPROT[i],"_mch",sep="") #"mch"# #cases when the mCherry overlap towar GFP but not the other way around
  } else {
    type[n+i]<-"mch"
  }
}
tempdata<-data.frame(type,GFP_LOD,GFP_effect,GFP_maxPOS,GFP_leftPOS,GFP_rightPOS
                     ,mCH_LOD,mCH_effect,mCH_maxPOS,mCH_leftPOS,mCH_rightPOS)
finalcumul2<-cbind(finalcumul,tempdata)
finalcumul3<-finalcumul2[order(finalcumul2$orf,finalcumul2$chromo,finalcumul2$GFP_maxPOS),]
summary(as.factor(finalcumul3$type))
summary(as.factor(finalcumul3$type[finalcumul3$GFP_LOD>4.5 | finalcumul3$mCH_LOD>4.5]))
summary(as.factor(finalcumul3$type[finalcumul3$GFP_LOD>10 | finalcumul3$mCH_LOD>10]))

#manual curration
finalcumul3[finalcumul3$type %in% c("oposite_mch","together_mch"),]
finalcumul3$type[finalcumul3$type %in% c("oposite_mch","together_mch") & finalcumul3$gene=="GPD1"]<-"mch"
finalcumul3$type[finalcumul3$type %in% c("oposite_mch","together_mch") & finalcumul3$gene=="ARO8"]<-"together"
finalcumul3$type[finalcumul3$type %in% c("oposite_mch","together_mch") & finalcumul3$gene=="CYC1"]<-"oposite"
finalcumul3$type[finalcumul3$type %in% c("oposite_mch","together_mch") & finalcumul3$gene=="MTD1"]<-c("together","mch")
finalcumul3$type[finalcumul3$type %in% c("oposite_mch","together_mch") & finalcumul3$gene=="RPS10A"]<-"together"
finalcumul3<-finalcumul3[!(finalcumul3$type == "gfp" & finalcumul3$gene=="MTD1" & finalcumul3$chromo==8),]

#finalcumul3$type[finalcumul3$type=="oposite_mch"]<-"mch"
#finalcumul3$type[finalcumul3$type=="together_mch"]<-"mch"
summary(as.factor(finalcumul3$type[finalcumul3$GFP_LOD>4.5 | finalcumul3$mCH_LOD>4.5]))
summary(as.factor(finalcumul3$type[finalcumul3$GFP_LOD>4.5 | finalcumul3$mCH_LOD>4.5]))["together"]/length(finalcumul3$type[finalcumul3$GFP_LOD>4.5 | finalcumul3$mCH_LOD>4.5])

finalcumul3[grep(pattern = "_",finalcumul3$type),]

x<-1:20
y<-1:20
for (i in x) {
  y[i]<-summary(as.factor(finalcumul3$type[finalcumul3$GFP_LOD>i | finalcumul3$mCH_LOD>i]))["together"]/length(finalcumul3$type[finalcumul3$GFP_LOD>i | finalcumul3$mCH_LOD>i])
}
plot(x,y,ylim=c(0,0.5),type="l",xlab="lod threshold",ylab = "percent of prot&RNA loci")

write.table(finalcumul3,file = "QTLfuseFinal.txt",sep="\t",col.names = T,row.names = F,quote = F)

cexallF<-c()
colF<-c()
for (i in 1:nrow(finalcumul3)) {
  cexallF[i]<-max(0.4,min(3,max(finalcumul3$GFP_LOD[i],finalcumul3$mCH_LOD[i])/20))
  if (finalcumul3$type[i]=="gfp") {
    colF[i]<-rgb(0,0.8,0.5)
  } else if (finalcumul3$type[i]=="mch") {
    colF[i]<-rgb(1,0,0)
  } else if (finalcumul3$type[i]=="oposite") {
    colF[i]<-rgb(0,0,1)
  } else if (finalcumul3$type[i]=="together") {
    colF[i]<-rgb(0,0,0)
  } else {
    colF[i]<-rgb(0.8,0.8,0.2,0.9)
  }
}
finalcumul3$cexallF<-cexallF
finalcumul3$colF<-colF

finalcumul45<-finalcumul3[finalcumul3$GFP_LOD>4.5 | finalcumul3$mCH_LOD>4.5,]

genesummary<-dcast(finalcumul45,gene~type,length)
genesummary<-genesummary[c(1,2,5,6,10,4,7,9,3,8),c(1,5,2,3,4)]
barplot(t(as.matrix(genesummary[,-1])),col=c("#404040","#33ab82ff","red","blue"),ylim=c(0,20),names.arg = genesummary$gene)

finalcumul45_14false<-finalcumul45[(finalcumul45$chromo==14 & abs(finalcumul45$mCH_maxPOS-450000)<100000),]
finalcumul45_14<-finalcumul45[!(finalcumul45$chromo==14 & abs(finalcumul45$mCH_maxPOS-450000)<100000),]

genesummary<-dcast(finalcumul45_14,gene~type,length)
genesummary<-genesummary[c(1,2,5,6,10,4,7,9,3,8),c(1,5,2,3,4)]
svg("allgeneCum_summary.svg",width = 6,height = 4)
barplot(t(as.matrix(genesummary[,-1])),col=c("#404040","#33ab82ff","red","blue"),ylim=c(0,20),names.arg = genesummary$gene)
dev.off()
summary(as.factor(finalcumul45$type))

summary(as.factor(finalcumul45_14$type))


svg("prot-mRNA-QTLeffect.svg",width = 6,height = 6)

plot(finalcumul45_14$GFP_effect,finalcumul45_14$mCH_effect,
     xlim=c(-1,1),ylim=c(-1,1),col=finalcumul45_14$colF,cex=finalcumul45_14$cexallF,pch=1,
     xlab="AF(BY/all) difference (high-low) - GFP gate", ylab="AF(BY/all) difference (high-low) - mCherry gate",
     main="QTL effect: protein vs mRNA")
abline(v=0)
abline(h=0)
abline(a=0,b=1)
reg<-lm(finalcumul45_14$mCH_effect~finalcumul45_14$GFP_effect)
abline(reg,col="grey",lty=2)
summary(reg)
cor(finalcumul45_14$GFP_effect,finalcumul45_14$mCH_effect,method = "p")
points(finalcumul45_14false$GFP_effect,finalcumul45_14false$mCH_effect,
     xlim=c(-1,1),ylim=c(-1,1),col="grey",cex=finalcumul45_14false$cexallF,pch=1)
sum(abs(finalcumul45_14$GFP_effect[finalcumul45_14$type=="together"]))/sum(abs(finalcumul45_14$GFP_effect[finalcumul45_14$type!="mch"]))
sum(abs(finalcumul45_14$mCH_effect[finalcumul45_14$type=="together"]))/sum(abs(finalcumul45_14$mCH_effect[finalcumul45_14$type!="gfp"]))
sum(sign(finalcumul45_14$mCH_effect)==sign(finalcumul45_14$GFP_effect))/nrow(finalcumul45_14)

dev.off()

genesummary2<-genesummary
sim_GFPtomCH<-c()
sim_mCHtoGFP<-c()
sign_compa<-c()
for (i in 1:nrow(genesummary)) {
  tempcum<-finalcumul45_14[finalcumul45_14$gene==genesummary$gene[i],]
  sim_GFPtomCH[i]<-sum(abs(tempcum$GFP_effect[tempcum$type=="together"]))/sum(abs(tempcum$GFP_effect[tempcum$type!="mch"]))
  sim_mCHtoGFP[i]<-sum(abs(tempcum$mCH_effect[tempcum$type=="together"]))/sum(abs(tempcum$mCH_effect[tempcum$type!="gfp"]))
  sign_compa[i]<-sum(sign(tempcum$mCH_effect)==sign(tempcum$GFP_effect))/nrow(tempcum)
  print(paste(sum(sign(tempcum$mCH_effect)==sign(tempcum$GFP_effect)),"/",nrow(tempcum)))
}
genesummary2$sim_GFPtomCH<-sim_GFPtomCH
genesummary2$sim_mCHtoGFP<-sim_mCHtoGFP
genesummary2$sign_compa<-sign_compa


svg("prot-mRNA-QTLeffect.svg",width = 6,height = 6)
finalcumul45_142<-finalcumul45_14[order(finalcumul45_14$cexallF),]
finalcumul45_142$colF2<-finalcumul45_142$colF
finalcumul45_142$colF2[finalcumul45_142$colF=="#FF0000"]<-"#FF000066"
finalcumul45_142$colF2[finalcumul45_142$colF=="#000000"]<-"#00000066"
finalcumul45_142$colF2[finalcumul45_142$colF=="#00CC80"]<-"#00CC8066"
finalcumul45_142$colF2[finalcumul45_142$colF=="#0000FF"]<-"#0000FF66"
plot(finalcumul45_142$GFP_effect,finalcumul45_142$mCH_effect,
     xlim=c(-1,1),ylim=c(-1,1),col=finalcumul45_142$colF2,cex=finalcumul45_142$cexallF,pch=16,
     xlab="AF(BY/all) difference (high-low) - GFP gate", ylab="AF(BY/all) difference (high-low) - mCherry gate",
     main="QTL effect: protein vs mRNA")
points(finalcumul45_142$GFP_effect,finalcumul45_142$mCH_effect,
       col=finalcumul45_142$colF,cex=finalcumul45_142$cexallF,pch=1)
abline(v=0)
abline(h=0)
abline(a=0,b=1)
reg<-lm(finalcumul45_14$mCH_effect~finalcumul45_14$GFP_effect)
abline(reg,col="grey",lty=2)
summary(reg)
cor.test(finalcumul45_14$GFP_effect,finalcumul45_14$mCH_effect,method = "p")
cor.test(finalcumul45_14$GFP_effect,finalcumul45_14$mCH_effect,method = "s")
#sum(sign(finalcumul45_14$GFP_effect)==sign(finalcumul45_14$mCH_effect))/nrow(finalcumul45_14)
binom.test(sum(sign(finalcumul45_14$GFP_effect)==sign(finalcumul45_14$mCH_effect)), nrow(finalcumul45_14), alternative="greater")
points(finalcumul45_14false$GFP_effect,finalcumul45_14false$mCH_effect,
       xlim=c(-1,1),ylim=c(-1,1),col="grey",cex=finalcumul45_14false$cexallF,pch=1)
dev.off()

summary(as.factor(finalcumul45_142$type))

finalcumul45_test<-finalcumul45_14[finalcumul45_14$type=="gfp",]
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "p")
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "s")
#sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect))/nrow(finalcumul45_test)
binom.test(sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect)), nrow(finalcumul45_test), alternative="greater")


finalcumul45_test<-finalcumul45_14[finalcumul45_14$type=="mch",]
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "p")
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "s")
#sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect))/nrow(finalcumul45_test)
binom.test(sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect)), nrow(finalcumul45_test), alternative="greater")


finalcumul45_test<-finalcumul45_14[finalcumul45_14$type%in%c("gfp","mch"),]
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "p")
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "s")
#sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect))/nrow(finalcumul45_test)
binom.test(sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect)), nrow(finalcumul45_test), alternative="greater")

finalcumul45_test<-finalcumul45_14[finalcumul45_14$type%in%c("together"),]
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "p")
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "s")
#sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect))/nrow(finalcumul45_test)
binom.test(sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect)), nrow(finalcumul45_test), alternative="greater")

finalcumul45_test<-finalcumul45_14[finalcumul45_14$type%in%c("gfp","mch","oposite"),]
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "p")
cor.test(finalcumul45_test$GFP_effect,finalcumul45_test$mCH_effect,method = "s")
#sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect))/nrow(finalcumul45_test)
binom.test(sum(sign(finalcumul45_test$GFP_effect)==sign(finalcumul45_test$mCH_effect)), nrow(finalcumul45_test), alternative="greater")

