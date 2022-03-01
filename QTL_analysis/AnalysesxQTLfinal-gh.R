# Author: Christian Brion - 2020 - UMN
#
#-merge QTL data across replicates
#-save merged QTL data in txt file
#-save merged data in R files
#-provide figures for tagged gene
#-compare QTL curve to previous QTL analysis and save output summary text file
#

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
xpQTL2<-merge(xpQTL,geneInfo,by.x=1,by.y=1,all.x =T,all.y=F)
eQTL_OI<-final_eQTL[final_eQTL$chr.x == chromosome[10] & final_eQTL$maxQTLpos > 130000  & final_eQTL$maxQTLpos < 170000,] #TIF2
pQTL_OI<-xpQTL2[gsub(pattern = "chr", "", xpQTL2$chromosome) == "10" & xpQTL2$peakPosition > 130000  & xpQTL$peakPosition < 170000,] #TIF2
QTG_OI<-geneQTG[geneQTG$chr == 10 & geneQTG$posup > 130000  & geneQTG$posup < 170000,] #TIF2
final_eQTL$chr.x2<-as.numeric(factor(final_eQTL$chr.x,levels = chromosome))
final_eQTLACT1<-final_eQTL[final_eQTL$name=="ACT1",]
final_eQTLtrans<-final_eQTL[final_eQTL$cis==F,]
xpQTL2ACT1<-xpQTL2[xpQTL2$geneName=="ACT1",]


#table with summary of all the QTL experiment
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


#first table to provide the complete in the loop
initialtable <- read.table(paste(alignmentDir,"repQTLresults/YDL022W_GPD1_a_QTLs.txt",sep=""), stringsAsFactors=FALSE, head=TRUE)
initialtable <- initialtable[-(1:nrow(initialtable)),]
compa_GFP<-list()
compa_mCH<-list()

medianC<-c()

#Run for each genes (including NA)
i=2
for (i in 1:length(levels(as.factor(experimentFile2$exp)))) { #loop every 5min
  
  #get info for the gene
  experiment<-levels(as.factor(experimentFile2$exp))[i]
  orf<-experimentFile2$orf[experimentFile2$exp==experiment][1]
  gene<-experimentFile2$gene[experimentFile2$exp==experiment][1]
  print(paste(orf,experiment,date(),paste(i,"/",length(levels(as.factor(experimentFile2$exp))),sep=""),sep=" - "))
  allREPexp<-experimentFile2[experimentFile2$exp==experiment,]
  repinfo<-levels(as.factor(allREPexp$globrep))
  
  #import data
  j=1
  counts<-list()
  multialldata<-list()
  qtlsresult<-list()
  for (j in 1:length(repinfo)) {
    reptemp<-repinfo[j]
    #getdata
    tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment,reptemp, "_1_.*RData$"), full.names=TRUE)
    load(file = tempFile)
    counts[[reptemp]][["G1"]]<-theseCounts
    medianC<-c(medianC,median(theseCounts$Coverage[theseCounts$Coverage!=0],na.rm = T))
    tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment,reptemp, "_2_.*RData$"), full.names=TRUE)
    load(file = tempFile)
    counts[[reptemp]][["G2"]]<-theseCounts
    medianC<-c(medianC,median(theseCounts$Coverage[theseCounts$Coverage!=0],na.rm = T))
    tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment,reptemp, "_3_.*RData$"), full.names=TRUE)
    load(file = tempFile)
    counts[[reptemp]][["G3"]]<-theseCounts
    medianC<-c(medianC,median(theseCounts$Coverage[theseCounts$Coverage!=0],na.rm = T))
    tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment,reptemp, "_4_.*RData$"), full.names=TRUE)
    load(file = tempFile)
    counts[[reptemp]][["G4"]]<-theseCounts
    medianC<-c(medianC,median(theseCounts$Coverage[theseCounts$Coverage!=0],na.rm = T))
    tempFile <- dir(paste(alignmentDir,"AlleleFreqresults/",sep=""), pattern=paste0(experiment,reptemp, "_5_.*RData$"), full.names=TRUE)
    load(file = tempFile)
    counts[[reptemp]][["G5"]]<-theseCounts
    medianC<-c(medianC,median(theseCounts$Coverage[theseCounts$Coverage!=0],na.rm = T))
    
    #load multipooldata
    tempFile <- dir(paste(alignmentDir,"QTLresults/",sep=""), pattern=paste0(experiment,reptemp, "_multi.*RData$"), full.names=TRUE)
    load(file = tempFile)
    multialldata[[reptemp]][["multiPeaksGFP"]]<-multiPeaksGFP
    multialldata[[reptemp]][["multiPeaksmCH"]]<-multiPeaksmCH
    multialldata[[reptemp]][["multipoolOutputGFP"]]<-multipoolOutputGFP
    multialldata[[reptemp]][["multipoolOutputmCH"]]<-multipoolOutputmCH
    
    #load QTL results
    tempFile <- dir(paste(alignmentDir,"QTLresults/",sep=""), pattern=paste0(experiment,reptemp, "_QTLs.txt$"), full.names=TRUE)
    qtltemps<-read.table(tempFile, stringsAsFactors=FALSE, head=TRUE, sep="\t")
    qtlsresult[[reptemp]]<-qtltemps
  }
  
  #fuse results
  exportpeak<-qtlsresult[[1]]
  if (length(repinfo)>1) {
    for (j in 2:length(repinfo)) {
      reptemp<-repinfo[j]
      exportpeak<-rbind(exportpeak,qtlsresult[[j]])
    }
  }
  
  G1<-counts[[1]][["G1"]]
  #get average effect for each QTL
  for (n in 1:nrow(exportpeak)) {
    allLOD<-c()
    alldelta<-c()
    allLODother<-c()
    alldeltaother<-c()
    for (j in 1:length(repinfo)) {
      reptemp<-repinfo[j]
      if (exportpeak$fluorecence[n] =="gfp") {
        HLcurve<-counts[[reptemp]][["G5"]]$roll-counts[[reptemp]][["G4"]]$roll
        HLcurveOther<-counts[[reptemp]][["G3"]]$roll-counts[[reptemp]][["G2"]]$roll
        mp<-multialldata[[reptemp]][["multipoolOutputGFP"]][[exportpeak$chromo[n]]][[2]]
        mpo<-multialldata[[reptemp]][["multipoolOutputmCH"]][[exportpeak$chromo[n]]][[2]]
      } else {
        HLcurveOther<-counts[[reptemp]][["G5"]]$roll-counts[[reptemp]][["G4"]]$roll
        HLcurve<-counts[[reptemp]][["G3"]]$roll-counts[[reptemp]][["G2"]]$roll
        mpo<-multialldata[[reptemp]][["multipoolOutputGFP"]][[exportpeak$chromo[n]]][[2]]
        mp<-multialldata[[reptemp]][["multipoolOutputmCH"]][[exportpeak$chromo[n]]][[2]]
      }
      HLcurve<-HLcurve[G1$chr==chromosome[exportpeak$chromo[n]]]
      HLcurve2<-HLcurve[!(is.na(HLcurve))]
      postemp<-G1$pos[G1$chr==chromosome[exportpeak$chromo[n]]]
      postemp<-postemp[!(is.na(HLcurve))]
      HLcurveOther<-HLcurveOther[G1$chr==chromosome[exportpeak$chromo[n]]]
      HLcurveOther2<-HLcurveOther[!(is.na(HLcurveOther))]
      postempOther<-G1$pos[G1$chr==chromosome[exportpeak$chromo[n]]]
      postempOther<-postempOther[!(is.na(HLcurveOther))]
      bestpos<-abs(postemp-exportpeak$maxIndex[n])==min(abs(postemp-exportpeak$maxIndex[n]))
      bestposOther<-abs(postempOther-exportpeak$maxIndex[n])==min(abs(postempOther-exportpeak$maxIndex[n]))
      allLOD[j]<-mp[mp[,1]==exportpeak$maxIndex[n],2]
      alldelta[j]<-HLcurve2[bestpos][1]
      allLODother[j]<-mpo[mpo[,1]==exportpeak$maxIndex[n],2]
      alldeltaother[j]<-HLcurveOther2[bestposOther][1]
    }
    if (sign(exportpeak$deltaHight_Low[n]) != sign(mean(alldelta))) {
      exportpeak$maxValue[n]<-0
    } else {
      exportpeak$maxValue[n]<-abs(mean(allLOD*sign(alldelta)))
    }
    exportpeak$deltaHight_Low[n]<-mean(alldelta)
    exportpeak$deltaOther[n]<-mean(alldeltaother)
    exportpeak$LOD_Other[n]<-abs(mean(allLODother*sign(alldeltaother)))
  }
  
  #merge peaks in each chanels if they overlap
  exportpeak_merged<-exportpeak[0,]
  exportpeakdestroyed<-exportpeak
  numbAcross<-c()
  range<-150000
  for (n in 1:nrow(exportpeak)) {
    temppeak<-exportpeakdestroyed[exportpeakdestroyed$fluorecence==exportpeak$fluorecence[n] &
                                    sign(exportpeakdestroyed$deltaHight_Low)==sign(exportpeak$deltaHight_Low[n]) &
                                    exportpeakdestroyed$chromo==exportpeak$chromo[n] &
                                    exportpeakdestroyed$maxIndex<exportpeak$maxIndex[n]+range/2 & 
                                    exportpeakdestroyed$maxIndex>exportpeak$maxIndex[n]-range/2,]
    if (nrow(temppeak)>0) {
      temppeak2<-temppeak[1,]
      temppeak2$maxIndex<-mean(temppeak$maxIndex)
      temppeak2$maxValue<-abs(mean(temppeak$maxValue*sign(temppeak$deltaHight_Low)))
      temppeak2$leftIndex<-mean(temppeak$leftIndex)
      temppeak2$rightIndex<-mean(temppeak$rightIndex)
      temppeak2$deltaHight_Low<-mean(temppeak$deltaHight_Low)
      temppeak2$LOD_Other<-abs(mean(temppeak$LOD_Other*sign(temppeak$deltaOther)))
      temppeak2$deltaOther<-mean(temppeak$deltaOther)
    } else {
      temppeak2<-exportpeak[n,]
    }
    temppeak<-exportpeakdestroyed[exportpeakdestroyed$fluorecence==temppeak2$fluorecence[1] &
                                    sign(exportpeakdestroyed$deltaHight_Low)==sign(temppeak2$deltaHight_Low[1]) &
                                    exportpeakdestroyed$chromo==temppeak2$chromo[1] &
                                    exportpeakdestroyed$maxIndex<temppeak2$maxIndex[1]+range/2 & 
                                    exportpeakdestroyed$maxIndex>temppeak2$maxIndex[1]-range/2,]
    exportpeakdestroyed<-exportpeakdestroyed[!(exportpeakdestroyed$fluorecence==temppeak2$fluorecence[1] &
                                                 sign(exportpeakdestroyed$deltaHight_Low)==sign(temppeak2$deltaHight_Low[1]) &
                                                 exportpeakdestroyed$chromo==temppeak2$chromo[1] &
                                                 exportpeakdestroyed$maxIndex<temppeak2$maxIndex[1]+range/2 & 
                                                 exportpeakdestroyed$maxIndex>temppeak2$maxIndex[1]-range/2),]
    if (nrow(temppeak)>0) {
      numbAcross[length(numbAcross)+1]<-nrow(temppeak)
      temppeak2<-temppeak[1,]
      temppeak2$maxIndex<-mean(temppeak$maxIndex)
      temppeak2$maxValue<-abs(mean(temppeak$maxValue*sign(temppeak$deltaHight_Low)))
      temppeak2$leftIndex<-mean(temppeak$leftIndex)
      temppeak2$rightIndex<-mean(temppeak$rightIndex)
      temppeak2$deltaHight_Low<-mean(temppeak$deltaHight_Low)
      temppeak2$LOD_Other<-abs(mean(temppeak$LOD_Other*sign(temppeak$deltaOther)))
      temppeak2$deltaOther<-mean(temppeak$deltaOther)
      exportpeak_merged<-rbind(exportpeak_merged,temppeak2)
    }
  }
  exportpeak_merged$numbAcross<-numbAcross
  exportpeak_merged_signi<-exportpeak_merged[exportpeak_merged$maxValue>3,]
  
  #define the type of RNAvProt. QTL before mergin them
  for (n in 1:nrow(exportpeak_merged_signi)) {
    if (exportpeak_merged_signi$LOD_Other[n]<3) {
      exportpeak_merged_signi$effectRNAPROT[n]<-"specific"
    } else {
      if (sign(exportpeak_merged_signi$deltaHight_Low[n])==sign(exportpeak_merged_signi$deltaOther[n])) {
        exportpeak_merged_signi$effectRNAPROT[n]<-"together"
      } else {
        exportpeak_merged_signi$effectRNAPROT[n]<-"oposite"
      }
    }
  }
  
  
  #merge peaks in accros chanels if they overlap, using green QTL first
  exportpeak_merged_signi_cumul<-exportpeak_merged_signi[0,c(1,2,9,18)]
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
  temp_mCH<-exportpeak_merged_signi[exportpeak_merged_signi$fluorecence=="mch",]
  temp_gfp<-exportpeak_merged_signi[exportpeak_merged_signi$fluorecence=="gfp",]
  for (n in 1:nrow(temp_gfp)) {
    exportpeak_merged_signi_cumul<-rbind(exportpeak_merged_signi_cumul,temp_gfp[n,c(1,2,9,18)])
    
    GFP_LOD[n]<-temp_gfp$maxValue[n]
    GFP_effect[n]<-temp_gfp$deltaHight_Low[n]
    GFP_maxPOS[n]<-temp_gfp$maxIndex[n]
    GFP_leftPOS[n]<-temp_gfp$leftIndex[n]
    GFP_rightPOS[n]<-temp_gfp$rightIndex[n]
    I_mCH<-(1:nrow(temp_mCH))[temp_mCH$orf==temp_gfp$orf[n] & temp_mCH$chromo==temp_gfp$chromo[n] & temp_mCH$rightIndex>min(temp_gfp$leftIndex[n],temp_gfp$maxIndex[n]) & temp_mCH$leftIndex<max(temp_gfp$rightIndex[n],temp_gfp$maxIndex[n])]
    if (length(I_mCH)==1) {
      mCH_LOD[n]<-temp_mCH$maxValue[I_mCH]
      mCH_effect[n]<-temp_mCH$deltaHight_Low[I_mCH]
      mCH_maxPOS[n]<-temp_mCH$maxIndex[I_mCH]
      mCH_leftPOS[n]<-temp_mCH$leftIndex[I_mCH]
      mCH_rightPOS[n]<-temp_mCH$rightIndex[I_mCH]
      if (sign(mCH_effect[n])==sign(GFP_effect[n])) {
        type[n]<-"together"
      } else {
        type[n]<-"oposite"
      }
      temp_mCH<-temp_mCH[-I_mCH,]
    } else if (length(I_mCH)>1) {
      print(temp_mCH[I_mCH,])
      break()
    } else {
      mCH_LOD[n]<-temp_gfp$LOD_Other[n]
      mCH_effect[n]<-temp_gfp$deltaOther[n]
      mCH_maxPOS[n]<-temp_gfp$maxIndex[n]
      mCH_leftPOS[n]<-temp_gfp$leftIndex[n]
      mCH_rightPOS[n]<-temp_gfp$rightIndex[n]
      type[n]<-paste("gfp",sep="")
    }
  }
  tempdata<-data.frame(type,GFP_LOD,GFP_effect,GFP_maxPOS,GFP_leftPOS,GFP_rightPOS
                       ,mCH_LOD,mCH_effect,mCH_maxPOS,mCH_leftPOS,mCH_rightPOS)
  
  #merge peaks in accros chanels if they overlap, using red QTL second
  exportpeak_merged_signi_cumul2<-cbind(exportpeak_merged_signi_cumul,tempdata)
  n<-nrow(tempdata)
  for (j in 1:nrow(temp_mCH)) {
    exportpeak_merged_signi_cumul<-rbind(exportpeak_merged_signi_cumul,temp_mCH[j,c(1,2,9,18)])
    mCH_LOD[n+j]<-temp_mCH$maxValue[j]
    mCH_effect[n+j]<-temp_mCH$deltaHight_Low[j]
    mCH_maxPOS[n+j]<-temp_mCH$maxIndex[j]
    mCH_leftPOS[n+j]<-temp_mCH$leftIndex[j]
    mCH_rightPOS[n+j]<-temp_mCH$rightIndex[j]
    GFP_LOD[n+j]<-temp_mCH$LOD_Other[j]
    GFP_effect[n+j]<-temp_mCH$deltaOther[j]
    GFP_maxPOS[n+j]<-temp_mCH$maxIndex[i]
    GFP_leftPOS[n+j]<-temp_mCH$leftIndex[j]
    GFP_rightPOS[n+j]<-temp_mCH$rightIndex[j]
    if (temp_mCH$effectRNAPROT[j]!="specific") {
      type[n+j]<-paste(temp_mCH$effectRNAPROT[j],"_mch",sep="") #"mch"#
    } else {
      type[n+j]<-"mch"
    }
  }
  
  tempdata<-data.frame(type,GFP_LOD,GFP_effect,GFP_maxPOS,GFP_leftPOS,GFP_rightPOS
                       ,mCH_LOD,mCH_effect,mCH_maxPOS,mCH_leftPOS,mCH_rightPOS)
  
  exportpeak_merged_signi_cumul2<-cbind(exportpeak_merged_signi_cumul,tempdata)
  exportpeak_merged_signi_cumul3<-exportpeak_merged_signi_cumul2[order(exportpeak_merged_signi_cumul2$orf,exportpeak_merged_signi_cumul2$chromo,exportpeak_merged_signi_cumul2$GFP_maxPOS),]
  
  
  write.table(exportpeak_merged_signi,paste(alignmentDir,"repQTLresults/",paste(orf,experiment, sep="_"), "_QTLs.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)
  write.table(exportpeak_merged_signi_cumul3,paste(alignmentDir,"repQTLresults/",paste(orf,experiment, sep="_"), "_QTLs_cumul.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)
  
  
  initialtable<-rbind(initialtable,exportpeak_merged_signi)
  
  #compare data to known frank QTL
  
  # chro<-c()
  # pos<-c()
  # meandAFg<-c()
  # meanLODg<-c()
  # meandAFm<-c()
  # meanLODm<-c()
  # n=1
  # for (c in 1:length(chromosome)) {
  #   for (p in 1:length(counts[[repinfo[1]]][["G5"]]$roll[counts[[repinfo[1]]][["G5"]]$chr==chromosome[c]])) {
  #     chro[n]<-c
  #     pos[n]<-counts[[repinfo[1]]][["G5"]]$pos[counts[[repinfo[1]]][["G5"]]$chr==chromosome[c]][p]
  #     G5<-c()
  #     G4<-c()
  #     LODg<-c()
  #     G3<-c()
  #     G2<-c()
  #     LODm<-c()
  #     for (j in 1:length(repinfo)) {
  #       G5[j]<-counts[[repinfo[j]]][["G5"]]$roll[counts[[repinfo[1]]][["G5"]]$chr==chromosome[c] & counts[[repinfo[1]]][["G5"]]$pos==counts[[repinfo[1]]][["G5"]]$pos[counts[[repinfo[1]]][["G5"]]$chr==chromosome[c]][p]]
  #       G4[j]<-counts[[repinfo[j]]][["G4"]]$roll[counts[[repinfo[1]]][["G4"]]$chr==chromosome[c] & counts[[repinfo[1]]][["G4"]]$pos==counts[[repinfo[1]]][["G5"]]$pos[counts[[repinfo[1]]][["G5"]]$chr==chromosome[c]][p]]
  #       G3[j]<-counts[[repinfo[j]]][["G3"]]$roll[counts[[repinfo[1]]][["G3"]]$chr==chromosome[c] & counts[[repinfo[1]]][["G3"]]$pos==counts[[repinfo[1]]][["G5"]]$pos[counts[[repinfo[1]]][["G5"]]$chr==chromosome[c]][p]]
  #       G2[j]<-counts[[repinfo[j]]][["G2"]]$roll[counts[[repinfo[1]]][["G2"]]$chr==chromosome[c] & counts[[repinfo[1]]][["G2"]]$pos==counts[[repinfo[1]]][["G5"]]$pos[counts[[repinfo[1]]][["G5"]]$chr==chromosome[c]][p]]
  #       val<-ceiling(pos[n]/100)
  #       LODg[j]<-multialldata[[repinfo[j]]][["multipoolOutputGFP"]][[c]][[2]][val,2]
  #       LODm[j]<-multialldata[[repinfo[j]]][["multipoolOutputmCH"]][[c]][[2]][val,2]
  #     }
  #     meandAFg[n]<-mean(G5, na.rm=T)-mean(G4,na.rm=T)
  #     meanLODg[n]<-mean(LODg, na.rm=T)
  #     meandAFm[n]<-mean(G3, na.rm=T)-mean(G2,na.rm=T)
  #     meanLODm[n]<-mean(LODm, na.rm=T)
  #     n=n+1
  #   }
  #   print(c)
  # } #15seg
  # mergeQTLplot<-data.frame(chro,pos,meandAFg,meanLODg,meandAFm,meanLODm)
  # save(mergeQTLplot,file = paste("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/190208-AllQTLprocess/repQTLresults/","mergeQTLplot_",orf,"_",gene,".RData",sep=""))
  # 
  load(paste("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/190208-AllQTLprocess/repQTLresults/","mergeQTLplot_",orf,"_",gene,".RData",sep="")) #mergeQTLplot
  
  #GFP comparizon
  chr<-c()
  brionLOD<-c()
  briondeltaAF<-c()
  brionmaxpos<-c()
  albertLOD<-c()
  albertdeltaAF<-c()
  albertmaxpos<-c()
  correspond<-c()
  repnum<-c()
  n=0
  load(paste("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/190208-AllQTLprocess/frank_QTL_data/pQTLfrank_",orf,"_",gene,".RData",sep="")) #pQTLfrank
  tempbrion<-exportpeak_merged_signi[exportpeak_merged_signi$fluorecence == "gfp",]
  tempalbert<-xpQTL2[xpQTL2$gene==orf,]
  tempalbert$chromosome<-as.numeric(gsub(pattern = "chr", "", tempalbert$chromosome))
  for (j in 1:nrow(tempbrion)) {
    chr[j]<-tempbrion$chromo[j]
    brionLOD[j]<-tempbrion$maxValue[j]
    briondeltaAF[j]<-tempbrion$deltaHight_Low[j]
    brionmaxpos[j]<-tempbrion$maxIndex[j]
    repnum[j]<-tempbrion$numbAcross[j]
    I_albert<-(1:nrow(tempalbert))[tempalbert$chromosome==tempbrion$chromo[j]
                                   & ((tempalbert$X2LODIntervalRight>min(tempbrion$leftIndex[j],tempbrion$maxIndex[j]) & tempalbert$X2LODIntervalLeft<max(tempbrion$rightIndex[j],tempbrion$maxIndex[j]))
                                      | abs(tempalbert$peakPosition-tempbrion$maxIndex[j])<20000 )]
    if (length(I_albert)==1) {
      albertLOD[j]<-tempalbert$LOD[I_albert]
      albertdeltaAF[j]<-tempalbert$alleleFrequencyDifference[I_albert]
      albertmaxpos[j]<-tempalbert$peakPosition[I_albert]
      if (sign(albertdeltaAF[j])==sign(tempbrion$deltaHight_Low[j])) {
        correspond[j]<-"together"
      } else {
        correspond[j]<-"oposite"
      }
      tempalbert<-tempalbert[-I_albert,]
    } else if (length(I_albert)>1) {
      print(tempalbert[I_albert,])
      break()
    } else {
      pQTLfranktemp<-pQTLfrank[pQTLfrank$chr==chr[j],]
      pQTLfranktemp<-pQTLfranktemp[order(abs(pQTLfranktemp$pos-brionmaxpos[j])),]
      albertLOD[j]<-0
      albertdeltaAF[j]<-mean(pQTLfranktemp$dAF[1:5],na.rm=T)
      albertmaxpos[j]<-pQTLfranktemp$pos[1]
      correspond[j]<-paste("unique_brion")
    }
  }
  n=j
  if(nrow(tempalbert)>0) {
    for (j in 1:nrow(tempalbert)) {
      chr[n+j]<-tempalbert$chromosome[j]
      albertLOD[n+j]<-tempalbert$LOD[j]
      albertdeltaAF[n+j]<-tempalbert$alleleFrequencyDifference[j]
      albertmaxpos[n+j]<-tempalbert$peakPosition[j]
      pQTLmytemp<-mergeQTLplot[mergeQTLplot$chro==chr[n+j],]
      pQTLmytemp<-pQTLmytemp[order(abs(pQTLmytemp$pos-albertmaxpos[n+j])),]
      brionLOD[n+j]<-pQTLmytemp$meanLODg[1]
      briondeltaAF[n+j]<-mean(pQTLmytemp$meandAFg[1:5],na.rm=T)
      brionmaxpos[n+j]<-pQTLmytemp$pos[1]
      repnum[n+j]<-0
      correspond[n+j]<-paste("unique_albert")
    }
  }
  genes<-rep(orf,length(correspond))
  names<-rep(gene,length(correspond))
  tempcompaGFP<-data.frame(genes,names,chr,brionmaxpos,briondeltaAF,brionLOD,repnum,albertmaxpos,albertdeltaAF,albertLOD,correspond)
  compa_GFP[[orf]]<-tempcompaGFP
  sum(sign(tempcompaGFP$briondeltaAF)==sign(tempcompaGFP$albertdeltaAF))/nrow(tempcompaGFP)
  
  
  #mCH comparizon
  chr<-c()
  brionLOD<-c()
  briondeltaAF<-c()
  brionmaxpos<-c()
  albertLOD<-c()
  albertR<-c()
  albertmaxpos<-c()
  correspond<-c()
  repnum<-c()
  n=0
  load(paste("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/190208-AllQTLprocess/frank_QTL_data/eQTLfrank_",orf,"_",gene,".RData",sep="")) #eQTLfrank
  tempbrion<-exportpeak_merged_signi[exportpeak_merged_signi$fluorecence == "mch",]
  tempalbert<-final_eQTLtrans[final_eQTLtrans$gene==orf,]
  for (j in 1:nrow(tempbrion)) {
    chr[j]<-tempbrion$chromo[j]
    brionLOD[j]<-tempbrion$maxValue[j]
    briondeltaAF[j]<-tempbrion$deltaHight_Low[j]
    brionmaxpos[j]<-tempbrion$maxIndex[j]
    repnum[j]<-tempbrion$numbAcross[j]
    I_albert<-(1:nrow(tempalbert))[tempalbert$chr.x2==tempbrion$chromo[j]
                                   & ((tempalbert$maxQTLright>min(tempbrion$leftIndex[j],tempbrion$maxIndex[j]) & tempalbert$maxQTLleft<max(tempbrion$rightIndex[j],tempbrion$maxIndex[j]))
                                      | abs(tempalbert$maxQTLpos-tempbrion$maxIndex[j])<20000 )]
    if (length(I_albert)>0) {
      albertLOD[j]<-tempalbert$LOD[I_albert[1]]
      albertR[j]<- -tempalbert$r[I_albert[1]]
      albertmaxpos[j]<-tempalbert$maxQTLpos[I_albert[1]]
      if (sign(albertR[j])==sign(tempbrion$deltaHight_Low[j])) {
        correspond[j]<-"together"
      } else {
        correspond[j]<-"oposite"
      }
      tempalbert<-tempalbert[-I_albert[1],]
    # } else if (length(I_albert)>1) {
    #   print(tempalbert[I_albert,])
    #   break()
    } else {
      eQTLfranktemp<-eQTLfrank[eQTLfrank$chr==chr[j],]
      eQTLfranktemp<-eQTLfranktemp[order(abs(eQTLfranktemp$pos-brionmaxpos[j])),]
      albertLOD[j]<-0
      albertR[j]<- -mean(eQTLfranktemp$r_val[1:5],na.rm=T)
      albertmaxpos[j]<-eQTLfranktemp$pos[1]
      correspond[j]<-paste("unique_brion")
    }
  }
  n=j
  if(nrow(tempalbert)>0) {
    for (j in 1:nrow(tempalbert)) {
      chr[n+j]<-tempalbert$chr.x2[j]
      albertLOD[n+j]<-tempalbert$LOD[j]
      albertR[n+j] <- -tempalbert$r[j]
      albertmaxpos[n+j] <- tempalbert$maxQTLpos[j]
      eQTLmytemp<-mergeQTLplot[mergeQTLplot$chro==chr[n+j],]
      eQTLmytemp<-eQTLmytemp[order(abs(eQTLmytemp$pos-albertmaxpos[n+j])),]
      brionLOD[n+j] <- eQTLmytemp$meanLODm[1]
      briondeltaAF[n+j] <- mean(eQTLmytemp$meandAFm[1:5],na.rm=T)
      brionmaxpos[n+j]<-eQTLmytemp$pos[1]
      repnum[n+j]<-0
      correspond[n+j]<-paste("unique_albert")
    }
  }
  genes<-rep(orf,length(correspond))
  names<-rep(gene,length(correspond))
  tempcompamCH<-data.frame(genes,names,chr,brionmaxpos,briondeltaAF,brionLOD,repnum,albertmaxpos,albertR,albertLOD,correspond)
  compa_mCH[[orf]]<-tempcompamCH
  sum(sign(tempcompamCH$briondeltaAF)==sign(tempcompamCH$albertR))/nrow(tempcompamCH)
  
  #}
  
  #now ploting in PDF
  pdf(file = paste(alignmentDir,"repQTLresults/",paste(orf,experiment, sep="_"), "_plot.pdf", sep=""), width=11, height=16)
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
  colreds<-c(rgb(1,0,0),rgb(0.9,0,0),rgb(0.8,0,0),rgb(0.7,0,0),rgb(1,0.1,0.1),rgb(1,0.2,0.2),rgb(1,0.3,0.3),rgb(1,0.4,0.4))
  colgreens<-c(rgb(0,0.8,0.5),rgb(0,0.75,0.45),rgb(0,0.7,0.4),rgb(0,0.64,0.35),rgb(0,0.85,0.55),rgb(0,0.9,0.6),rgb(0,1,0.65),rgb(0.5,1,0.7))
  
  plot(gcoords, rep(0, length(gcoords)), ylim=c(-1,1), main = paste(orf,experiment, sep="_"), xaxt='n', xlab="chromosome", ylab="AF(BY/all) difference (high-low)", type="n")
  #abline(v = getGcoords("chrV", 33466, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # CAN1 is on chr05, pos 33466
  #abline(v = getGcoords("chrIII", 198671, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # MATALPHA is on chr03, 198671
  abline(v = getGcoords(geneInfo[orf, "chr"], mean(as.numeric(geneInfo[orf, c("start", "end")])), sepBetweenChr, sgd.table=sgd_table), lwd=2, col="purple")
  for (repn in 1:length(repinfo)) {
    reptemp<-repinfo[repn]
    Gdata<-counts[[reptemp]]
    for (j in unique(G1[,"chr"])){
      points(gcoords[G1$chr == j], Gdata[["G3"]]$roll[G1$chr == j]-Gdata[["G2"]]$roll[G1$chr == j], type="l", lwd=2,col=colreds[repn],lty=1)
      points(gcoords[G1$chr == j], Gdata[["G5"]]$roll[G1$chr == j]-Gdata[["G4"]]$roll[G1$chr == j], type="l", lwd=2,col=colgreens[repn],lty=1)
    }
  }
  for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
  }
  abline(h = 0, lty=2, col="light blue")
  # these difference thresholds from the nullSorts:
  abline(h = AFThres, lty = 2, lwd=2, col = c(1,1,1,0.8))
  abline(h = -AFThres, lty = 2, lwd=2, col = c(1,1,1,0.8))
  axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)
  
  ylimMaxGFP<-c()
  ylimMaxmch<-c()
  for (repn in 1:length(repinfo)) {
    reptemp<-repinfo[repn]
    ylimMaxGFP[repn] <- max(c(multiThres, sapply(multialldata[[reptemp]][["multipoolOutputGFP"]], function(x){max(x[[2]][,2])}))) + 1
    ylimMaxmch[repn]<- max(c(multiThres, sapply(multialldata[[reptemp]][["multipoolOutputmCH"]], function(x){max(x[[2]][,2])}))) + 1
  }
  ylimMax <- max(ylimMaxGFP,ylimMaxmch)
  ylimMax <- min(ylimMax,70)
  plot(gcoords, rep(0, length(gcoords)), main = paste(orf,experiment, sep="_"), type="n", xaxt='n', xlab="chromosome", ylab="Multipool LOD", ylim=c(0, ylimMax))
  abline(v = getGcoords(geneInfo[orf, "chr"], mean(as.numeric(geneInfo[orf, c("start", "end")])), sepBetweenChr, sgd.table=sgd_table), lwd=2, col="purple")
  for (repn in 1:length(repinfo)) {
    reptemp<-repinfo[repn]
    Mdata<-multialldata[[reptemp]]
    for (j in 1:16){
      points(getGcoords(paste0("chr", as.roman(j)), Mdata[["multipoolOutputGFP"]][[j]][[2]][,1], sepBetweenChr, sgd.table=sgd_table), Mdata[["multipoolOutputGFP"]][[j]][[2]][,2], type="l", lwd=2, col=colgreens[repn])
      points(getGcoords(paste0("chr", as.roman(j)), Mdata[["multipoolOutputmCH"]][[j]][[2]][,1], sepBetweenChr, sgd.table=sgd_table), Mdata[["multipoolOutputmCH"]][[j]][[2]][,2], type="l", lwd=2, col=colreds[repn])
    }
  }
  # add stars at avearge peaks
  for (n in 1:nrow(exportpeak_merged_signi)) {
    thisPeakPlotPos <- getGcoords(paste0("chr", as.roman(exportpeak_merged_signi$chromo[n])), exportpeak_merged_signi$maxIndex[n], sepBetweenChr, sgd.table=sgd_table)
    if (exportpeak_merged_signi$fluorecence[n]=="gfp") {colstar<-"*"} else {colstar<-"+"}
    text(thisPeakPlotPos, exportpeak_merged_signi$maxValue[n], labels=colstar, col="black", cex=3)
  }
  
  abline(h = 0, lty=2, col="light blue")
  for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
  }
  axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)
  abline(h = multiThres, lty=2, col="black", lwd=2)
  dev.off()
  
  #plot curve on JPG
  typeQTLs<-c("gfp","mch","together","oposite","oposite_mch","together_mch")
  colQTLs<-c(rgb(0,0.8,0.5),rgb(1,0,0),rgb(0,0,0),rgb(0,0,1),rgb(0,0.7,1),rgb(1,0.6,0))
  jpeg(filename = paste(alignmentDir,"repQTLresults/",paste(orf,experiment, sep="_"), "_plot.jpeg", sep=""),width = 1800,height = 1000,quality = 99)
  plot(gcoords, rep(0, length(gcoords)), ylim=c(-1,1), main = paste(orf,experiment, sep="_"), xaxt='n', xlab="chromosome", ylab="AF(BY/all) difference (high-low)", type="n",cex.main=2,cex.axis=2,cex.lab=2)
  #abline(v = getGcoords("chrV", 33466, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # CAN1 is on chr05, pos 33466
  #abline(v = getGcoords("chrIII", 198671, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # MATALPHA is on chr03, 198671
  abline(v = getGcoords(geneInfo[orf, "chr"], mean(as.numeric(geneInfo[orf, c("start", "end")])), sepBetweenChr, sgd.table=sgd_table), lwd=3, col="purple")
  for (repn in 1:length(repinfo)) {
    reptemp<-repinfo[repn]
    Gdata<-counts[[reptemp]]
    for (j in unique(G1[,"chr"])){
      points(gcoords[G1$chr == j], Gdata[["G3"]]$roll[G1$chr == j]-Gdata[["G2"]]$roll[G1$chr == j], type="l", lwd=3,col=colreds[repn],lty=1)
      points(gcoords[G1$chr == j], Gdata[["G5"]]$roll[G1$chr == j]-Gdata[["G4"]]$roll[G1$chr == j], type="l", lwd=3,col=colgreens[repn],lty=1)
    }
  }
  for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
  }
  abline(h = 0, lty=2, col="light blue")
  # these difference thresholds from the nullSorts:
  abline(h = AFThres, lty = 2, lwd=3, col = c(1,1,1,0.8))
  abline(h = -AFThres, lty = 2, lwd=3, col = c(1,1,1,0.8))
  exportpeak_merged_signi_cumul4<-exportpeak_merged_signi_cumul3[exportpeak_merged_signi_cumul3$GFP_LOD>4.5 | exportpeak_merged_signi_cumul3$mCH_LOD > 4.5,]
  for (j in 1:nrow(exportpeak_merged_signi_cumul4)){
    if (exportpeak_merged_signi_cumul4$numbAcross[j]==1 & !(exportpeak_merged_signi_cumul4$gene[j] %in% c("RPS10A","TDH3","NA"))) {pchtem<-12} else {pchtem<-15}
    colQTL<-colQTLs[typeQTLs==exportpeak_merged_signi_cumul4$type[j]]
    ypos<-(max(abs(exportpeak_merged_signi_cumul4$GFP_effect[j]),abs(exportpeak_merged_signi_cumul4$mCH_effect[j]))+0.1)*sign(exportpeak_merged_signi_cumul4$GFP_effect[j]+exportpeak_merged_signi_cumul4$mCH_effect[j])
    points(getGcoords(chromosome[exportpeak_merged_signi_cumul4$chromo[j]],mean(exportpeak_merged_signi_cumul4$GFP_maxPOS[j],exportpeak_merged_signi_cumul4$mCH_maxPOS[j]), sepBetweenChr, sgd.table=sgd_table),ypos,col=colQTL,pch=pchtem,cex=3)
  }
  axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F,cex.axis=2)
  dev.off()
  
  #plot curve on svg
  svg(file = paste(alignmentDir,"repQTLresults/",paste(orf,experiment, sep="_"), "_plot.svg", sep=""), width=8, height=4.6)
  plot(gcoords, rep(0, length(gcoords)), ylim=c(-1,1), main = paste(orf,experiment, sep="_"), xaxt='n', xlab="chromosome", ylab="AF(BY/all) difference (high-low)", type="n")
  #abline(v = getGcoords("chrV", 33466, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # CAN1 is on chr05, pos 33466
  #abline(v = getGcoords("chrIII", 198671, sepBetweenChr, sgd.table=sgd_table), lwd=2, col="darkgrey") # MATALPHA is on chr03, 198671
  abline(v = getGcoords(geneInfo[orf, "chr"], mean(as.numeric(geneInfo[orf, c("start", "end")])), sepBetweenChr, sgd.table=sgd_table), lwd=3, col="purple")
  for (repn in 1:length(repinfo)) {
    reptemp<-repinfo[repn]
    Gdata<-counts[[reptemp]]
    for (j in unique(G1[,"chr"])){
      #points(gcoords[G1$chr == j], Gdata[["G3"]]$roll[G1$chr == j]-Gdata[["G2"]]$roll[G1$chr == j], type="l", lwd=3,col=colreds[repn],lty=1)
      points(gcoords[G1$chr == j][(1:1000)*10], Gdata[["G3"]]$roll[G1$chr == j][(1:1000)*10]-Gdata[["G2"]]$roll[G1$chr == j][(1:1000)*10], type="l", lwd=2.1,col=colreds[repn],lty=1)
      #points(gcoords[G1$chr == j], Gdata[["G5"]]$roll[G1$chr == j]-Gdata[["G4"]]$roll[G1$chr == j], type="l", lwd=3,col=colgreens[repn],lty=1)
      points(gcoords[G1$chr == j][(1:1000)*10], Gdata[["G5"]]$roll[G1$chr == j][(1:1000)*10]-Gdata[["G4"]]$roll[G1$chr == j][(1:1000)*10], type="l", lwd=2.1,col=colgreens[repn],lty=1)
    }
  }
  for (j in chrCutoffs){
    abline(v = j, lty=2, col="light blue")
  }
  abline(v = {gcoords[G1$chr == "chrXVI"][length(gcoords[G1$chr == "chrXVI"])] + sepBetweenChr/2}, lty=2, col="light blue")
  abline(h = 0, lty=2, col="light blue")
  # these difference thresholds from the nullSorts:
  abline(h = AFThres, lty = 2, lwd=1.5, col = c(1,1,1,0.8))
  abline(h = -AFThres, lty = 2, lwd=1.5, col = c(1,1,1,0.8))
  axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)
  dev.off()
  
}

medianC2<-medianC[c(1:(5*19),(5*19+6):130)]

write.table(initialtable,paste(alignmentDir,"All_QTL_fused.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)


compa_GFPdf<-compa_GFP[[names(compa_GFP)[1]]][0,]
compa_mCHdf<-compa_mCH[[names(compa_mCH)[1]]][0,]
for (j in names(compa_GFP)) {
  compa_GFPdf<-rbind(compa_GFPdf,compa_GFP[[j]])
  compa_mCHdf<-rbind(compa_mCHdf,compa_mCH[[j]])
}
compa_GFPdf$channel<-rep("gfp",nrow(compa_GFPdf))
compa_mCHdf$channel<-rep("mch",nrow(compa_mCHdf))

write.table(compa_GFPdf,paste(alignmentDir,"compa_GFPdf.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)
write.table(compa_mCHdf,paste(alignmentDir,"compa_mCHdf.txt", sep=""),quote = F,sep = "\t",row.names = F,col.names = T)
