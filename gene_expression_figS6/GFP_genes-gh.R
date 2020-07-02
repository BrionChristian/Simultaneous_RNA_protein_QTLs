#import, merge and plot GFP tag and RNA sequencing expression data from Albert, Bloom et al. 2018 and Huh et al 2003

setwd("YOUR/WORKING/DIRECTORY")

#importation of the DATA (gene expression level)
GFPdata<-read.table(paste("GFPlevels.txt",sep=""),header = T,sep = "\t",na.strings = c("NaN","NA"),quote = "\"",stringsAsFactors = F) #GFP level from Huh et al 2003
UTRdata<-read.table(paste("UTRfile.txt",sep=""),header = T,sep = "\t",na.strings = c("Not available ","NA"),quote = "\"",stringsAsFactors = F) #Size of UTR based on NRA sequencing
GFPdata$Name[GFPdata$ORF=="YKL134C"]<-"OCT1"
keep<-colnames(GFPdata)
media<-"SD" #"YEPD""SD" #which media to use the GFP data from
if (media=="SD") {keep2<-keep[c(1,2,5,6,16)]} else {keep2<-keep[c(1,2,3,4,15)]}
GFPdataSD<-GFPdata[,colnames(GFPdata)%in%keep2]
colnames(GFPdataSD)<-c("ORF","Name","GFP","error","Reason_Absent")

genesINFO<-read.table(paste("geneINFOsc.txt", sep=""), sep="\t", head=TRUE, stringsAsFactors=FALSE, quote = "") #get gene annotation
genesINFO2<-genesINFO[genesINFO$type=="ORF",]
expDATA<-read.table(paste("difExp.txt", sep=""), sep="\t", head=TRUE, stringsAsFactors=FALSE, quote = "") #get expression data from my own RNAseq analysis
expDATA2<-expDATA[,c(1,2)]
expDATAFrank<-read.table(paste("meanEXP.txt", sep=""), sep="\t", head=TRUE, stringsAsFactors=FALSE, quote = "") #get expression data from Albert, Bloom, et al 2018
GOI <- read.table("GOI.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE)
GOI <- rbind(GOI,c("TDH3","YGR192C","GFP"))


#merge the data in one big data.frame
DATA<-merge(genesINFO2,GFPdataSD,by.x=1,by.y=1,all.x=T,all.y=F)
notther<-GFPdataSD[!(GFPdataSD$ORF%in%genesINFO2$sysname),]
for (i in 1:nrow(notther)) {
  if (length(grep(notther$ORF[i],DATA$alias))>0) {
    n<-(1:nrow(DATA))[grepl(notther$ORF[i],DATA$alias)]
    if (length(n)!=1) {print(n)}
    DATA[n,ncol(DATA)-(3:0)]<-notther[i,2:5]
  }
}
DATA2<-merge(DATA,expDATA2,by.x=1,by.y=1,all.x=T,all.y=F)
notther2<-expDATA2[!(expDATA2$Gene%in%DATA$sysname),]
DATA3<-merge(DATA2,expDATAFrank,by.x=1,by.y=1,all.x=T,all.y=F)
notther2<-expDATAFrank[!(expDATAFrank$ORF%in%DATA2$sysname),]
noquant<-DATA3[is.na(DATA3$Reason_Absent) & is.na(DATA3$logmeanRNA),]

#prepare data for ploting
base<-log2(20)
dataquant<-DATA3[!(is.na(DATA3$Reason_Absent) & is.na(DATA3$logmeanRNA)),]
dataquant$logGFP<-log2(dataquant$GFP)
#dataquant$logRNA<-log2(dataquant$meanexp)
dataquant$logRNA<-dataquant$logmeanRNA
dataquant$logGFP[is.na(dataquant$logGFP)]<-0
dataquant$logGFP[grepl("_Low",dataquant$Reason_Absent)]<-base
dataquant<-dataquant[order(dataquant$logRNA),]
dataquant<-dataquant[order(dataquant$logGFP),]
n_noquant<-sum(dataquant$logGFP==0)
n_toolow<-sum(dataquant$logGFP==base)
dataquant$logGFP[dataquant$logGFP==0]<-base
dataquant$logGFP<-dataquant$logGFP-base
dataquant[dataquant$name.x%in%c("ACT1","GPD1"),]
ACT1RNAthr<-dataquant$logRNA[dataquant$name.x%in%c("ACT1")]-1
sum(dataquant$logRNA>ACT1RNAthr,na.rm = T)/nrow(dataquant)
RNAfact<-1.3

#ploting
pdf(file = paste("GFP-RNAthreshold_",media,".pdf",sep=""),width = 11,height = 6)
plot(0,0,typ="n",xlim=c(0,nrow(dataquant)),ylim=c(0,12),xaxt="n",yaxt="n",main=media)
x<-c((1:nrow(dataquant))[dataquant$logGFP>0],max((1:nrow(dataquant))[dataquant$logGFP>0]),min((1:nrow(dataquant))[dataquant$logGFP>0]))
y<-c(dataquant$logGFP[dataquant$logGFP>0],0,0)
obj<-data.frame(x,y)
polygon(obj,col="#34ab83ff",border = NA)
points(1:nrow(dataquant),dataquant$logRNA/RNAfact,pch=16,col=rgb(1,0,0,0.2),cex=0.5)
posplot<-c()
GFPGOI<-c()
RNAGOI<-c()
for(i in 1:nrow(GOI)) {
  posplot[i]<-(1:nrow(dataquant))[dataquant$sysname==GOI$name[i]]
  GFPGOI[i]<-dataquant$GFP[dataquant$sysname==GOI$name[i]]
  RNAGOI[i]<-dataquant$logRNA[dataquant$sysname==GOI$name[i]]
}
GOI$posplot<-posplot
GOI$GFP<-GFPGOI
GOI$RNA<-RNAGOI
GOI<-GOI[order(posplot),]

thr<-log2(GOI$GFP[GOI$gene=="GPD1"]*1.5)-base
#abline(h=thr)
thr2<-(GOI$RNA[GOI$gene=="GPD1"]+log2(1.5))/RNAfact
thr2<-ACT1RNAthr/RNAfact
abline(h=thr2,col="red")
points((1:nrow(dataquant))[dataquant$logRNA/RNAfact>thr2],(dataquant$logRNA/RNAfact)[dataquant$logRNA/RNAfact>thr2],pch=16,col=rgb(1,0,0,1),cex=0.7)

for(i in 1:nrow(GOI)) {
  points(c(GOI$posplot[i],GOI$posplot[i]),c(0,log2(GOI$GFP[i])-base),typ="l",col="blue")
  text(GOI$posplot[i],rep(c(0.3,0.5),12)[i],cex=0.5,labels = GOI$gene[i])
}
points(GOI$posplot,GOI$RNA/RNAfact,pch=16,col=rgb(0,0,1,1),cex=1)
axis(1,c(0,n_noquant,n_noquant+n_toolow,nrow(dataquant)),labels = c(0,n_noquant,n_toolow,nrow(dataquant)-n_noquant-n_toolow))
axis(2,log2(c(20,80,300,1300,5000,20000,80000))-base,labels = c("20","80","300","1.3e3","5e3","20e3","80e3"))
axis(4,log2(c(1,20,250,3000,50000))/RNAfact,labels = c("1","20","250","3e3","50e3"))
dev.off()

#quantification
2^(thr2*RNAfact)
nRNAhight<-sum(dataquant$logRNA>=thr2*RNAfact,na.rm = T)
nRNAhight/nrow(dataquant)*100
n_toolow/(nrow(dataquant)-n_noquant)*100


##==============================
#investing link between UTR size and expression level (not used in the paper)

library(Biostrings)
SCgenome <- readDNAStringSet("sacCer3.fa")
nameCHR<-c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrM")

DATA4<-merge(DATA3,UTRdata[,-1],by.x=1,by.y=1,all.x=T,all.y=F)

i=29
UTR3min<-c()
UTR3max<-c()
GC<-c()
for (i in 1:nrow(DATA4)) {
  if (is.na(DATA4$UTR3[i])) {
    UTR3min[i]<-NA
    UTR3max[i]<-NA
    GC[i]<-NA
  } else if (DATA4$UTR3[i]==0) {
    UTR3min[i]<-NA
    UTR3max[i]<-NA
    GC[i]<-NA
  } else {
    if (DATA4$strand[i]=="W") {
      UTR3min[i]<-DATA4$End[i]+1
      UTR3max[i]<-DATA4$End[i]+DATA4$UTR3[i]
    } else {
      UTR3min[i]<-DATA4$End[i]-DATA4$UTR3[i]
      UTR3max[i]<-DATA4$End[i]-1
    }
    UTR3seq<-subseq(SCgenome[nameCHR[DATA4$chr[i]]], UTR3min[i], UTR3max[i])[[1]]
    GC[i]<-sum(alphabetFrequency(UTR3seq)[c("G","C")])/sum(alphabetFrequency(UTR3seq))
  }
}

size<-DATA4$posdown-DATA4$posup
sizel<-log2(size)

expYl<-log2(DATA4$meanexp)
#expYl<-DATA4$logmeanRNA

UTRlog<-DATA4$UTR3
UTRlog[UTRlog<10]<-NA
UTRlog<-log2(UTRlog)
plot(expYl,UTRlog,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(UTRlog~expYl)
abline(reg)
summary(reg)

plot(GC,UTRlog,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(UTRlog~GC)
abline(reg)
summary(reg)

plot(GC,expYl,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(expYl~GC)
abline(reg)
summary(reg)

plot(sizel,expYl,pch=16,cex=0.6,col=rgb(0,0,0,0.15))
reg<-lm(expYl[expYl>5]~sizel[expYl>5])
abline(reg)
summary(reg)
expYl_ZCorr<-expYl-((sizel*reg$coefficients[2])+reg$coefficients[1])+mean(expYl,na.rm=T)
plot(sizel,expYl_ZCorr,pch=16,cex=0.6,col=rgb(0,0,1,0.1))

plot(sizel,UTRlog,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(UTRlog~sizel)
abline(reg)
summary(reg)

plot(GC,expYl_ZCorr,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(expYl_ZCorr~GC)
abline(reg)
summary(reg)
expYl_ZCorr_CGcorr<-expYl_ZCorr-((GC*reg$coefficients[2])+reg$coefficients[1])+mean(expYl_ZCorr,na.rm=T)

plot(expYl_ZCorr,UTRlog,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(UTRlog~expYl_ZCorr)
abline(reg)
summary(reg)

plot(UTRlog,expYl_ZCorr,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(expYl_ZCorr~UTRlog)
abline(reg)
summary(reg)

reg<-lm(UTRlog~GC+expYl_ZCorr+GC*expYl_ZCorr)
summary(reg)

plot(expYl_ZCorr_CGcorr,UTRlog,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(UTRlog~expYl_ZCorr_CGcorr)
abline(reg)
summary(reg)

plot(expYl_ZCorr,UTRlog,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(UTRlog~expYl_ZCorr)
abline(reg)
summary(reg)

plot(expYl,UTRlog,pch=16,cex=0.6,col=rgb(0,0,0,0.2))
reg<-lm(UTRlog~expYl_ZCorr)
abline(reg)
summary(reg)

reg<-lm(UTRlog~GC+expYl_ZCorr_CGcorr+GC*expYl_ZCorr_CGcorr)
summary(reg)


#============================
#preparing final table for export

DATAf<-merge(genesINFO2,GFPdata,by.x=1,by.y=1,all.x=T,all.y=F)
notther<-GFPdata[!(GFPdata$ORF%in%genesINFO2$sysname),]
for (i in 1:nrow(notther)) {
  if (length(grep(notther$ORF[i],DATAf$alias))>0) {
    n<-(1:nrow(DATAf))[grepl(notther$ORF[i],DATAf$alias)]
    if (length(n)!=1) {print(n)}
    DATAf[n,ncol(DATAf)-((ncol(notther)-1):0)]<-notther[i,2:ncol(notther)]
  }
}

DATAf2<-merge(DATAf,expDATA2,by.x=1,by.y=1,all.x=T,all.y=F)
notther2<-expDATA2[!(expDATA2$Gene%in%DATAf$sysname),]

DATAf3<-merge(DATAf2,expDATAFrank,by.x=1,by.y=1,all.x=T,all.y=F)
notther2<-expDATAFrank[!(expDATAFrank$ORF%in%DATAf2$sysname),]

keepcol<-colnames(DATAf3)[c(1,2,3,4,5,6,9,10,12,14,24,25,31)]

DATAf4<-DATAf3[,colnames(DATAf3)%in%keepcol]
DATAf4<-DATAf4[,c(1,2,3,4,5,6,7,13,9,11,10,12,8)]
colnames(DATAf4)

write.table(DATAf4,file = "expre_mRNA_GFP_allgene.txt",quote = F,sep = "\t",na = "",col.names = T,row.names = F)
