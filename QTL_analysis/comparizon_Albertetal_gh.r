library(plyr)
library(reshape2)


alignmentDir <- "C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/190208-AllQTLprocess/"
setwd(alignmentDir)


compa_GFPdf<-read.table(paste(alignmentDir,"compa_GFPdf.txt",sep=""), stringsAsFactors=FALSE, head=TRUE, na.strings = "")
compa_mCHdf<-read.table(paste(alignmentDir,"compa_mCHdf.txt",sep=""), stringsAsFactors=FALSE, head=TRUE, na.strings = "")


#GFP comparizon
min(compa_GFPdf$albertLOD[compa_GFPdf$albertLOD>0])
summary(as.factor(compa_GFPdf$correspond))
compa_GFPdftrue<-compa_GFPdf[!(compa_GFPdf$names %in% c("RPS10A","CTS1","NA","TDH3")),]
summary(as.factor(compa_GFPdftrue$correspond))
compa_GFPdftrue45<-compa_GFPdftrue[(compa_GFPdftrue$brionLOD>4.5 | compa_GFPdftrue$albertLOD>4.5),]
summary(as.factor(compa_GFPdftrue45$correspond))
gfpsimBrion_Albert<-sum(abs(compa_GFPdftrue45$briondeltaAF[compa_GFPdftrue45$correspond=="together"]))/
  sum(abs(compa_GFPdftrue45$briondeltaAF[compa_GFPdftrue45$correspond%in%c("together","oposite","unique_brion")]))
gfpsimAlbert_Brion<-sum(abs(compa_GFPdftrue45$albertdeltaAF[compa_GFPdftrue45$correspond=="together"]))/
  sum(abs(compa_GFPdftrue45$albertdeltaAF[compa_GFPdftrue45$correspond%in%c("together","oposite","unique_albert")]))
sum(sign(compa_GFPdftrue45$albertdeltaAF)==sign(compa_GFPdftrue45$briondeltaAF))/nrow(compa_GFPdftrue45)
boxplot(compa_GFPdftrue45$brionLOD~compa_GFPdftrue45$correspond,ylim=c(0,100))
boxplot(compa_GFPdftrue45$albertLOD~compa_GFPdftrue45$correspond,ylim=c(0,100))

compa_GFPdftrue45$cexall<-compa_GFPdftrue45$brionLOD
compa_GFPdftrue45$pchall<-compa_GFPdftrue45$brionLOD
for (i in 1:nrow(compa_GFPdftrue45)) {
  if (compa_GFPdftrue45$correspond[i] %in% c("together","oposite")) {
    x<-max(compa_GFPdftrue45$brionLOD[i],compa_GFPdftrue45$albertLOD[i])
    compa_GFPdftrue45$pchall[i]<-16
  } else if (compa_GFPdftrue45$correspond[i] %in% c("unique_albert")) {
    x<-compa_GFPdftrue45$albertLOD[i]
    compa_GFPdftrue45$pchall[i]<-1
  } else if (compa_GFPdftrue45$correspond[i] %in% c("unique_brion")) {
    x<-compa_GFPdftrue45$brionLOD[i]
    compa_GFPdftrue45$pchall[i]<-1
  }
  compa_GFPdftrue45$cexall[i]<-max(0.4,min(3,x/20))
}

svg(filename = "GFPvsAlbert.svg",width = 5.5,height = 5.7)
plot(compa_GFPdftrue45$albertdeltaAF,compa_GFPdftrue45$briondeltaAF,cex=compa_GFPdftrue45$cexall,xlim=c(-1,1),ylim=c(-1,1),col="#33ab82ff",pch=compa_GFPdftrue45$pchall)
points(compa_GFPdftrue45$albertdeltaAF[compa_GFPdftrue45$pchall==16],compa_GFPdftrue45$briondeltaAF[compa_GFPdftrue45$pchall==16],cex=compa_GFPdftrue45$cexall[compa_GFPdftrue45$pchall==16],xlim=c(-1,1),ylim=c(-1,1),col="black",pch=1)
abline(h=0)
abline(v=0)
abline(a=0,b=1)
reg<-lm(compa_GFPdftrue45$briondeltaAF~compa_GFPdftrue45$albertdeltaAF)
abline(reg,col="#33ab82ff",lty=2)
summary(reg)
cor(compa_GFPdftrue45$albertdeltaAF,compa_GFPdftrue45$briondeltaAF,method = "p")
dev.off()


compa_GFPdftrue2<-compa_GFPdf[!(compa_GFPdf$names %in% c("NA","TDH3")),]
summary(as.factor(compa_GFPdftrue$correspond))
compa_GFPdftrue245<-compa_GFPdftrue2[(compa_GFPdftrue2$brionLOD>4.5 | compa_GFPdftrue2$albertLOD>4.5),]
summary(as.factor(compa_GFPdftrue245$correspond))

write.table(x = compa_GFPdftrue245,file = "processGFPcompa.txt",quote = F,sep = "\t",row.names = F,col.names = T)


min(compa_mCHdf$albertLOD[compa_mCHdf$albertLOD>0])
summary(as.factor(compa_mCHdf$correspond))
compa_mCHdftrue<-compa_mCHdf[!(compa_mCHdf$names %in% c("NA","TDH3")),]
summary(as.factor(compa_mCHdftrue$correspond))
compa_mCHdftrue45<-compa_mCHdftrue[(compa_mCHdftrue$brionLOD>4.5 | compa_mCHdftrue$albertLOD>2.8),]
summary(as.factor(compa_mCHdftrue45$correspond))
mchsimBrion_Albert<-sum(abs(compa_mCHdftrue45$briondeltaAF[compa_mCHdftrue45$correspond=="together"]))/
  sum(abs(compa_mCHdftrue45$briondeltaAF[compa_mCHdftrue45$correspond%in%c("together","oposite","unique_brion")]))
mchsimAlbert_Brion<-sum(abs(compa_mCHdftrue45$albertR[compa_mCHdftrue45$correspond=="together"]))/
  sum(abs(compa_mCHdftrue45$albertR[compa_mCHdftrue45$correspond%in%c("together","oposite","unique_albert")]))
sum(sign(compa_mCHdftrue45$albertR)==sign(compa_mCHdftrue45$briondeltaAF))/nrow(compa_mCHdftrue45)
boxplot(compa_mCHdftrue45$brionLOD~compa_mCHdftrue45$correspond,ylim=c(0,100))
boxplot(compa_mCHdftrue45$albertLOD~compa_mCHdftrue45$correspond,ylim=c(0,100))

compa_mCHdftrue45$cexall<-compa_mCHdftrue45$brionLOD
compa_mCHdftrue45$pchall<-compa_mCHdftrue45$brionLOD
for (i in 1:nrow(compa_mCHdftrue45)) {
  if (compa_mCHdftrue45$correspond[i] %in% c("together","oposite")) {
    x<-max(compa_mCHdftrue45$brionLOD[i],compa_mCHdftrue45$albertLOD[i])
    compa_mCHdftrue45$pchall[i]<-16
  } else if (compa_mCHdftrue45$correspond[i] %in% c("unique_albert")) {
    x<-compa_mCHdftrue45$albertLOD[i]
    compa_mCHdftrue45$pchall[i]<-1
  } else if (compa_mCHdftrue45$correspond[i] %in% c("unique_brion")) {
    x<-compa_mCHdftrue45$brionLOD[i]
    compa_mCHdftrue45$pchall[i]<-1
  }
  compa_mCHdftrue45$cexall[i]<-max(0.4,min(3,x/20))
}

plot(compa_mCHdftrue45$albertR,compa_mCHdftrue45$briondeltaAF,cex=compa_mCHdftrue45$cexall,xlim=c(-1,1),ylim=c(-1,1),col="red",pch=compa_mCHdftrue45$pchall)
points(compa_mCHdftrue45$albertR[compa_mCHdftrue45$pchall==16],compa_mCHdftrue45$briondeltaAF[compa_mCHdftrue45$pchall==16],cex=compa_mCHdftrue45$cexall[compa_mCHdftrue45$pchall==16],col="black",pch=1)
abline(h=0)
abline(v=0)
abline(a=0,b=1)
reg<-lm(compa_mCHdftrue45$briondeltaAF~compa_mCHdftrue45$albertR)
abline(reg,col="red",lty=2)
summary(reg)


compa_mCHdftrue45_14<-compa_mCHdftrue45[!(compa_mCHdftrue45$chr==14 & abs(compa_mCHdftrue45$brionmaxpos-450000)<100000),]
#compa_mCHdftrue45_14<-compa_mCHdftrue45[!(compa_mCHdftrue45$chr==14),]
summary(as.factor(compa_mCHdftrue45_14$correspond))
mchsimBrion_Albert_14<-sum(abs(compa_mCHdftrue45_14$briondeltaAF[compa_mCHdftrue45_14$correspond=="together"]))/
  sum(abs(compa_mCHdftrue45_14$briondeltaAF[compa_mCHdftrue45_14$correspond%in%c("together","oposite","unique_brion")]))
mchsimAlbert_Brion_14<-sum(abs(compa_mCHdftrue45_14$albertR[compa_mCHdftrue45_14$correspond=="together"]))/
  sum(abs(compa_mCHdftrue45_14$albertR[compa_mCHdftrue45_14$correspond%in%c("together","oposite","unique_albert")]))
sum(sign(compa_mCHdftrue45_14$albertR)==sign(compa_mCHdftrue45_14$briondeltaAF))/nrow(compa_mCHdftrue45_14)
boxplot(compa_mCHdftrue45_14$brionLOD~compa_mCHdftrue45_14$correspond,ylim=c(0,100))
boxplot(compa_mCHdftrue45_14$albertLOD~compa_mCHdftrue45_14$correspond,ylim=c(0,100))


svg(filename = "mCHvsAlbert.svg",width = 5.5,height = 5.7)
plot(compa_mCHdftrue45_14$albertR,compa_mCHdftrue45_14$briondeltaAF,cex=compa_mCHdftrue45_14$cexall,xlim=c(-1,1),ylim=c(-1,1),col="red",pch=compa_mCHdftrue45_14$pchall)
points(compa_mCHdftrue45_14$albertR[compa_mCHdftrue45_14$pchall==16],compa_mCHdftrue45_14$briondeltaAF[compa_mCHdftrue45_14$pchall==16],cex=compa_mCHdftrue45_14$cexall[compa_mCHdftrue45_14$pchall==16],col="black",pch=1)
abline(h=0)
abline(v=0)
abline(a=0,b=1)
reg<-lm(compa_mCHdftrue45_14$briondeltaAF~compa_mCHdftrue45_14$albertR)
abline(reg,col="red",lty=2)
summary(reg)
cor(compa_mCHdftrue45_14$albertR,compa_mCHdftrue45_14$briondeltaAF,method = "p")
compamCHbrionTrueFALSE<-compa_mCHdftrue45[(compa_mCHdftrue45$chr==14 & abs(compa_mCHdftrue45$brionmaxpos-450000)<100000),]
#compamCHbrionTrueFALSE<-compa_mCHdftrue45[(compa_mCHdftrue45$chr==14),]
points(compamCHbrionTrueFALSE$albertR,compamCHbrionTrueFALSE$briondeltaAF,pch=1,cex=compamCHbrionTrueFALSE$cexall,col="grey")
dev.off()

write.table(x = compa_mCHdftrue45,file = "processmCHcompa.txt",quote = F,sep = "\t",row.names = F,col.names = T)
