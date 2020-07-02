setwd("C:/Users/brion/Dropbox/Gdrive/MN_postdoc/diverR/200131-FlowSegPlot")
heri<-read.table("Heritability.txt",header = T,sep = "\t", quote = "", na.strings = "",stringsAsFactors = F)
heri<-read.table("HeritabilitySc.txt",header = T,sep = "\t", quote = "", na.strings = "",stringsAsFactors = F)
heriD<-read.table("HeritabilityScD.txt",header = T,sep = "\t", quote = "", na.strings = "",stringsAsFactors = F)
correl<-read.table("CorrelFluo.txt",header = T,sep = "\t", quote = "", na.strings = "",stringsAsFactors = F)

genes<-levels(as.factor(correl$strains))
genesO<-genes[c(1,2,5,7,12,4,9,11,3,10,8,6)]

correlF<-correl[1:12,]
for (i in 1:length(genesO)) {
  correlF$strains[i]<-genesO[i]
  correlF$correl[i]<-mean(correl$correl[correl$strains==genesO[i]])
  correlF$r2[i]<-mean(correl$r2[correl$strains==genesO[i]])
  correlF$correlc[i]<-mean(correl$correlc[correl$strains==genesO[i]])
  correlF$r2c[i]<-mean(correl$r2c[correl$strains==genesO[i]])
}
heriF<-heri[1:44,]
heriDF<-heriD[1:22,]
for (i in 1:length(genesO[-12])) {
  for (j in 1:2) {
    tempheri<-heri[heri$chanel==c("mCH","GFP")[j],]
    for (k in 1:2) {
      n<-((j-1)*22)+((k-1)*11)+i
      temp2heri<-tempheri[tempheri$pop==c("high","low")[k],]
      heriF$GOI[n]<-genesO[i]
      heriF$chanel[n]<-c("mCH","GFP")[j]
      heriF$pop[n]<-c("high","low")[k]
      heriF$mCHeffect[n]<-mean(temp2heri$mCHeffect[temp2heri$GOI==genesO[i]])
      heriF$mCHpval[n]<-mean(temp2heri$mCHpval[temp2heri$GOI==genesO[i]])
      heriF$GFPeffect[n]<-mean(temp2heri$GFPeffect[temp2heri$GOI==genesO[i]])
      heriF$GFPpval[n]<-mean(temp2heri$GFPpval[temp2heri$GOI==genesO[i]])
    }
    n<-((j-1)*11)+i
    temp2heri<-heriD[heriD$chanel==c("mCH","GFP")[j],]
    heriDF$GOI[n]<-genesO[i]
    heriDF$chanel[n]<-c("mCH","GFP")[j]
    heriDF$mCHeffect[n]<-mean(temp2heri$mCHeffect[temp2heri$GOI==genesO[i]])
    heriDF$mCHpval[n]<-mean(temp2heri$mCHpval[temp2heri$GOI==genesO[i]])
    heriDF$GFPeffect[n]<-mean(temp2heri$GFPeffect[temp2heri$GOI==genesO[i]])
    heriDF$GFPpval[n]<-mean(temp2heri$GFPpval[temp2heri$GOI==genesO[i]])
  }
}
heri$pch_mch<-rep(16,nrow(heri))
heri$pch_mch[heri$mCHpval>0.00001]<-1
heri$pch_GFP<-rep(16,nrow(heri))
heri$pch_GFP[heri$GFPpval>0.00001]<-1
heriD$pch_mch<-rep(16,nrow(heriD))
heriD$pch_mch[heriD$mCHpval>0.00001]<-1
heriD$pch_GFP<-rep(16,nrow(heriD))
heriD$pch_GFP[heriD$GFPpval>0.00001]<-1

svg(filename = "correl_r2.svg",width = 8,height = 5.5)
h<-barplot(sqrt(correlF$r2),ylim=c(0,1),col=rgb(0.5,0.6,0.75,0.2))
for (i in 1:12) {
  tempcorrelplot<-correl[correl$strains==genesO[i],]
  points(rep(h[i],nrow(tempcorrelplot)),sqrt(tempcorrelplot$r2),pch=16,col=rgb(0.5,0.6,0.75,0.5))
}
dev.off()
svg(filename = "correl_r2c.svg",width = 8,height = 5.5)
h<-barplot(sqrt(correlF$r2c),ylim=c(0,1),col=rgb(0.7,0.7,0.7,1))
for (i in 1:12) {
  tempcorrelplot<-correl[correl$strains==genesO[i],]
  points(rep(h[i],nrow(tempcorrelplot)),sqrt(tempcorrelplot$r2c),pch=16,col=rgb(0.1,0.1,0.1,0.7))
}
dev.off()

chan<-c(rep("mCH",4),rep("GFP",4))
fluo<-c(rep(c("mCHeffect","GFPeffect"),4))
gate<-c(rep("high",2),rep("low",2),rep("high",2),rep("low",2))
colplot<-c(rep(c("red","#34ab83ff"),4))
j=1
for (j in 1:8) {
  coluValue<-(1:ncol(heriF))[colnames(heriF)==fluo[j]]
  if (fluo[j]=="mCHeffect") {pch_temp<-13} else {pch_temp<-14}
  #svg(filename = paste("heri",chan[j],fluo[j],gate[j],".svg",sep="_"),width = 5,height = 5.5)
  h<-barplot(heriF[heriF$chanel==chan[j] & heriF$pop==gate[j],coluValue],ylim=c(-1,1),col=colplot[j])
  for (i in 1:11) {
    tempheriplot<-heri[heri$chanel==chan[j] & heri$pop==gate[j] & heri$GOI==genesO[i],]
    points(rep(h[i],nrow(tempheriplot)),tempheriplot[,coluValue],pch=tempheriplot[,pch_temp])
  }
  #dev.off()
}

for (j in 1:8) {
  coluValue<-(1:ncol(heriF))[colnames(heriF)==fluo[j]]
  if (fluo[j]=="mCHeffect") {pch_temp<-13} else {pch_temp<-14}
  svg(filename = paste("heri",chan[j],fluo[j],gate[j],".svg",sep="_"),width = 5,height = 5.5)
  h<-barplot(sqrt(abs(heriF[heriF$chanel==chan[j] & heriF$pop==gate[j],coluValue]))*sign(heriF[heriF$chanel==chan[j] & heriF$pop==gate[j],coluValue])
             ,ylim=c(-1,1),col=colplot[j],yaxt="n")
  for (i in 1:11) {
    tempheriplot<-heri[heri$chanel==chan[j] & heri$pop==gate[j] & heri$GOI==genesO[i],]
    points(rep(h[i],nrow(tempheriplot)),sqrt(abs(tempheriplot[,coluValue]))*sign(tempheriplot[,coluValue]),pch=tempheriplot[,pch_temp])
  }
  axis(side = 2,at = c(sqrt(c(c(1:4,6:9,11:14,16:20)*0.05)),-sqrt(c(c(1:4,6:9,11:14,16:20)*0.05))),labels = rep(NA,20),col.ticks = rgb(0,0,0,0.3))
  axis(side = 2,at = c(sqrt(c((1:4)*0.25)),0,-sqrt(c((1:4)*0.25))),labels = rep(NA,20),col.ticks = rgb(0,0,0,1))
  dev.off()
}



chan<-c(rep("mCH",2),rep("GFP",2))
fluo<-c(rep(c("mCHeffect","GFPeffect"),2))
colplot<-c(rep(c("red","#34ab83ff"),2))
j=1
for (j in 1:4) {
  coluValue<-(1:ncol(heriDF))[colnames(heriDF)==fluo[j]]
  if (fluo[j]=="mCHeffect") {pch_temp<-13} else {pch_temp<-14}
  #svg(filename = paste("heridetla",chan[j],fluo[j],".svg",sep="_"),width = 7,height = 5.5)
  h<-barplot(heriDF[heriDF$chanel==chan[j],coluValue],ylim=c(-1.5,1.5),col=colplot[j])
  for (i in 1:11) {
    tempheriplot<-heriD[heriD$chanel==chan[j] & heriD$GOI==genesO[i],]
    points(rep(h[i],nrow(tempheriplot)),tempheriplot[,coluValue],pch=tempheriplot[,pch_temp])
  }
  #dev.off()
}

for (j in 1:4) {
  coluValue<-(1:ncol(heriDF))[colnames(heriDF)==fluo[j]]
  if (fluo[j]=="mCHeffect") {pch_temp<-13} else {pch_temp<-14}
  svg(filename = paste("heridetla",chan[j],fluo[j],".svg",sep="_"),width = 7,height = 5.5)
  h<-barplot(sqrt(abs(heriDF[heriDF$chanel==chan[j],coluValue]))*sign(heriDF[heriDF$chanel==chan[j],coluValue])
             ,ylim=c(-sqrt(0.25),sqrt(1.5)),col=colplot[j],yaxt="n")
  for (i in 1:11) {
    tempheriplot<-heriD[heriD$chanel==chan[j] & heriD$GOI==genesO[i],]
    points(rep(h[i],nrow(tempheriplot)),sqrt(abs(tempheriplot[,coluValue]))*sign(tempheriplot[,coluValue]),pch=tempheriplot[,pch_temp])
  }
  axis(side = 2,at = c(sqrt(c(c(1:4,6:9,11:14,16:19,21:24,26:29)*0.05)),-sqrt(c(c(1:4)*0.05))),labels = rep(NA,20),col.ticks = rgb(0,0,0,0.3))
  axis(side = 2,at = c(sqrt(c((1:6)*0.25)),0,-sqrt(0.25)),labels = rep(NA,20),col.ticks = rgb(0,0,0,1))
  dev.off()
}

