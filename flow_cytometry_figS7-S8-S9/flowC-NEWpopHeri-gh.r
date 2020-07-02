library(flowCore)
library(reshape2)
library(plyr)
library(data.table)
par(mfrow = c(1, 1))

setwd("YOUR/WORKING/DIRECTORY")
file<-all[1]


infoplate<-read.csv("181119matSample.txt", sep="\t",header = T)

import_dataO <- function(file){
  x = unlist(strsplit(basename(file), "\\."))[1]
  name = unlist(strsplit(basename(x), "\\_"))[4]
  #print(paste("Reading", condition, id, time, sep=" "))  
  name <- as.numeric(name)
  data <- read.FCS(file, transformation=FALSE)
  print(paste("Reading", x, nrow(data), sep=" "))  
  mdata <- data.frame(exprs(data))
  mdata$gene <- infoplate$gene[infoplate$sample==name]
  mdata$rep <- infoplate$rep[infoplate$sample==name]
  mdata$tec <- infoplate$teck[infoplate$sample==name]
  mdata$gate <- infoplate$gate[infoplate$sample==name]
  mdata$sample <- infoplate$sample[infoplate$sample==name]
  return(mdata)
}

all = list.files(path="sort181119", pattern = "Specimen", full = TRUE)
data_all = ldply(all, import_dataO) #plyr library
colnames(data_all)[colnames(data_all)=="eGFP.A"]<-"GFP.A"

data_all$strain<-paste(data_all$gene,data_all$rep,sep="_")
data_all$name<-paste(data_all$gene,data_all$rep,data_all$tec,data_all$gate,sep="_")

boxplot(log(data_all$SSC.A)~data_all$name)


#only single
randdata40000<-data_all[sample(nrow(data_all),40000,replace = FALSE),]

plot(log(randdata40000$FSC.A),log(randdata40000$SSC.A),pch=16,col=rgb(0,0,1,0.05),cex=0.4,ylim=c(9,12.5),xlim=c(10,12.5))
gateTOPx<-c(10.8,10.9,11.8,12.3)
gateTOPy<-c(9.8,10.2,11.2,11.2)
gateBOTx<-c(10.8,11.2,12,12.3)
gateBOTy<-c(9.8,9.8,10.66,11.2)
points(c(gateTOPx,rev(gateBOTx)),c(gateTOPy,rev(gateBOTy)),typ="l",col="grey",lwd=2)


subdata<-data_all
for (i in 1:(length(gateTOPx)-1)) {
  coor<-lm(gateTOPy[i:(i+1)]~gateTOPx[i:(i+1)])
  subdata<-subdata[log(subdata$SSC.A) < coor$coefficients[2]*log(subdata$FSC.A)+coor$coefficients[1],]
}
for (i in 1:(length(gateBOTx)-1)) {
  coor<-lm(gateBOTy[i:(i+1)]~gateBOTx[i:(i+1)])
  subdata<-subdata[log(subdata$SSC.A) > coor$coefficients[2]*log(subdata$FSC.A)+coor$coefficients[1],]
}
subdata<-subdata[subdata$FSC.A>0 & subdata$SSC.A>0 & subdata$mCherry.A>0 & subdata$GFP.A>0,]
subdata<-subdata[!is.na(subdata$gene),]

genes<-levels(as.factor(subdata$gene))
genes<-genes[c(5,6,10,1,7,3,8,2,11)]
typetest<-c("whole pop mCH","whole pop GFP")
mat<-matrix(c(1,2,3,1,4,5),nrow = 2,ncol = 3,byrow = T)

subdata$mCHc<-subdata$mCherry.A
subdata$GFPc<-subdata$mCherry.A


strains<-levels(as.factor(subdata$strain))
strains <- strains[c(4:7,8,12,1,9,3,11,2,10)]

strain <- strains[1]
plot(log(subdata$mCherry.A[subdata$strain==strain])~subdata$FSC.A[subdata$strain==strain])
#test<-loess(log(subdata$mCherry.A[subdata$strain==strain])~subdata$FSC.A[subdata$strain==strain]) #15s
#plot(unlist(test$residuals)~subdata$FSC.A[subdata$strain==strain])

plot(log(subdata$GFP.A[subdata$strain==strain])~subdata$FSC.A[subdata$strain==strain])
#test2<-loess(log(subdata$GFP.A[subdata$strain==strain])~subdata$FSC.A[subdata$strain==strain])
#plot(unlist(test2$residuals)~subdata$FSC.A[subdata$strain==strain])



correl<-c()
r2<-c()
correlc<-c()
r2c<-c()
strains<-levels(as.factor(subdata$strain))
strains <- strains[c(4:7,8,12,1,9,3,11,2,10)]

pdf("181120_Flow_effectSize.pdf", width = 12, height = 12)
par(mfcol = c(2, 2))
i=1
date()
#subdataG1<-subdata[subdata$gate==1,]
for (i in 1:length(strains)) {
  strain<-strains[i]
  test<-loess(log(subdata$mCherry.A[subdata$strain==strain])~subdata$FSC.A[subdata$strain==strain])
  subdata$mCHc[subdata$strain==strain]<-unlist(test$residuals)+mean(log(subdata$mCherry.A[subdata$strain==strain]))
  test2<-loess(log(subdata$GFP.A[subdata$strain==strain])~subdata$FSC.A[subdata$strain==strain])
  subdata$GFPc[subdata$strain==strain]<-unlist(test2$residuals)+mean(log(subdata$GFP.A[subdata$strain==strain]))
  reg<-lm(log(subdata$mCherry.A[subdata$strain==strain])~log(subdata$GFP.A[subdata$strain==strain]))
  reg2<-lm(unlist(test$residuals)~unlist(test2$residuals))
  
  plot(log(subdata$mCherry.A[subdata$strain==strain])~log(subdata$GFP.A[subdata$strain==strain]),pch=16,cex=0.6,col=rgb(0,0,0,0.12),xlim=c(2,12),ylim=c(4,10),main=paste(strain),xlab = "GFP",ylab = "mCherry")
  reg<-lm(log(subdata$mCherry.A[subdata$strain==strain])~log(subdata$GFP.A[subdata$strain==strain]))
  summary(reg)
  abline(reg,col="red")
  text(x = 3, y = 9, paste("s:",round(reg$coefficients[[2]],3)))
  text(x = 3, y = 8, paste("r:",round(summary(reg)$r.squared,3)))
  
  plot(subdata$mCHc[subdata$strain==strain]~subdata$GFPc[subdata$strain==strain],pch=16,cex=0.6,col=rgb(0,0,0,0.12),xlim=c(2,12),ylim=c(4,10),main=paste(strain,"corrected"),xlab = "GFP",ylab = "mCherry")
  reg2<-lm(subdata$mCHc[subdata$strain==strain]~subdata$GFPc[subdata$strain==strain])
  summary(reg2)
  abline(reg2,col="red")
  text(x = 3, y = 9, paste("s:",round(reg2$coefficients[[2]],3)))
  text(x = 3, y = 8, paste("r2:",round(summary(reg2)$r.squared,3)))
  
  
  correl[i]<-reg$coefficients[[2]]
  correlc[i]<-reg2$coefficients[[2]]
  r2[i]<-summary(reg)$r.squared
  r2c[i]<-summary(reg2)$r.squared
}
date() #10min

dev.off()


#===================
#plot for paper

setwd("/home/christian/Dropbox/Gdrive/MN_postdoc/diverR/forPaper")

subdataplot<-subdata[subdata$strain=="GPD1_1",]
subdataplot10000<-subdataplot[sample(nrow(subdataplot),2000,replace = FALSE),]

transp<-0.3
cexs<-0.7
svg(filename = "correlGPD1-raw.svg",width = 5.5,height = 5)
plot(log(subdataplot10000$mCherry.A)~log(subdataplot10000$GFP.A),pch=16,cex=cexs,col=rgb(0,0,0,transp),xlim=c(5,10),ylim=c(5,10),main=paste(strain),xlab = "GFP",ylab = "mCherry")
reg<-lm(log(subdataplot$mCherry.A)~log(subdataplot$GFP.A))
#summary(reg)
abline(reg,col="red")
text(x = 5.5, y = 9, paste("s:",round(reg$coefficients[[2]],3)))
text(x = 5.5, y = 8, paste("r2:",round(summary(reg)$r.squared,3)))
dev.off()
svg(filename = "loessGPD1-mCH-raw.svg",width = 4,height = 3.5)
plot((subdataplot10000$FSC.A/1000),log(subdataplot10000$mCherry.A),pch=16,cex=cexs,col=rgb(1,0,0,transp),xlim=c(50,200),ylim=c(5,10),main=paste(strain),xlab = "FSC.A",ylab = "mCherry")
dev.off()
svg(filename = "loessGPD1-mCH-correct.svg",width = 4,height = 3.5)
plot((subdataplot10000$FSC.A/1000),subdataplot10000$mCHc,pch=16,cex=cexs,col=rgb(1,0,0,transp),xlim=c(50,200),ylim=c(5,10),main=paste(strain),xlab = "FSC.A",ylab = "mCherry corrected")
dev.off()
svg(filename = "loessGPD1-GFP-raw.svg",width = 4,height = 3.5)
plot(subdataplot10000$FSC.A/1000,log(subdataplot10000$GFP.A),pch=16,cex=cexs,col=rgb(0.2,0.67,0.5,transp),xlim=c(50,200),ylim=c(5,10),main=paste(strain),xlab = "FSC.A",ylab = "GFP")
dev.off()
svg(filename = "loessGPD1-GFP-correct.svg",width = 4,height = 3.5)
plot(subdataplot10000$FSC.A/1000,subdataplot10000$GFPc,pch=16,cex=cexs,col=rgb(0.2,0.67,0.5,transp),xlim=c(50,200),ylim=c(5,10),main=paste(strain),xlab = "FSC.A",ylab = "GFP corrected")
dev.off()
svg(filename = "correlGPD1-correct.svg",width = 5.5,height = 5)
plot(subdataplot10000$mCHc~subdataplot10000$GFPc,pch=16,cex=cexs,col=rgb(0,0,0,transp),xlim=c(5,10),ylim=c(5,10),main=paste(strain,"corrected"),xlab = "GFP",ylab = "mCherry")
reg2<-lm(subdataplot$mCHc~subdataplot$GFPc)
#summary(reg2)
abline(reg2,col="red")
text(x = 5.5, y = 9, paste("s:",round(reg2$coefficients[[2]],3)))
text(x = 5.5, y = 8, paste("r2:",round(summary(reg2)$r.squared,3)))
dev.off()


transp<-0.3
cexs<-0.7
svg(filename = "correlGPD1-corfus.svg",width = 5.5,height = 5)
plot(log(subdataplot10000$mCherry.A)~log(subdataplot10000$GFP.A),pch=16,cex=cexs,col=rgb(0.5,0.6,0.75,transp),xlim=c(5,10),ylim=c(5,10),main=paste(strain),xlab = "GFP",ylab = "mCherry")
reg<-lm(log(subdataplot$mCherry.A)~log(subdataplot$GFP.A))
#summary(reg)
abline(reg,col=rgb(0.5,0.6,0.75))
points(subdataplot10000$mCHc~subdataplot10000$GFPc,pch=16,cex=cexs,col=rgb(0.2,0.2,0.2,transp),xlim=c(5,10),ylim=c(5,10),main=paste(strain,"corrected"),xlab = "GFP",ylab = "mCherry")
reg2<-lm(subdataplot$mCHc~subdataplot$GFPc)
#summary(reg2)
abline(reg2,col=rgb(0.2,0.2,0.2))
dev.off()

#=========================================
setwd("/home/christian/Dropbox/Gdrive/MN_postdoc/diverR")

sizecor<-data.frame(strains,correl,r2,correlc,r2c)
write.table(sizecor,"FLOWsizecor_new.txt",sep="\t",row.names = FALSE, quote = F)




GOI<-as.character(factor(rep(strains,4),levels=strains)[order(factor(rep(strains,4),levels=strains))])

final<-data.frame(GOI,stringsAsFactors = F)
final$chanel<-rep(c("mCH","mCH","GFP","GFP"),length(GOI)/4)
final$pop<-rep(c("high","low"),length(GOI)/4)
final$atlchanel<-rep("all",length(GOI))
final$mCHeffect<-rep(0,length(GOI))
final$mCHpval<-rep(0,length(GOI))
final$GFPeffect<-rep(0,length(GOI))
final$GFPpval<-rep(0,length(GOI))
final$finalPV<-rep(0,length(GOI))
final$heritability<-rep("",length(GOI))
final$synergy<-rep("",length(GOI))
finalc<-final
finalcD<-finalc[finalc$pop=="high",]

pdf("181120_Flow_sortAll.pdf", width = 22, height = 12)
par(mfcol = c(3, 4))

n=1
for (strain in strains) {
  plotdata<-subdata[subdata$strain==strain,]
  gene<-strain
  
  for (i in 1:2) {
    m<-mat[i,]
    plot(log(plotdata$GFP.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[1]]),pch=16,col=rgb(0,0,0,0.2),cex=0.4,xlim=c(2,12),ylim=c(2,12),xlab="GFP",ylab="mCherry",main=paste(gene,typetest[i],sep=" "))
    points(log(plotdata$GFP.A[plotdata$gate==m[2]]),log(plotdata$mCherry.A[plotdata$gate==m[2]]),pch=16,col=rgb(0,0,1,0.2),cex=0.4)
    points(log(plotdata$GFP.A[plotdata$gate==m[3]]),log(plotdata$mCherry.A[plotdata$gate==m[3]]),pch=16,col=rgb(0.9,0.3,0.1,0.2),cex=0.4)
    d1 <- density(log(plotdata$GFP.A[plotdata$gate==m[1]]),na.rm = T) 
    d2 <- density(log(plotdata$GFP.A[plotdata$gate==m[2]]),na.rm = T) # returns the density data
    d3 <- density(log(plotdata$GFP.A[plotdata$gate==m[3]]),na.rm = T) # returns the density data
    plot(d1, lwd = 2, main="", col = rgb(0,0,0,1), xlim=c(2,12),xlab="GFP",ylim=c(0,1.1*max(c(d1$y,d2$y,d3$y)))) # plots the results 
    points(d2,typ = "l", lwd = 2, col = rgb(0,0,1,1))
    points(d3,typ = "l", lwd = 2, col = rgb(0.9,0.3,0.1,1))
    testup<-t.test(log(plotdata$GFP.A[plotdata$gate==m[1]]),log(plotdata$GFP.A[plotdata$gate==m[3]]))
    testdown<-t.test(log(plotdata$GFP.A[plotdata$gate==m[1]]),log(plotdata$GFP.A[plotdata$gate==m[2]]))
    legend("topleft",legend = c(paste(typetest[i],"up",round(testup$estimate[2]-testup$estimate[1],digits = 3),round(-log10(testup$p.value),0),sep=" "),
                                paste(typetest[i],"down",round(testdown$estimate[2]-testdown$estimate[1],digits = 3),round(-log10(testdown$p.value),0),sep=" ")),
           lty=1,col=c(rgb(0.9,0.3,0.1,1),rgb(0,0,1,1)))
    final$GFPeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$estimate[2]-testup$estimate[1]
    final$GFPpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$p.value
    final$GFPeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$estimate[2]-testdown$estimate[1]
    final$GFPpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$p.value
    
    d1 <- density(log(plotdata$mCherry.A[plotdata$gate==m[1]]),na.rm = T) 
    d2 <- density(log(plotdata$mCherry.A[plotdata$gate==m[2]]),na.rm = T) # returns the density data
    d3 <- density(log(plotdata$mCherry.A[plotdata$gate==m[3]]),na.rm = T) # returns the density data
    plot(d1, lwd = 2, main="", col = rgb(0,0,0,1), xlim=c(2,12),xlab="mCherry",ylim=c(0,1.1*max(c(d1$y,d2$y,d3$y)))) # plots the results 
    points(d2,typ = "l", lwd = 2, col = rgb(0,0,1,1))
    points(d3,typ = "l", lwd = 2, col = rgb(0.9,0.3,0.1,1))
    testup<-t.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[3]]))
    testdown<-t.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[2]]))
    legend("topleft",legend = c(paste(typetest[i],"up",round(testup$estimate[2]-testup$estimate[1],digits = 3),round(-log10(testup$p.value),0),sep=" "),
                                paste(typetest[i],"down",round(testdown$estimate[2]-testdown$estimate[1],digits = 3),round(-log10(testdown$p.value),0),sep=" ")),
           lty=1,col=c(rgb(0.9,0.3,0.1,1),rgb(0,0,1,1)))
    final$mCHeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$estimate[2]-testup$estimate[1]
    final$mCHpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$p.value
    final$mCHeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$estimate[2]-testdown$estimate[1]
    final$mCHpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$p.value
    
    testup<-t.test(plotdata$GFPc[plotdata$gate==m[1]],plotdata$GFPc[plotdata$gate==m[3]])
    testdown<-t.test(plotdata$GFPc[plotdata$gate==m[1]],plotdata$GFPc[plotdata$gate==m[2]])
    finalc$GFPeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$estimate[2]-testup$estimate[1]
    finalc$GFPpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$p.value
    finalc$GFPeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$estimate[2]-testdown$estimate[1]
    finalc$GFPpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$p.value
    
    
    testup<-t.test(plotdata$mCHc[plotdata$gate==m[1]],plotdata$mCHc[plotdata$gate==m[3]])
    testdown<-t.test(plotdata$mCHc[plotdata$gate==m[1]],plotdata$mCHc[plotdata$gate==m[2]])
    finalc$mCHeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$estimate[2]-testup$estimate[1]
    finalc$mCHpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$p.value
    finalc$mCHeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$estimate[2]-testdown$estimate[1]
    finalc$mCHpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$p.value
    
    testup<-t.test(plotdata$GFPc[plotdata$gate==m[2]],plotdata$GFPc[plotdata$gate==m[3]])
    finalcD$GFPeffect[finalcD$GOI==gene & finalcD$chanel==final$chanel[i*2] & finalcD$pop=="high" & finalcD$atlchanel==final$atlchanel[i*2]]<-testup$estimate[2]-testup$estimate[1]
    finalcD$GFPpval[finalcD$GOI==gene & finalcD$chanel==final$chanel[i*2] & finalcD$pop=="high" & finalcD$atlchanel==final$atlchanel[i*2]]<-testup$p.value
    
    testup<-t.test(plotdata$mCHc[plotdata$gate==m[2]],plotdata$mCHc[plotdata$gate==m[3]])
    finalcD$mCHeffect[finalcD$GOI==gene & finalcD$chanel==final$chanel[i*2] & finalcD$pop=="high" & finalcD$atlchanel==final$atlchanel[i*2]]<-testup$estimate[2]-testup$estimate[1]
    finalcD$mCHpval[finalcD$GOI==gene & finalcD$chanel==final$chanel[i*2] & finalcD$pop=="high" & finalcD$atlchanel==final$atlchanel[i*2]]<-testup$p.value
    
  }
}
dev.off()


for (i in 1:nrow(final)) {
  if (final$pop[i]=="high") {m=1} else {m=-1}
  if (final$chanel[i]=="mCH") {
    cha<-5
    oth<-7
  } else {
    cha<-7
    oth<-5
  }
  final$finalPV[i]<-final[i,cha+1]
  if (final[i,cha]*m<0 | final[i,cha+1] > 0.01) {final$heritability[i]<-"nul"}
  else if (final[i,cha+1] > 0.001) {final$heritability[i]<-"average"}
  else if (final[i,cha+1] <= 0.001) {final$heritability[i]<-"good"}
  if (final[i,oth+1] > 0.01) {final$synergy[i]<-"none"}
  else if (final[i,oth]*m<0) {final$synergy[i]<-"inverted"}
  else if (final[i,oth]*m>0) {final$synergy[i]<-"similar"}
}

plot(final$mCHeffect,finalc$mCHeffect)
abline(a=0,b=1)
abline(v=0,col="grey")
abline(h=0,col="grey")
plot(final$GFPeffect,finalc$GFPeffect)
abline(a=0,b=1)
abline(v=0,col="grey")
abline(h=0,col="grey")
plot(log(-log10(final$mCHpval)),log(-log10(finalc$mCHpval)))
abline(a=0,b=1)
abline(v=log(-log10(0.00001)),col="red")
abline(h=log(-log10(0.00001)),col="red")
plot(log(-log10(final$GFPpval)),log(-log10(finalc$GFPpval)))
abline(a=0,b=1)
abline(v=log(-log10(0.00001)),col="red")
abline(h=log(-log10(0.00001)),col="red")


write.table(final,"finalFLOW_new.txt",sep="\t",row.names = FALSE, quote = F)
write.table(finalc,"finalFLOW_sc_new.txt",sep="\t",row.names = FALSE, quote = F)
write.table(finalcD,"finalFLOW_scD_new.txt",sep="\t",row.names = FALSE, quote = F)
#final2<-read.table("finalFLOW_new.txt",sep="\t",stringsAsFactors = F, header = T)
#finalF<-rbind(final,final2)
finalF<-final



barplot(as.matrix(dcast(finalF,heritability~chanel+finalF$GOI)[c(3,1,2),-1]))

coltyp<-c(rep(rgb(0,0.5,0.3,1),7),rep(rgb(1,0,0,1),7),rep(rgb(0,0.37,0.25,1),7),rep(rgb(0.6,0,0,1),7))
boxplot(-log10(finalF$finalPV)~finalF$GOI+finalF$chanel,ylim=c(0,300),col=coltyp)
points(c(0,29),c(3,3),typ="l",col="red")
points(c(0,29),c(2,2),typ="l",col="grey")

barplot(as.matrix(dcast(finalF,heritability~atlchanel+chanel)[c(3,1,2),-1]))

coltyp<-c(rep(rgb(0,0.5,0.3,1),2),rep(rgb(1,0,0,1),2),rep(rgb(0,0.37,0.25,1),2),rep(rgb(0.6,0,0,1),2))
boxplot(-log10(finalF$finalPV)~finalF$atlchanel+finalF$chanel,ylim=c(0,300),col=coltyp)
points(c(0.5,8.5),c(3,3),typ="l",col="red")
points(c(0.5,8.5),c(2,2),typ="l",col="grey")

coltyp<-c(rgb(0.75,0.45,0.15,1),rgb(1,1,1,1),rgb(0.2,0.8,1,1))
barplot(as.matrix(dcast(finalF,synergy~atlchanel+finalF$GOI)[,-1]),col=coltyp)

barplot(as.matrix(dcast(finalF,synergy~atlchanel)[,-1]),col=coltyp)


#===================
#plot for paper density

setwd("/home/christian/Dropbox/Gdrive/MN_postdoc/diverR/forPaper")

subdataplot<-subdata[subdata$strain=="GPD1_1",]

m<-mat[1,]
d1 <- density(log(subdataplot$GFP.A[subdataplot$gate==m[1]]),na.rm = T) 
d2 <- density(log(subdataplot$GFP.A[subdataplot$gate==m[2]]),na.rm = T) # returns the density data
d3 <- density(log(subdataplot$GFP.A[subdataplot$gate==m[3]]),na.rm = T) # returns the density data
svg(filename = "density-GPF1-GFP-mchpop.svg",width = 5,height = 3.5)
plot(d1, lwd = 2, main="", col = rgb(0,0,0,1), xlim=c(5,10),xlab="GFP",ylim=c(0,1.1*max(c(d1$y,d2$y,d3$y)))) # plots the results 
points(d2,typ = "l", lwd = 2, col = rgb(0,0,1,1))
points(d3,typ = "l", lwd = 2, col = rgb(0.9,0.3,0.1,1))
testup<-t.test(log(subdataplot$GFP.A[subdataplot$gate==m[1]]),log(subdataplot$GFP.A[subdataplot$gate==m[3]]))
testdown<-t.test(log(subdataplot$GFP.A[subdataplot$gate==m[1]]),log(subdataplot$GFP.A[subdataplot$gate==m[2]]))
abline(v=c(testup$estimate[1],testup$estimate[2],testdown$estimate[2]))
#legend("topleft",legend = c(paste(typetest[1],"up",round(testup$estimate[2]-testup$estimate[1],digits = 3),round(-log10(testup$p.value),0),sep=" "),
#                            paste(typetest[1],"down",round(testdown$estimate[2]-testdown$estimate[1],digits = 3),round(-log10(testdown$p.value),0),sep=" ")),
#       lty=1,col=c(rgb(0.9,0.3,0.1,1),rgb(0,0,1,1)))
dev.off()

d1 <- density(log(plotdata$mCherry.A[plotdata$gate==m[1]]),na.rm = T) 
d2 <- density(log(plotdata$mCherry.A[plotdata$gate==m[2]]),na.rm = T) # returns the density data
d3 <- density(log(plotdata$mCherry.A[plotdata$gate==m[3]]),na.rm = T) # returns the density data
svg(filename = "density-GPF1-mCH-mchpop.svg",width = 5,height = 3.5)
plot(d1, lwd = 2, main="", col = rgb(0,0,0,1), xlim=c(5,9),xlab="mCherry",ylim=c(0,1.1*max(c(d1$y,d2$y,d3$y)))) # plots the results 
points(d2,typ = "l", lwd = 2, col = rgb(0,0,1,1))
points(d3,typ = "l", lwd = 2, col = rgb(0.9,0.3,0.1,1))
testup<-t.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[3]]))
testdown<-t.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[2]]))
abline(v=c(testup$estimate[1],testup$estimate[2],testdown$estimate[2]))
dev.off()

m<-mat[2,]
d1 <- density(log(subdataplot$GFP.A[subdataplot$gate==m[1]]),na.rm = T) 
d2 <- density(log(subdataplot$GFP.A[subdataplot$gate==m[2]]),na.rm = T) # returns the density data
d3 <- density(log(subdataplot$GFP.A[subdataplot$gate==m[3]]),na.rm = T) # returns the density data
svg(filename = "density-GPF1-GFP-GFPpop.svg",width = 5,height = 3.5)
plot(d1, lwd = 2, main="", col = rgb(0,0,0,1), xlim=c(5,10),xlab="GFP",ylim=c(0,1.1*max(c(d1$y,d2$y,d3$y)))) # plots the results 
points(d2,typ = "l", lwd = 2, col = rgb(0,0,1,1))
points(d3,typ = "l", lwd = 2, col = rgb(0.9,0.3,0.1,1))
testup<-t.test(log(subdataplot$GFP.A[subdataplot$gate==m[1]]),log(subdataplot$GFP.A[subdataplot$gate==m[3]]))
testdown<-t.test(log(subdataplot$GFP.A[subdataplot$gate==m[1]]),log(subdataplot$GFP.A[subdataplot$gate==m[2]]))
abline(v=c(testup$estimate[1],testup$estimate[2],testdown$estimate[2]))
dev.off()
#legend("topleft",legend = c(paste(typetest[1],"up",round(testup$estimate[2]-testup$estimate[1],digits = 3),round(-log10(testup$p.value),0),sep=" "),
#                            paste(typetest[1],"down",round(testdown$estimate[2]-testdown$estimate[1],digits = 3),round(-log10(testdown$p.value),0),sep=" ")),
#       lty=1,col=c(rgb(0.9,0.3,0.1,1),rgb(0,0,1,1)))

d1 <- density(log(plotdata$mCherry.A[plotdata$gate==m[1]]),na.rm = T) 
d2 <- density(log(plotdata$mCherry.A[plotdata$gate==m[2]]),na.rm = T) # returns the density data
d3 <- density(log(plotdata$mCherry.A[plotdata$gate==m[3]]),na.rm = T) # returns the density data
svg(filename = "density-GPF1-mCH-GFPpop.svg",width = 5,height = 3.5)
plot(d1, lwd = 2, main="", col = rgb(0,0,0,1), xlim=c(5,9),xlab="mCherry",ylim=c(0,1.1*max(c(d1$y,d2$y,d3$y)))) # plots the results 
points(d2,typ = "l", lwd = 2, col = rgb(0,0,1,1))
points(d3,typ = "l", lwd = 2, col = rgb(0.9,0.3,0.1,1))
testup<-t.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[3]]))
testdown<-t.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[2]]))
abline(v=c(testup$estimate[1],testup$estimate[2],testdown$estimate[2]))
dev.off()



