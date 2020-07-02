library(flowCore)
library(reshape2)
library(plyr)
library(data.table)
par(mfrow = c(1, 1))

setwd("YOUR/WORKING/DIRECTORY")
file<-all[1]

import_dataO <- function(file){
  x = unlist(strsplit(basename(file), "\\."))[1]
  name = unlist(strsplit(basename(x), "\\_"))[2]
  gate = unlist(strsplit(basename(x), "\\_"))[3]
  #print(paste("Reading", condition, id, time, sep=" "))  
  data <- read.FCS(file, transformation=FALSE)
  print(paste("Reading", x, nrow(data), sep=" "))  
  mdata <- data.frame(exprs(data))
  mdata$gene <- name
  mdata$gate <- as.numeric(gate)
  return(mdata)
}
import_data1 <- function(file){
  x = unlist(strsplit(basename(file), "\\."))[1]
  name = unlist(strsplit(basename(x), "\\_"))[3]
  gate = unlist(strsplit(basename(x), "\\_"))[4]
  #print(paste("Reading", condition, id, time, sep=" "))  
  data <- read.FCS(file, transformation=FALSE)
  print(paste("Reading", x, nrow(data), sep=" "))  
  mdata <- data.frame(exprs(data))
  mdata$gene <- name
  mdata$gate <- as.numeric(gate)
  return(mdata)
}

all = list.files(path="180419-CHB-Sort32", pattern = ".fcs", full = TRUE)
data_all = ldply(all, import_dataO) #plyr library
all2 = list.files(path="Sort_3-1-180328", pattern = ".fcs", full = TRUE)
data_all2 = ldply(all2, import_data1) #plyr library
colnames(data_all)[colnames(data_all)=="eGFP.A"]<-"GFP.A"
colnames(data_all2)[colnames(data_all2)=="eGFP.A"]<-"GFP.A"


#only single
plot(log(data_all$FSC.A),log(data_all$SSC.A),pch=16,col=rgb(0,0,1,0.01),cex=0.4,ylim=c(9,12.5),xlim=c(10,12.5))
gateTOPx<-c(10,10.5,11.6,11.9)
gateTOPy<-c(9.7,10.3,11.4,11.2)
gateBOTx<-c(10,10.4,11.3,11.9)
gateBOTy<-c(9.7,9.6,10.4,11.2)
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
#plot(log(subdata$FSC.A),log(subdata$SSC.A),pch=16,col=rgb(0,0,1,0.20),cex=0.4,ylim=c(8,12))
#plot(subdata$FSC.A,subdata$SSC.A,pch=16,col=rgb(0,0,1,0.01),cex=0.4)#,ylim=c(0,260),xlim=c(0,260))

subdata<-subdata[subdata$FSC.A>0 & subdata$SSC.A>0 & subdata$mCherry.A>0 & subdata$GFP.A>0,]



#only single
plot(log(data_all2$FSC.A),log(data_all2$SSC.A),pch=16,col=rgb(0,0,1,0.01),cex=0.4,ylim=c(9,12.5),xlim=c(10,12.5))
gateTOPx<-c(10,10.5,11.6,11.9)
gateTOPy<-c(9.7,10.3,11.4,11.2)
gateBOTx<-c(10,10.4,11.3,11.9)
gateBOTy<-c(9.7,9.6,10.4,11.2)
points(c(gateTOPx,rev(gateBOTx)),c(gateTOPy,rev(gateBOTy)),typ="l",col="grey",lwd=2)


subdata2<-data_all2
for (i in 1:(length(gateTOPx)-1)) {
  coor<-lm(gateTOPy[i:(i+1)]~gateTOPx[i:(i+1)])
  subdata2<-subdata2[log(subdata2$SSC.A) < coor$coefficients[2]*log(subdata2$FSC.A)+coor$coefficients[1],]
}
for (i in 1:(length(gateBOTx)-1)) {
  coor<-lm(gateBOTy[i:(i+1)]~gateBOTx[i:(i+1)])
  subdata2<-subdata2[log(subdata2$SSC.A) > coor$coefficients[2]*log(subdata2$FSC.A)+coor$coefficients[1],]
}
plot(log(subdata2$FSC.A),log(subdata2$SSC.A),pch=16,col=rgb(0,0,1,0.20),cex=0.4,ylim=c(8,12))
plot(subdata2$FSC.A,subdata2$SSC.A,pch=16,col=rgb(0,0,1,0.01),cex=0.4)#,ylim=c(0,260),xlim=c(0,260))

subdata2<-subdata2[subdata2$FSC.A>0 & subdata2$SSC.A>0 & subdata2$mCherry.A>0 & subdata2$GFP.A>0,]

subdata<-rbind(subdata,subdata2)



strains<-levels(as.factor(subdata$gene))
typetest<-c("whole pop mCH","whole pop GFP","fixed mCH","fixed GFP","fixeeffect")
mat<-matrix(c(1,2,3,1,4,5,6,7,8,9,10,11,1,6,9),nrow = 5,ncol = 3,byrow = T)

lastlist <- function(x) {x[length(x)]}
subdata$CRISPR<-unlist(lapply(strsplit(subdata$gene,""), lastlist))

genename<-function(x) {
  paste(strsplit(x,"")[[1]][-length(strsplit(x,"")[[1]])],sep="",collapse = "")
}
subdata$GOI<-unlist(lapply(subdata$gene,genename))

subdata$mCHc<-subdata$mCherry.A
subdata$GFPc<-subdata$mCherry.A

subdata<-subdata[!is.na(subdata$gate),]

subdata_temp<-subdata
subdata<-subdata_temp
subdata<-subdata[subdata$gate==1,]

gene <- strains[1]
plot(log(subdata$mCherry.A[subdata$gene==gene])~subdata$FSC.A[subdata$gene==gene])
test<-loess(log(subdata$mCherry.A[subdata$gene==gene])~subdata$FSC.A[subdata$gene==gene])
plot(unlist(test$residuals)~subdata$FSC.A[subdata$gene==gene])

plot(log(subdata$GFP.A[subdata$gene==gene])~subdata$FSC.A[subdata$gene==gene])
test2<-loess(log(subdata$GFP.A[subdata$gene==gene])~subdata$FSC.A[subdata$gene==gene])
plot(unlist(test2$residuals)~subdata$FSC.A[subdata$gene==gene])


correl<-c()
r2<-c()
correlc<-c()
r2c<-c()
strains<-levels(as.factor(subdata$gene))
strains <- strains[c(5,6,10,11,3,4,7,8,12,13,1,2,9)]

pdf("180504_Flow_effectSize.pdf", width = 12, height = 12)
par(mfcol = c(2, 2))
subdata<-subdata_temp
#subdata<-subdata[subdata$gate%in%c(1,2,3,4,5),]
for (i in 1:length(strains)) {
  print(paste(i,strains[i],date()))
  gene<-strains[][i]
  test<-loess(log(subdata$mCherry.A[subdata$gene==gene])~subdata$FSC.A[subdata$gene==gene])
  subdata$mCHc[subdata$gene==gene]<-unlist(test$residuals)+mean(log(subdata$mCherry.A[subdata$gene==gene]))
  test2<-loess(log(subdata$GFP.A[subdata$gene==gene])~subdata$FSC.A[subdata$gene==gene])
  subdata$GFPc[subdata$gene==gene]<-unlist(test2$residuals)+mean(log(subdata$GFP.A[subdata$gene==gene]))
  reg<-lm(log(subdata$mCherry.A[subdata$gene==gene])~log(subdata$GFP.A[subdata$gene==gene]))
  reg2<-lm(unlist(test$residuals)~unlist(test2$residuals))
  
  plot(log(subdata$mCherry.A[subdata$gene==gene])~log(subdata$GFP.A[subdata$gene==gene]),pch=16,cex=0.6,col=rgb(0,0,0,0.12),xlim=c(2,12),ylim=c(4,10),main=paste(gene),xlab = "GFP",ylab = "mCherry")
  reg<-lm(log(subdata$mCherry.A[subdata$gene==gene])~log(subdata$GFP.A[subdata$gene==gene]))
  summary(reg)
  abline(reg,col="red")
  text(x = 3, y = 9, paste("s:",round(reg$coefficients[[2]],3)))
  text(x = 3, y = 8, paste("r:",round(summary(reg)$r.squared,3)))
  
  plot(subdata$mCHc[subdata$gene==gene]~subdata$GFPc[subdata$gene==gene],pch=16,cex=0.6,col=rgb(0,0,0,0.12),xlim=c(2,12),ylim=c(4,10),main=paste(gene,"corrected"),xlab = "GFP",ylab = "mCherry")
  reg2<-lm(subdata$mCHc[subdata$gene==gene]~subdata$GFPc[subdata$gene==gene])
  summary(reg2)
  abline(reg2,col="red")
  text(x = 3, y = 9, paste("s:",round(reg2$coefficients[[2]],3)))
  text(x = 3, y = 8, paste("r2:",round(summary(reg2)$r.squared,3)))
  
  
  correl[i]<-reg$coefficients[[2]]
  correlc[i]<-reg2$coefficients[[2]]
  r2[i]<-summary(reg)$r.squared
  r2c[i]<-summary(reg2)$r.squared
}
print(date()) # 12min

gene <- "NAi"
plot(log(subdata$mCherry.A[subdata$gene==gene])~subdata$FSC.A[subdata$gene==gene],pch=16,cex=0.6,col=rgb(0,0,0,0.12),ylim=c(4,10),xlab="FSC",ylab="mCherry",main="correction example (NAi)")
test<-loess(log(subdata$mCherry.A[subdata$gene==gene])~subdata$FSC.A[subdata$gene==gene])
plot(unlist(test$residuals)+mean(log(subdata$mCherry.A[subdata$gene==gene]))~subdata$FSC.A[subdata$gene==gene],pch=16,cex=0.6,col=rgb(0,0,0,0.12),ylim=c(4,10),xlab="FSC",ylab="corrected mCherry (loess)")
dev.off()

sizecor<-data.frame(strains,correl,r2,correlc,r2c)
write.table(sizecor,"FLOWsizecor_12.txt",sep="\t",row.names = FALSE, quote = F)

#subdata<-subdata_temp


lastlist <- function(x) {x[length(x)]}
GOI<-rep(strains,8)[order(rep(strains,8))]
CRISPR<-unlist(lapply(strsplit(basename(GOI),""), lastlist))
final<-data.frame(GOI,CRISPR,stringsAsFactors = F)
final$chanel<-rep(c("mCH","mCH","GFP","GFP"),length(GOI)/4)
final$pop<-rep(c("high","low","high","low"),length(GOI)/4)
final$atlchanel<-rep(c(rep("all",4),rep("fix",4)),length(GOI)/8)
final$mCHeffect<-rep(0,length(GOI))
final$mCHpval<-rep(0,length(GOI))
final$GFPeffect<-rep(0,length(GOI))
final$GFPpval<-rep(0,length(GOI))
final$finalPV<-rep(0,length(GOI))
final$heritability<-rep("",length(GOI))
final$synergy<-rep("",length(GOI))
finalc<-final
finalcD<-finalc[finalc$pop=="high",]

pdf("180417_Flow_sortAll.pdf", width = 22, height = 12)
par(mfcol = c(3, 4))

n=1
for (gene in strains) {
  plotdata<-subdata[subdata$gene==gene,]
  for (i in 1:4) {
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
    if (i != 5) {
      testup<-t.test(log(plotdata$GFP.A[plotdata$gate==m[1]]),log(plotdata$GFP.A[plotdata$gate==m[3]]))
      testdown<-t.test(log(plotdata$GFP.A[plotdata$gate==m[1]]),log(plotdata$GFP.A[plotdata$gate==m[2]]))
      legend("topleft",legend = c(paste(typetest[i],"up",round(testup$estimate[2]-testup$estimate[1],digits = 3),round(-log10(testup$p.value),0),sep=" "),
                                  paste(typetest[i],"down",round(testdown$estimate[2]-testdown$estimate[1],digits = 3),round(-log10(testdown$p.value),0),sep=" ")),
             lty=1,col=c(rgb(0.9,0.3,0.1,1),rgb(0,0,1,1)))
      final$GFPeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$estimate[2]-testup$estimate[1]
      final$GFPpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$p.value
      final$GFPeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$estimate[2]-testdown$estimate[1]
      final$GFPpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$p.value
    } else if (i == 5) {
      testup<-var.test(log(plotdata$GFP.A[plotdata$gate==m[1]]),log(plotdata$GFP.A[plotdata$gate==m[3]]))
      testdown<-var.test(log(plotdata$GFP.A[plotdata$gate==m[1]]),log(plotdata$GFP.A[plotdata$gate==m[2]]))
      legend("topleft",legend = c(paste(typetest[i],"mCH",round(testup$estimate,digits = 3),round(-log10(testup$p.value),0),sep=" "),
                                  paste(typetest[i],"GFP",round(testdown$estimate,digits = 3),round(-log10(testdown$p.value),0),sep=" ")),
             lty=1,col=c(rgb(0.9,0.3,0.1,1),rgb(0,0,1,1)))
    }
    
    d1 <- density(log(plotdata$mCherry.A[plotdata$gate==m[1]]),na.rm = T) 
    d2 <- density(log(plotdata$mCherry.A[plotdata$gate==m[2]]),na.rm = T) # returns the density data
    d3 <- density(log(plotdata$mCherry.A[plotdata$gate==m[3]]),na.rm = T) # returns the density data
    plot(d1, lwd = 2, main="", col = rgb(0,0,0,1), xlim=c(2,12),xlab="mCherry",ylim=c(0,1.1*max(c(d1$y,d2$y,d3$y)))) # plots the results 
    points(d2,typ = "l", lwd = 2, col = rgb(0,0,1,1))
    points(d3,typ = "l", lwd = 2, col = rgb(0.9,0.3,0.1,1))
    if (i != 5) {
      testup<-t.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[3]]))
      testdown<-t.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[2]]))
      legend("topleft",legend = c(paste(typetest[i],"up",round(testup$estimate[2]-testup$estimate[1],digits = 3),round(-log10(testup$p.value),0),sep=" "),
                                  paste(typetest[i],"down",round(testdown$estimate[2]-testdown$estimate[1],digits = 3),round(-log10(testdown$p.value),0),sep=" ")),
             lty=1,col=c(rgb(0.9,0.3,0.1,1),rgb(0,0,1,1)))
      final$mCHeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$estimate[2]-testup$estimate[1]
      final$mCHpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="high" & final$atlchanel==final$atlchanel[i*2]]<-testup$p.value
      final$mCHeffect[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$estimate[2]-testdown$estimate[1]
      final$mCHpval[final$GOI==gene & final$chanel==final$chanel[i*2] & final$pop=="low" & final$atlchanel==final$atlchanel[i*2]]<-testdown$p.value
    } else if (i == 5) {
      testup<-var.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[3]]))
      testdown<-var.test(log(plotdata$mCherry.A[plotdata$gate==m[1]]),log(plotdata$mCherry.A[plotdata$gate==m[2]]))
      legend("topleft",legend = c(paste(typetest[i],"mCH",round(testup$estimate,digits = 3),round(-log10(testup$p.value),0),sep=" "),
                                  paste(typetest[i],"GFP",round(testdown$estimate,digits = 3),round(-log10(testdown$p.value),0),sep=" ")),
             lty=1,col=c(rgb(0.9,0.3,0.1,1),rgb(0,0,1,1)))
    }
    
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
    cha<-6
    oth<-8
  } else {
    cha<-8
    oth<-6
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


write.table(final,"finalFLOW3-12.txt",sep="\t",row.names = FALSE, quote = F)
write.table(finalc,"finalFLOW3-12_sc.txt",sep="\t",row.names = FALSE, quote = F)
write.table(finalcD,"finalFLOW3-12_scD.txt",sep="\t",row.names = FALSE, quote = F)
# final2<-read.table("finalFLOW3-1.txt",sep="\t",stringsAsFactors = F, header = T)
# finalF<-rbind(final,final2)
finalF<-final



barplot(as.matrix(dcast(finalF,heritability~chanel+finalF$CRISPR+finalF$GOI)[c(3,1,2),-1]))

coltyp<-c(rep(rgb(0,0.5,0.3,1),7),rep(rgb(1,0,0,1),7),rep(rgb(0,0.37,0.25,1),7),rep(rgb(0.6,0,0,1),7))
boxplot(-log10(finalF$finalPV)~finalF$GOI+finalF$chanel+finalF$CRISPR,ylim=c(0,300),col=coltyp)
points(c(0,29),c(3,3),typ="l",col="red")
points(c(0,29),c(2,2),typ="l",col="grey")

barplot(as.matrix(dcast(finalF,heritability~finalF$CRISPR+atlchanel+chanel)[c(3,1,2),-1]))

coltyp<-c(rep(rgb(0,0.5,0.3,1),2),rep(rgb(1,0,0,1),2),rep(rgb(0,0.37,0.25,1),2),rep(rgb(0.6,0,0,1),2))
boxplot(-log10(finalF$finalPV)~finalF$atlchanel+finalF$chanel+finalF$CRISPR,ylim=c(0,300),col=coltyp)
points(c(0.5,8.5),c(3,3),typ="l",col="red")
points(c(0.5,8.5),c(2,2),typ="l",col="grey")

coltyp<-c(rgb(0.75,0.45,0.15,1),rgb(1,1,1,1),rgb(0.2,0.8,1,1))
barplot(as.matrix(dcast(finalF,synergy~finalF$CRISPR+atlchanel+finalF$GOI)[,-1]),col=coltyp)

barplot(as.matrix(dcast(finalF,synergy~finalF$CRISPR+atlchanel)[,-1]),col=coltyp)

