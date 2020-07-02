#Import data, process calibation and quantify RNA for RNA seq data (from 2017/02/20)

setwd("YOUR/WORKING/DIRECTORY")

#import data
qPCRdata<-read.table("qPCRdata.txt",header = T,sep = "\t", quote = "", na.strings = "NaN",stringsAsFactors = F)
qPCRdataDNA<-qPCRdata[qPCRdata$type=="DNA",]

col=c("#34ab83ff","#6a6a6aff","red","#6a6a6aff")

#fit calibation curves for the 4 paires of primers
#ax+b
a<-c()
b<-c()
plot(0,0,xlim=c(0,2),ylim=c(15,30),typ="n",xlab="DNA cons (log2)",ylab="Cq")
for (i in 1:4) {
  points(log2(qPCRdataDNA$DNAcons[qPCRdataDNA$Primer==i]),qPCRdataDNA$Cq[qPCRdataDNA$Primer==i],col=col[i])
  reg<-lm(qPCRdataDNA$Cq[qPCRdataDNA$Primer==i]~log2(qPCRdataDNA$DNAcons[qPCRdataDNA$Primer==i]))
  abline(reg,lwd=2,col=col[i])
  a[i]<-reg$coefficients[2]
  b[i]<-reg$coefficients[1]
}
primerCoef<-data.frame(a,b)
primerCoef$eff<-(2^-(1/primerCoef$a)-1)*100

#aplly calibation coef on qPCR data
qPCRdataExp<-qPCRdata[qPCRdata$type=="RT",]
qPCRdataExp$eqConslog2<-rep(0,nrow(qPCRdataExp))
qPCRdataExp<-qPCRdataExp[order(qPCRdataExp$Primer),]
qPCRdataExp<-qPCRdataExp[order(qPCRdataExp$RNAsample),]
for (i in 1:4) {
  qPCRdataExp$eqConslog2[qPCRdataExp$Primer==i]<-(qPCRdataExp$Cq[qPCRdataExp$Primer==i]-primerCoef$b[i])/primerCoef$a[i]
}
qPCRdataExp$eqCons<-2^qPCRdataExp$eqConslog2
qPCRdataExp$eqCons[is.na(qPCRdataExp$eqCons)]<-0

col2<-rep(col,2)[order(rep(1:4,2))]
barplot(qPCRdataExp$eqCons,col=col2)#,ylim=c(0,25)

barplot(qPCRdataExp$eqCons[qPCRdataExp$RNAsample%in%c(7)],ylim=c(0,0.4),col=col2)

qPCRdataExp2<-qPCRdataExp[qPCRdataExp$RNAsample%in%c(3,7),][order(qPCRdataExp$Primer[qPCRdataExp$RNAsample%in%c(3,7)]),]

col3<-rep(col,4)[order(rep(1:4,4))]

#export barplot
svg(filename = "qpcrCUT.svg",width = 4.3,height = 4)
barplot(qPCRdataExp2$eqCons,ylim=c(0,0.4),col=col3)
dev.off()
