#importing qPCR data (2019-07-31) and fluorescence, re-runing calibation data (with ACT1 primer this time) and processing data compared to ACT1 mRNA quantification

setwd("YOUR/WORKING/DIRECTORY")

#importing data
data<-read.table(paste("allDATA.txt",sep=""),header = T,sep = "\t",na.strings = c("NaN","NA")) #estradiol and YAK1 data (with ACT1 calibation)
qPCRdata<-read.table("qPCRdata.txt",header = T,sep = "\t", quote = "", na.strings = "NaN",stringsAsFactors = F) #other primer calibation data (frome 2017-02-20)
qPCRdataDNA<-qPCRdata[qPCRdata$type=="DNA",]


#run primer calibation and generate the plot
svg(filename = "qPCReff.svg",width = 6.5,height = 6)
col=c("#34ab83ff","#6a6a6aff","red","#6a6a6aff")
#ax+b
a<-c()
b<-c()
plot(0,0,xlim=c(-2,4),ylim=c(15,30),typ="n",xlab="DNA cons (log2)",ylab="Cq")
for (i in 1:4) {
  points(log2(qPCRdataDNA$DNAcons[qPCRdataDNA$Primer==i]),qPCRdataDNA$Cq[qPCRdataDNA$Primer==i],col=col[i])
  reg<-lm(qPCRdataDNA$Cq[qPCRdataDNA$Primer==i]~log2(qPCRdataDNA$DNAcons[qPCRdataDNA$Primer==i]))
  abline(reg,lwd=2,col=col[i])
  a[i]<-reg$coefficients[2]
  b[i]<-reg$coefficients[1]
}
dataDNA<-data[data$test=="DNA" & data$cond>0.4,]
points(log2(dataDNA$cond),dataDNA$ACT1,col="blue")
reg<-lm(dataDNA$ACT1~log2(dataDNA$cond))
abline(reg,lwd=2,col="blue")
a[i+1]<-reg$coefficients[2]
b[i+1]<-reg$coefficients[1]
primerCoef<-data.frame(a,b)
primerCoef$eff<-(2^-(1/primerCoef$a)-1)*100
dev.off()

#qPCR correction, multiple correction tested
maxGFPr<-data$GFPr[data$Row=="A" & data$Col==9]
maxgRNA<-data$gRNA[data$Row=="A" & data$Col==9]
maxACT1<-data$ACT1[data$Row=="A" & data$Col==9]
# data$ACT1[is.na(data$ACT1)]<-maxACT1
# data$gRNA[is.na(data$gRNA)]<-maxgRNA
# data$GFPr[is.na(data$GFPr)]<-maxGFPr

plot(data$cons,-data$ACT1)
plot(log(data$cons),-data$ACT1)
plot(data$ACT1,data$GFPr)
regM<-lm(data$GFPr ~ data$ACT1)
plot(data$ACT1,data$gRNA)
regG<-lm(data$gRNA ~ data$ACT1)

data$ACT1cons<-2^((data$ACT1-b[5])/a[5]) #keep this one but then use ACT1 as reference for log2-ratio
data$mRcons<-2^((data$GFPr-b[1])/a[1])
data$gRcons<-2^((data$gRNA-b[3])/a[3])

plot(data$cons,data$ACT1cons,ylim=c(0,1))
plot(data$cons,data$mRcons,ylim=c(0,1))
plot(data$cons,data$gRcons,ylim=c(0,1))

avMax<-mean(c(maxGFPr,maxgRNA,maxACT1))
data$iACT1 <- -data$ACT1+avMax
averageACT1 <- mean(data$iACT1,na.rm=T)
data$cGFPr <- -data$GFPr+avMax -data$iACT1 +averageACT1
data$cgRNA <- -data$gRNA+avMax -data$iACT1 +averageACT1


plot(data$cGFPr,data$iACT1)
plot(data$cgRNA,data$iACT1)
plot(data$cGFPr~data$test)
plot(data$cgRNA~data$test)

#florescence plot and correction
ODblanc<-mean(data$OD[data$test=="Blanc"])
GFPblanc<-mean(data$GFP[data$test=="Blanc"])
mCHblanc<-mean(c(data$mCH1[data$test=="Blanc"],data$mCH2[data$test=="Blanc"],data$mCH2[data$test=="Blanc"]))
data$ODc<-data$OD-ODblanc
data$GFPc<-(data$GFP-GFPblanc)/data$ODc
data$mCH1c<-(data$mCH1-mCHblanc)/data$ODc
data$mCH2c<-(data$mCH2-mCHblanc)/data$ODc
data$mCH3c<-(data$mCH3-mCHblanc)/data$ODc
data$mCHcM<-rowMeans(data[,colnames(data)%in%c("mCH1","mCH2","mCH3")],na.rm=T)


#Estradiol qPCR curves (figure2)
dataE<-data[data$test=="estra",]
dataE<-dataE[dataE$cons>5,]

hist(dataE$cons,breaks = 30)

plot(dataE$ACT1cons,dataE$mRcons,typ="n",main="qPCR")
text(dataE$ACT1cons,dataE$mRcons,labels = dataE$cond)

plot(dataE$ACT1cons,dataE$gRcons,typ="n",main="qPCR")
text(dataE$ACT1cons,dataE$gRcons,labels = dataE$cond)

plot(dataE$ACT1cons,dataE$gRcons/dataE$mRcons,typ="n",main="qPCR")
text(dataE$ACT1cons,dataE$gRcons/dataE$mRcons,labels = dataE$cond)

pdf(file = "estraQ_GFPvmCH.pdf",width = 6.5,height = 6)
plot(dataE$GFPc,dataE$mCHcM,typ="n",main="flurescenses",ylim=c(400,1000),xlim=c(0,25000))
text(dataE$GFPc,dataE$mCHcM,labels = dataE$cond)
dev.off()

dataE$lGFP<-dataE$GFPc
dataE$lGFP[dataE$lGFP<1]<-1
dataE$lGFP<-log2(dataE$lGFP)
dataE$lmch<-dataE$mCHcM
dataE$lmch[dataE$lmch<1]<-1
dataE$lmch<-log2(dataE$lmch)

plot(dataE$lGFP,dataE$lmch,typ="n",main="flurescenses")#,ylim=c(400,1000),xlim=c(0,25000))
text(dataE$lGFP,dataE$lmch,labels = dataE$cond)

pdf(file = "estraQ_GFPvmRNA.pdf",width = 6.5,height = 6)
plot(log2(dataE$mRcons/dataE$ACT1cons),dataE$GFPc,typ="n",main="flurescenses",ylim=c(0,25000))
text(log2(dataE$mRcons/dataE$ACT1cons),dataE$GFPc,labels = dataE$cond,col="#34ab83ff")
dev.off()

plot(log2(dataE$mRcons/dataE$ACT1cons),dataE$lGFP,typ="n",main="flurescenses")#,ylim=c(0,25000))
text(log2(dataE$mRcons/dataE$ACT1cons),dataE$lGFP,labels = dataE$cond,col="#34ab83ff")

pdf(file = "estraQ_mCHvmRNA.pdf",width = 6.5,height = 6)
plot(log2(dataE$mRcons/dataE$ACT1cons),dataE$mCHcM,typ="n",main="flurescenses",ylim=c(400,1000))
text(log2(dataE$mRcons/dataE$ACT1cons),dataE$mCHcM,labels = dataE$cond,col="red")
reg<-lm(dataE$mCHcM[dataE$cond<4]~log2(dataE$mRcons/dataE$ACT1cons)[dataE$cond<4])
abline(reg)
dev.off()

plot(log2(dataE$mRcons/dataE$ACT1cons),dataE$lmch,typ="n",main="flurescenses")#,ylim=c(0,25000))
text(log2(dataE$mRcons/dataE$ACT1cons),dataE$lmch,labels = dataE$cond,col="red")
reg<-lm(dataE$lmch[dataE$cond<4]~log2(dataE$mRcons/dataE$ACT1cons)[dataE$cond<4])
abline(reg)

pdf(file = "estraQ_mCHvgRNA.pdf",width = 6.5,height = 6)
plot(log2(dataE$gRcons/dataE$ACT1cons),dataE$mCHcM,typ="n",main="flurescenses",ylim=c(400,1000))
text(log2(dataE$gRcons/dataE$ACT1cons),dataE$mCHcM,labels = dataE$cond,col="red")
reg<-lm(dataE$mCHcM[dataE$cond<4]~log2(dataE$gRcons/dataE$ACT1cons)[dataE$cond<4])
abline(reg)
dev.off()

plot(log2(dataE$gRcons/dataE$ACT1cons),dataE$lmch,typ="n",main="flurescenses")#,ylim=c(400,1000))
text(log2(dataE$gRcons/dataE$ACT1cons),dataE$lmch,labels = dataE$cond,col="red")
reg<-lm(dataE$lmch[dataE$cond<4]~log2(dataE$gRcons/dataE$ACT1cons)[dataE$cond<4])
abline(reg)


pdf(file = "estraQ_gRNAvmRNA.pdf",width = 6.5,height = 6)
plot(log2(dataE$mRcons/dataE$ACT1cons),log2(dataE$gRcons/dataE$ACT1cons),typ="n",main="qpcr")
text(log2(dataE$mRcons/dataE$ACT1cons),log2(dataE$gRcons/dataE$ACT1cons),labels = dataE$cond)
reg<-lm(log2(dataE$gRcons/dataE$ACT1cons)[dataE$cond<4]~log2(dataE$mRcons/dataE$ACT1cons)[dataE$cond<4])
abline(reg)
dev.off()

plot(log2(dataE$mRcons/dataE$ACT1cons),log2(dataE$gRcons/dataE$mRcons),typ="n",main="qpcr")
text(log2(dataE$mRcons/dataE$ACT1cons),log2(dataE$gRcons/dataE$mRcons),labels = dataE$cond)


svg(filename = "estraQ_mRNAgRNARatio.svg",width = 6.5,height = 6)
vRGBini=c(124,231,229)
vRGB=t(vRGBini%*%(t((6:0)+6)/12))/255
colestra=rgb(vRGB)
levest<-factor(dataE$cond,levels = c("0","0.5","1","2","4","8","16"))
boxplot(log2(dataE$gRcons/dataE$mRcons)~levest,main="qPCR",ylim=c(-3,-1),col=colestra[1:7])
points(as.numeric(levest)+runif(n = length(levest),min = -0.12,max = 0.12),log2(dataE$gRcons/dataE$mRcons),pch=16,cex=1)
dev.off()

#exploring YAK1 data in a PDF
pdf(file = "190731-ChB-ESTRA-YAK1-qPCR.pdf",width = 8,height = 6.5)

dataE<-data[data$test=="estra",]

plot(dataE$GFPc,dataE$mCHcM,typ="n",main="flurescenses")
text(dataE$GFPc,dataE$mCHcM,labels = dataE$cond)


plot(dataE$cGFPr,dataE$GFPc,typ="n",main="qPCR vs fluorescence: GFP")#cex=dataE$cond/3+0.5
text(dataE$cGFPr,dataE$GFPc,labels = dataE$cond)
plot(dataE$cgRNA,dataE$mCHcM,typ="n",main="qPCR vs fluorescence: gRNA/mCherry")#cex=dataE$cond/3+0.5
text(dataE$cgRNA,dataE$mCHcM,labels = dataE$cond)
plot(dataE$cGFPr,dataE$mCHcM,typ="n",main="qPCR(GFPseq) vs fluorescence(mCherry)")#cex=dataE$cond/3+0.5
text(dataE$cGFPr,dataE$mCHcM,labels = dataE$cond)

plot(dataE$cGFPr,dataE$cgRNA,typ="n",main="qPCRs")#cex=dataE$cond/3+0.5
text(dataE$cGFPr,dataE$cgRNA,labels = dataE$cond)

levels(data$test)

dataY<-data[data$test %in% c("GPD1GFP","GPD1GFPgRNA"),]
dataY$cat<-as.factor(paste(dataY$test,dataY$YAK1,sep="_"))

plot(dataY$cGFPr~dataY$cat,col="green",main="qPCR")
text(as.numeric(dataY$cat)+as.numeric(as.character(dataY$Clone))*0.07-0.4,dataY$cGFPr,labels = dataY$Clone,cex=0.7, col='black')

plot(dataY$cgRNA~dataY$cat,col="red",main="qPCR")
text(as.numeric(dataY$cat)+as.numeric(as.character(dataY$Clone))*0.07-0.4,dataY$cgRNA,labels = dataY$Clone,cex=0.7, col='black')

plot(log2(dataY$GFPc)~dataY$cat,col="green",main="flurescense (log)")
text(as.numeric(dataY$cat)+as.numeric(as.character(dataY$Clone))*0.07-0.4,log2(dataY$GFPc),labels = dataY$Clone,cex=0.7, col='black')

plot(log2(dataY$mCHcM)~dataY$cat,col="red",main="flurescense (log)")
text(as.numeric(dataY$cat)+as.numeric(as.character(dataY$Clone))*0.07-0.4,log2(dataY$mCHcM),labels = dataY$Clone,cex=0.7, col='black')

dataY$logGFP<-log2(dataY$GFPc)

#pvalue
dataYv<-data[data$test %in% c("GPD1GFP"),]
dataYv$cat<-as.factor(paste(dataYv$test,dataYv$YAK1,sep="_"))

t.test(dataYv$cGFPr~dataYv$cat)
t.test(log(dataYv$GFPc)~dataYv$cat)
dev.off()


#ploting in SVG format
dataY<-data[data$test %in% c("GPD1GFP","GPD1GFPgRNA"),]
dataY$cat<-as.factor(paste(dataY$test,dataY$YAK1,sep="_"))

widplotsvg<-3.2
#GFP yak1 variant
dataYtemp<-dataY[dataY$test=="GPD1GFP",]
dataYtemp$YAK1<-factor(as.character(dataYtemp$YAK1),levels=c("WT","VAR"))
dataYtemp$colg<-rep("blue",nrow(dataYtemp))
dataYtemp$colg[dataYtemp$YAK1=="VAR"]<-"#5c5cc1ff"

svg(filename = "YAK1_var_qPCRgfp.svg",width = widplotsvg,height = 4.5)
plot(log2(dataYtemp$mRcons/dataYtemp$ACT1cons)~dataYtemp$YAK1,main="qPCR",ylim=c(-3,2),col="#34ab83ff",cex=0.5)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.175,log2(dataYtemp$mRcons/dataYtemp$ACT1cons),cex=1,pch=16, col=dataYtemp$colg)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.175,log2(dataYtemp$mRcons/dataYtemp$ACT1cons),cex=1,pch=1, col="black")
test<-t.test(log2(dataYtemp$mRcons/dataYtemp$ACT1cons)~dataYtemp$YAK1)
print(paste(test$estimate[2]-test$estimate[1],test$p.value))
dev.off()

svg(filename = "YAK1_var_FLUOgfp.svg",width = widplotsvg,height = 4.5)
plot(log2(dataYtemp$GFPc)~dataYtemp$YAK1,main="fluo",ylim=c(11.6,12.2),col="#34ab83ff",cex=0.5)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.175,log2(dataYtemp$GFPc),cex=1,pch=16, col=dataYtemp$colg)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.175,log2(dataYtemp$GFPc),cex=1,pch=1, col="black")
test<-t.test(log2(dataYtemp$GFPc)~dataYtemp$YAK1)
print(paste(test$estimate[2]-test$estimate[1],test$p.value))
dev.off()

#GFP mCH yak1 deletion
dataYtemp<-dataY[dataY$test=="GPD1GFPgRNA",]
dataYtemp$YAK1<-factor(as.character(dataYtemp$YAK1),levels=c("WT","delta"))
dataYtemp$colg<-rep("blue",nrow(dataYtemp)) #"#34ab83ff"
dataYtemp$colg[dataYtemp$YAK1=="delta"]<-"white"
dataYtemp$colm<-rep("blue",nrow(dataYtemp)) #"red"
dataYtemp$colm[dataYtemp$YAK1=="delta"]<-"white"

svg(filename = "YAK1_dlta_qPCRgfp.svg",width = widplotsvg,height = 4.5)
plot(log2(dataYtemp$mRcons/dataYtemp$ACT1cons)~dataYtemp$YAK1,main="qPCR",ylim=c(-3,2),col="#34ab83ff",cex=0.5)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.35,log2(dataYtemp$mRcons/dataYtemp$ACT1cons),cex=1,pch=16, col=dataYtemp$colg)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.35,log2(dataYtemp$mRcons/dataYtemp$ACT1cons),cex=1,pch=1, col="black")
test<-t.test(log2(dataYtemp$mRcons/dataYtemp$ACT1cons)~dataYtemp$YAK1)
print(paste(test$estimate[2]-test$estimate[1],test$p.value))
dev.off()

svg(filename = "YAK1_dlta_FLUOgfp.svg",width = widplotsvg,height = 4.5)
plot(log2(dataYtemp$GFPc)~dataYtemp$YAK1,main="fluo",ylim=c(12.4,13.4),col="#34ab83ff",cex=0.5)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.35,log2(dataYtemp$GFPc),cex=1,pch=16, col=dataYtemp$colg)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.35,log2(dataYtemp$GFPc),cex=1,pch=1, col="black")
test<-t.test(log2(dataYtemp$GFPc)~dataYtemp$YAK1)
print(paste(test$estimate[2]-test$estimate[1],test$p.value))
dev.off()

svg(filename = "YAK1_dlta_qPCRmch.svg",width = widplotsvg,height = 4.5)
plot(log2(dataYtemp$gRcons/dataYtemp$ACT1cons)~dataYtemp$YAK1,main="qPCR",ylim=c(-5,-1),col="red",cex=0.5)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.35,log2(dataYtemp$gRcons/dataYtemp$ACT1cons),cex=1,pch=16, col=dataYtemp$colm)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.35,log2(dataYtemp$gRcons/dataYtemp$ACT1cons),cex=1,pch=1, col="black")
test<-t.test(log2(dataYtemp$gRcons/dataYtemp$ACT1cons)~dataYtemp$YAK1)
print(paste(test$estimate[2]-test$estimate[1],test$p.value))
dev.off()

svg(filename = "YAK1_dlta_FLUOmch.svg",width = widplotsvg,height = 4.5)
plot(log2(dataYtemp$mCHcM)~dataYtemp$YAK1,main="fluo",ylim=c(9.6,10.4),col="red",cex=0.5)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.35,log2(dataYtemp$mCHcM),cex=1,pch=16, col=dataYtemp$colm)
points(as.numeric(dataYtemp$YAK1)+as.numeric(as.character(dataYtemp$Clone))*0.07-0.35,log2(dataYtemp$mCHcM),cex=1,pch=1, col="black")
test<-t.test(log2(dataYtemp$mCHcM)~dataYtemp$YAK1)
print(paste(test$estimate[2]-test$estimate[1],test$p.value))
dev.off()

