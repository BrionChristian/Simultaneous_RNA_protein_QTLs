#this R script reads plate reader data and process growth curve and fluorescence ratio for tagged strains (GFP and mCherry)
#TDH3 and GPD1 tagged strains (data from 2017-03-09)

library(reshape2)
library(plyr)
par(mfrow = c(1, 1))


#functions used for reading plate reader data

timeC <- function(Strtime) { #trasnform time value into numeric (hours)
  val<-as.numeric(strsplit(Strtime,":")[[1]]) %*% c(3600,60,1)
  val
}

import_data <- function(file,matr){ #import read data and merfe them with matrix information
  x = unlist(strsplit(basename(file), "\\."))[1]
  date = unlist(strsplit(basename(x), "\\-"))[1]
  type = unlist(strsplit(basename(x), "\\-"))[2]
  data = read.table(file,header = T,sep = "\t", quote = "", na.strings = "")
  print(paste("Reading", x, nrow(data), sep=" "))  
  mdata <- melt(data, measure.vars=c(3:ncol(data)))
  colnames(mdata) <- c("time","temp","id","value")
  mdata$value<-as.numeric(as.character(mdata$value))
  mdata[,1] <- as.vector(mapply(timeC,as.character(mdata[,1])))
  mdata$date <- date
  mdata$type <- type
  mdata<-mdata[order(mdata$time),]
  mdata <- cbind(matr, mdata)
  return(mdata)
}

blancCorrectratio <- function(data) { #correct values from blanc wells
  BlancOD=mean(data$value[data$strain=="Blanc" & data$type=="OD"]) #mean
  BlancGFP=min(data$value[data$strain!="EMPTY" & data$type=="GFP"])
  BlancmCH=mean(data$value[data$strain=="Blanc" & data$type=="mCH"]) #mean
  datac <- data
  datac$value[datac$type=="OD"]<-datac$value[datac$type=="OD"]-BlancOD
  datac$value[datac$type=="GFP"]<-datac$value[datac$type=="GFP"]-BlancGFP
  datac$value[datac$type=="mCH"]<-datac$value[datac$type=="mCH"]-BlancmCH
  datac$value[datac$type=="OD" & datac$value < 0.001] <- 0.001
  datac$value[datac$type=="GFP" & datac$value < 0.01] <- 1
  datac$value[datac$type=="mCH" & datac$value < 0.01] <- 1
  
  datacf<-datac[datac$type=="OD",]
  datacf$GFP<-datac$value[datac$type=="GFP"]
  datacf$mCH<-datac$value[datac$type=="mCH"]
  datacf$rGFP<-datacf$GFP/datacf$value
  datacf$rmCH<-datacf$mCH/datacf$value
  
  datacf$uniqueID<-paste(datacf$strain,datacf$const,datacf$cond,sep="_")
  datacf$uniquePos<-paste(datacf$strain,datacf$const,datacf$cond,datacf$cond2,datacf$clone,datacf$id,sep="_")
  datacf$uniquest<-paste(datacf$strain,datacf$const,sep="_")
  
  
  return(datacf)
}

matrix_import <- function(file,name) { #import 96 wells matrix information and format them as data.frame
  matrix_t <- read.table(paste(file,name,sep=""),header = T,sep = "\t", quote = "", na.strings = "NA",stringsAsFactors = F)
  matrix <- melt(matrix_t,measure.vars = c(2:13))
  matrix$variable <- as.numeric(sub("X","",matrix$variable))
  matrix <- cbind(matrix,t(data.frame(lapply(strsplit(matrix$value,"[-,]"), function(u) transform(c(u,rep(NA,5-length(u))))))))
  matrix <- matrix[,-3]
  colnames(matrix) <- c("row","col","strain","const","cond","cond2","clone")
  matrix1<-matrix[order(matrix$row),]
}

smouth_func <- function(datacf) { #smooth curve to reduce point-to-point noise
  datacf$SmOD <- datacf$value
  datacf$SmGFP <- datacf$value
  datacf$SmmCH <- datacf$value
  for (id in levels(as.factor(datacf$uniquePos))) {
    x <- datacf$time[datacf$uniquePos == id & !(is.na(datacf$value))]/3600
    y <- log(datacf$value[datacf$uniquePos == id])[!is.na(datacf$value[datacf$uniquePos == id])]
    if (length(x) < 5) {next}
    sl6 = smooth.spline(x, y, spar=0.5) #0.7
    predict(sl6)$y
    datacf$SmOD[datacf$uniquePos == id & !(is.na(datacf$value))]<-exp(predict(sl6)$y)
    if (length(exp(predict(sl6)$y)) != length(datacf$SmOD[datacf$uniquePos == id & !(is.na(datacf$value))])) {print(id)}
    x <- datacf$time[datacf$uniquePos == id & !(is.na(datacf$GFP))]/3600
    y <- datacf$GFP[datacf$uniquePos == id][!is.na(datacf$GFP[datacf$uniquePos == id])]
    sl6 = smooth.spline(x, y, spar=0.5)
    predict(sl6)$y
    datacf$SmGFP[datacf$uniquePos == id & !(is.na(datacf$GFP))]<-predict(sl6)$y
    x <- datacf$time[datacf$uniquePos == id & !(is.na(datacf$mCH))]/3600
    y <- datacf$mCH[datacf$uniquePos == id][!is.na(datacf$mCH[datacf$uniquePos == id])]
    sl6 = smooth.spline(x, y, spar=0.5)
    predict(sl6)$y
    datacf$SmmCH[datacf$uniquePos == id & !(is.na(datacf$mCH))]<-predict(sl6)$y
  }
  datacf$SmGFP[datacf$SmGFP < 0.01] <- 1
  datacf$SmmCH[datacf$SmmCH < 0.01] <- 1
  datacf$rSmGFP<-datacf$SmGFP/datacf$SmOD
  datacf$rSmmCH<-datacf$SmmCH/datacf$SmOD
  return(datacf)
}



#####################################

setwd("YOUR/WORKING/DIRECTORY")

#import data, blank correct and smoothing
f1<-"170309-PR-dCAS9TDH3"
mat1<-"/170309Mat.txt"
matrix1 <- matrix_import(f1,mat1)
all1 = list.files(path=paste(getwd(),"/",f1,sep = ""), "170309-.*.txt", full = TRUE)
data1 = ldply(all1, import_data, matr=matrix1) #plyr library
datacf1<-blancCorrectratio(data1)
datacfs1<-smouth_func(datacf1)


#general visualisation across the plate
datacfs1_row<-9-as.numeric(factor(datacfs1$row,levels=c("A","B","C","D","E","F","G","H")))
plot(0,0,xlim=c(0,30*12),ylim=c(0,6*8),type = "n",xlab = "",ylab = "OD (log)",yaxt="n",main = "OD")
for (id in levels(as.factor(datacfs1$uniquePos))) {
  points(datacfs1$time[datacfs1$uniquePos == id]/3600+30*(datacfs1$col[datacfs1$uniquePos == id][1]-1),log(datacfs1$SmOD[datacfs1$uniquePos == id])+5+6*(datacfs1_row[datacfs1$uniquePos == id][1]-1),type = "l",col = "black")
}

datacfs1_row<-9-as.numeric(factor(datacfs1$row,levels=c("A","B","C","D","E","F","G","H")))
plot(0,0,xlim=c(0,30*12),ylim=c(0,20000*8),type = "n",xlab = "",ylab = "GFP (ratio)",yaxt="n",main = "GFP")
for (id in levels(as.factor(datacfs1$uniquePos))) {
  points(datacfs1$time[datacfs1$uniquePos == id]/3600+30*(datacfs1$col[datacfs1$uniquePos == id][1]-1),datacfs1$rSmGFP[datacfs1$uniquePos == id]+5+20000*(datacfs1_row[datacfs1$uniquePos == id][1]-1),type = "l",col = "black")
}

datacfs1_row<-9-as.numeric(factor(datacfs1$row,levels=c("A","B","C","D","E","F","G","H")))
plot(0,0,xlim=c(0,30*12),ylim=c(0,2000*8),type = "n",xlab = "",ylab = "mCH (ratio)",yaxt="n",main = "mCH")
for (id in levels(as.factor(datacfs1$uniquePos))) {
  points(datacfs1$time[datacfs1$uniquePos == id]/3600+30*(datacfs1$col[datacfs1$uniquePos == id][1]-1),datacfs1$rSmmCH[datacfs1$uniquePos == id]+5+2000*(datacfs1_row[datacfs1$uniquePos == id][1]-1),type = "l",col = "black")
}

#functions for getting ratio at specific phases of the growth
rowIndicesMidlog <- function(OD, numberOfTimePoints = 5){
  flank <- floor((numberOfTimePoints-1)/2)   # Translate to the amount flanking on each side (prevents even numbers)
  log_OD <- log(OD)
  log_OD[1:2] <- NA # Remove first two point for jumpy data
  log_OD[log_OD < -5] <- NA # Want to remove very small data (inaccurate)
  max_log <- max(log_OD, na.rm = TRUE)
  min_log <- max(log(0.05),min(log_OD, na.rm = TRUE))
  halfmax <- mean(c(max_log, min_log))
  #print(c(exp(max_log),exp(min_log)))
  abs_dist <- abs(log_OD - halfmax)
  max_point <- which(log_OD == max_log)[1]
  abs_dist[max_point:length(OD)] <- NA
  k <- which(abs_dist == min(abs_dist, na.rm = TRUE))[1]
  indices <- (k-flank):(k+flank)
  return(indices)
}
rowIndicesUser <- function(OD, userDefined = 0.3, numberOfTimePoints = 5){
  flank <- floor((numberOfTimePoints-1)/2)   # Translate to the amount flanking on each side (prevents even numbers)
  abs_dist <- abs(OD-userDefined)
  abs_dist[1:2] <- NA
  max_index <- min(which(OD==max(OD[3:length(OD)])))
  abs_dist[max_index:length(OD)] <- NA
  k <- which(abs_dist == min(abs_dist, na.rm = TRUE))[1]
  indices <- (k-flank):(k+flank)
  return(indices)
}
rowIndicesSat <- function(OD, numberOfTimePoints = 5){
  flank <- floor((numberOfTimePoints-1)/2)   # Translate to the amount flanking on each side (prevents even numbers)
  abs_dist <- OD
  abs_dist[1:2] <- NA # remove to avoid jumps
  k <- min(which(OD==max(OD[3:length(OD)], na.rm = TRUE)))[1] # k is timepoint at max OD
  if(k+flank > length(OD)) {k <- length(OD) - flank} # Make sure that range remains within points
  indices <- floor((k-flank):(k+flank))
  return(indices)
}
getGrowthRate <- function(OD, timepoint, initialcut = 4, percentlow = 0.15, percenthigh = 0.10) {
  ODcut<-OD[-c(1:initialcut)]
  timepointcut<-timepoint[-c(1:initialcut)]
  initialOD<-min(ODcut,na.rm = T)
  maxOD<-max(ODcut,na.rm = T)
  LodODcutUP<-log(maxOD)-percenthigh*(log(maxOD)-log(initialOD))
  LodODcutLOW<-log(initialOD)+percentlow*(log(maxOD)-log(initialOD))
  lODcut<-log(ODcut)
  lODcut2<-lODcut[lODcut>LodODcutLOW & lODcut<LodODcutUP]
  timepointcut2<-timepointcut[lODcut>LodODcutLOW & lODcut<LodODcutUP]
  reg<-lm(lODcut2~timepointcut2)
  growthrate<-reg$coefficients[2]*3600
  lag<-(log(initialOD)-reg$coefficients[1])/reg$coefficients[2]/3600
  return(list(growthrate=growthrate,lag=lag,initialOD=initialOD,maxOD=maxOD))
}

#get all fluorescence ratio across the plate at different growth point in a loop
levels(as.factor(datacfs1$uniquePos[!(datacfs1$strain %in% c("Blanc","EMPTY"))]))
wells<-levels(as.factor(datacfs1$uniquePos[!(datacfs1$strain %in% c("Blanc","EMPTY"))]))
Report<-c()
Genetested<-c()
Clone<-c()
Version<-c()
Background<-c()
ODini<-rep(NA,length(wells))
ODmax<-rep(NA,length(wells))
mu<-rep(NA,length(wells))
lag<-rep(NA,length(wells))
ODmid<-rep(NA,length(wells))
GFPmidR<-rep(NA,length(wells))
mCHmidR<-rep(NA,length(wells))
GFPuseR<-rep(NA,length(wells))
mCHuseR<-rep(NA,length(wells))
GFPsatR<-rep(NA,length(wells))
mCHsatR<-rep(NA,length(wells))
ODuse<-rep(NA,length(wells))
ODsat<-rep(NA,length(wells))
for (i in 1:length(wells)) {
  OD<-datacfs1$SmOD[datacfs1$uniquePos==wells[i]]
  timepoint<-datacfs1$time[datacfs1$uniquePos==wells[i]]
  Gdatatemp<-getGrowthRate(OD,timepoint)
  ODini[i]<-Gdatatemp$initialOD
  ODmax[i]<-Gdatatemp$maxOD
  mu[i]<-Gdatatemp$growthrate
  lag[i]<-Gdatatemp$lag
  Report[i]<-as.character(datacfs1$const[datacfs1$uniquePos==wells[i]][1])
  Genetested[i]<-as.character(datacfs1$cond[datacfs1$uniquePos==wells[i]][1])
  Version[i]<-datacfs1$cond2[datacfs1$uniquePos==wells[i]][1]
  Clone[i]<-as.character(datacfs1$clone[datacfs1$uniquePos==wells[i]][1])
  Background[i]<-as.character(datacfs1$strain[datacfs1$uniquePos==wells[i]][1])
  
  if (max(OD)>0.4) {
    ODm<-mean(datacfs1$SmOD[datacfs1$uniquePos==wells[i]][rowIndicesMidlog(OD)])
    GFPm<-mean(datacfs1$SmGFP[datacfs1$uniquePos==wells[i]][rowIndicesMidlog(OD)])
    mCHm<-mean(datacfs1$SmmCH[datacfs1$uniquePos==wells[i]][rowIndicesMidlog(OD)])
    ODu<-mean(datacfs1$SmOD[datacfs1$uniquePos==wells[i]][rowIndicesUser(OD)])
    GFPu<-mean(datacfs1$SmGFP[datacfs1$uniquePos==wells[i]][rowIndicesUser(OD)])
    mCHu<-mean(datacfs1$SmmCH[datacfs1$uniquePos==wells[i]][rowIndicesUser(OD)])
    ODs<-mean(datacfs1$SmOD[datacfs1$uniquePos==wells[i]][rowIndicesSat(OD)])
    GFPs<-mean(datacfs1$SmGFP[datacfs1$uniquePos==wells[i]][rowIndicesSat(OD)])
    mCHs<-mean(datacfs1$SmmCH[datacfs1$uniquePos==wells[i]][rowIndicesSat(OD)])
    ODmid[i]<-ODm
    ODuse[i]<-ODu
    ODsat[i]<-ODs
    GFPmidR[i]<-GFPm/ODm
    mCHmidR[i]<-mCHm/ODm
    GFPuseR[i]<-GFPu/ODu
    mCHuseR[i]<-mCHu/ODu
    GFPsatR[i]<-GFPs/ODs
    mCHsatR[i]<-mCHs/ODs
  }
}
values<-data.frame(wells,Report,Genetested,Version,Clone,Background,
                   ODini,ODmax,mu,lag,
                   ODmid,GFPmidR,mCHmidR,GFPuseR,mCHuseR,GFPsatR,mCHsatR,ODuse,ODsat)


GPD1data<-datacfs1[datacfs1$row=="H" & datacfs1$col==7 & datacfs1$value>0.04 & datacfs1$time<19*3600,]
TDH3data<-datacfs1[datacfs1$row=="G" & datacfs1$col==10 & datacfs1$value>0.04 & datacfs1$time<19*3600,]
TDH3Ndata<-datacfs1[datacfs1$row=="G" & datacfs1$col==11 & datacfs1$value>0.04 & datacfs1$time<19*3600,]
TDH3Ndata<-TDH3Ndata[-1,]

#parsing data in sub table for ploting
datagene<-rbind(GPD1data,TDH3data)

genes<-c("GPD1","TDH3")
type<-c(rep("mid",5),rep("user",5),rep("sat",5))
datagene2<-datagene[datagene$const==genes[1],]
ODv<-datagene2$value
gene<-rep(genes[1],5*3)
time<-c(datagene2$time[rowIndicesMidlog(ODv)],datagene2$time[rowIndicesUser(ODv)],datagene2$time[rowIndicesSat(ODv)])
OD<-c(datagene2$value[rowIndicesMidlog(ODv)],datagene2$value[rowIndicesUser(ODv)],datagene2$value[rowIndicesSat(ODv)])
GFP<-c(datagene2$GFP[rowIndicesMidlog(ODv)],datagene2$GFP[rowIndicesUser(ODv)],datagene2$GFP[rowIndicesSat(ODv)])
mCH<-c(datagene2$mCH[rowIndicesMidlog(ODv)],datagene2$mCH[rowIndicesUser(ODv)],datagene2$mCH[rowIndicesSat(ODv)])
valuesGPD1<-data.frame(gene,type,time,OD,GFP,mCH)

type<-c(rep("mid",5),rep("user",5),rep("sat",5))
datagene2<-datagene[datagene$const==genes[2],]
ODv<-datagene2$value
gene<-rep(genes[2],5*3)
time<-c(datagene2$time[rowIndicesMidlog(ODv)],datagene2$time[rowIndicesUser(ODv)],datagene2$time[rowIndicesSat(ODv)])
OD<-c(datagene2$value[rowIndicesMidlog(ODv)],datagene2$value[rowIndicesUser(ODv)],datagene2$value[rowIndicesSat(ODv)])
GFP<-c(datagene2$GFP[rowIndicesMidlog(ODv)],datagene2$GFP[rowIndicesUser(ODv)],datagene2$GFP[rowIndicesSat(ODv)])
mCH<-c(datagene2$mCH[rowIndicesMidlog(ODv)],datagene2$mCH[rowIndicesUser(ODv)],datagene2$mCH[rowIndicesSat(ODv)])
valuesTDH3<-data.frame(gene,type,time,OD,GFP,mCH)

type<-c(rep("mid",5),rep("user",5),rep("sat",5))
datagene2<-TDH3Ndata
ODv<-datagene2$value
gene<-rep("TDH3N",5*3)
time<-c(datagene2$time[rowIndicesMidlog(ODv)],datagene2$time[rowIndicesUser(ODv)],datagene2$time[rowIndicesSat(ODv)])
OD<-c(datagene2$value[rowIndicesMidlog(ODv)],datagene2$value[rowIndicesUser(ODv)],datagene2$value[rowIndicesSat(ODv)])
GFP<-c(datagene2$GFP[rowIndicesMidlog(ODv)],datagene2$GFP[rowIndicesUser(ODv)],datagene2$GFP[rowIndicesSat(ODv)])
mCH<-c(datagene2$mCH[rowIndicesMidlog(ODv)],datagene2$mCH[rowIndicesUser(ODv)],datagene2$mCH[rowIndicesSat(ODv)])
valuesTDH3N<-data.frame(gene,type,time,OD,GFP,mCH)


#now ploting (curve and barplot)

plot(0,0,xlim=c(0,20),ylim=c(0,0.7),type='n',main="GPD1",xlab="time (h)", ylab="OD (log)")
points(GPD1data$time/3600,GPD1data$value,cex=0.7)
points(valuesGPD1$time/3600,valuesGPD1$OD,pch=16)

svg(filename = "GPD1-OD.svg",width = 4.3,height = 4)
plot(0,0,xlim=c(0,20),ylim=c(-3.2,0),type='n',main="GPD1",xlab="time (h)", ylab="OD (log scale)",yaxt="n")
points(GPD1data$time/3600,log(GPD1data$value))
points(valuesGPD1$time[valuesTDH3$type=="user"]/3600,log(valuesGPD1$OD)[valuesTDH3$type=="user"],pch=16)
axis(2,at = log(c(0.05,0.1,0.2,0.4,0.8)),labels = c(0.05,0.1,0.2,0.4,0.8))
dev.off()

svg(filename = "GPD1-GFP.svg",width = 4.3,height = 4)
plot(0,0,xlim=c(0,20),ylim=c(0,3000),type='n',main="GPD1",xlab="time (h)", ylab="GFP (ratio)")
points(GPD1data$time/3600,GPD1data$GFP,col='#34ab83ff')
points(valuesGPD1$time[valuesTDH3$type=="user"]/3600,valuesGPD1$GFP[valuesTDH3$type=="user"],pch=16,col='#34ab83ff')
dev.off()

svg(filename = "GPD1-mCH.svg",width = 4.3,height = 4)
plot(0,0,xlim=c(0,20),ylim=c(0,1000),type='n',main="GPD1",xlab="time (h)", ylab="mCH (ratio)")
points(GPD1data$time/3600,GPD1data$mCH,col='red')
points(valuesGPD1$time[valuesTDH3$type=="user"]/3600,valuesGPD1$mCH[valuesTDH3$type=="user"],pch=16,col='red')
dev.off()

valuesTDH3F<-rbind(valuesTDH3,valuesTDH3N,valuesGPD1)
valuesTDH3F$GFPr<-valuesTDH3F$GFP/valuesTDH3F$OD
valuesTDH3F$mCHr<-valuesTDH3F$mCH/valuesTDH3F$OD

valuesTDH3F2<-valuesTDH3F[valuesTDH3F$type=="user",]

plot(0,0,xlim=c(0,20),ylim=c(0,0.7),type='n',main="TDH3",xlab="time (h)", ylab="OD")
points(TDH3data$time/3600,TDH3data$value,pch=1,cex=0.7,col='black')
points(TDH3Ndata$time/3600,TDH3Ndata$value,pch=1,cex=0.4,col='black')
points(GPD1data$time/3600,GPD1data$value,pch=2,cex=0.7,col='black')
points(valuesTDH3F2$time/3600,valuesTDH3F2$OD,pch=16,cex=1,col='black')

svg(filename = "TDH3-OD.svg",width = 4.3,height = 4)
plot(0,0,xlim=c(0,20),ylim=c(-3.2,0),type='n',main="TDH3",xlab="time (h)", ylab="OD (log scale)",yaxt="n")
points(TDH3data$time/3600,log(TDH3data$value),pch=1,cex=0.5,col='black')
points(TDH3Ndata$time/3600,log(TDH3Ndata$value),pch=3,cex=0.5,col='black')
points(GPD1data$time/3600,log(GPD1data$value),pch=2,cex=0.5,col='black')
points(valuesTDH3F2$time/3600,log(valuesTDH3F2$OD),pch=16,cex=0.8,col='black')
axis(2,at = log(c(0.05,0.1,0.2,0.4,0.8)),labels = c(0.05,0.1,0.2,0.4,0.8))
legend("bottomright",pch=c(1,1,2),pt.cex=(c(0.6,0.4,0.6)),legend = c("TDH3-GFP-gRNA","TDH3-GFP","GPD1-GFP-gRNA"))
dev.off()

svg(filename = "TDH3-GFP.svg",width = 4.3,height = 4)
plot(0,0,xlim=c(0,20),ylim=c(0,35000),type='n',main="TDH3",xlab="time (h)", ylab="GFP")
points(TDH3data$time/3600,TDH3data$GFP,pch=1,cex=0.5,col='#34ab83ff')
points(TDH3Ndata$time/3600,TDH3Ndata$GFP,pch=3,cex=0.5,col='#34ab83ff')
points(GPD1data$time/3600,GPD1data$GFP,pch=2,cex=0.5,col='#34ab83ff')
points(valuesTDH3F2$time/3600,valuesTDH3F2$GFP,pch=16,cex=0.8,col='#34ab83ff')
dev.off()

svg(filename = "TDH3-mCH.svg",width = 4.3,height = 4)
plot(0,0,xlim=c(0,20),ylim=c(0,1000),type='n',main="TDH3",xlab="time (h)", ylab="mCherry")
points(TDH3data$time/3600,TDH3data$mCH,pch=1,cex=0.5,col='red')
points(TDH3Ndata$time/3600,TDH3Ndata$mCH,pch=3,cex=0.5,col='red')
points(GPD1data$time/3600,GPD1data$mCH,pch=2,cex=0.5,col='red')
points(valuesTDH3F2$time/3600,valuesTDH3F2$mCH,pch=16,cex=0.8,col='red')
dev.off()

svg(filename = "BAR-TDH3-GFP.svg",width = 2.3,height = 4)
p<-barplot(height = c(mean(valuesTDH3F2$GFPr[valuesTDH3F2$gene=="TDH3"]),mean(valuesTDH3F2$GFPr[valuesTDH3F2$gene=="TDH3N"]),mean(valuesTDH3F2$GFPr[valuesTDH3F2$gene=="GPD1"]))
           ,names.arg = c("TDH3","TDH3N","GPD1"),ylim=c(0,100000),col='#34ab83ff')
points(rep(p,5)[order(rep(p,5))],valuesTDH3F2$GFPr,pch=16,cex=0.8)
dev.off()

svg(filename = "BAR-TDH3-mCH.svg",width = 2.3,height = 4)
p<-barplot(height = c(mean(valuesTDH3F2$mCHr[valuesTDH3F2$gene=="TDH3"]),mean(valuesTDH3F2$mCHr[valuesTDH3F2$gene=="TDH3N"]),mean(valuesTDH3F2$mCHr[valuesTDH3F2$gene=="GPD1"]))
           ,names.arg = c("TDH3","TDH3N","GPD1"),ylim=c(0,2000),col='red')
points(rep(p,5)[order(rep(p,5))],valuesTDH3F2$mCHr,pch=16,cex=0.8)
dev.off()
