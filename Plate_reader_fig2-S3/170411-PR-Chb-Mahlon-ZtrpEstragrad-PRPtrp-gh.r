#this R script reads plate reader data and process growth curve and fluorescence ratio for tagged strains (GFP and mCherry)
#Z3EV experiment with increase estradio concentration (data from 2017-04-11)

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
f1<-"folder_name"

#import data, blank correct and smoothing
mat1<-"/170412MAT.txt"
matrix1 <- matrix_import(f1,mat1)
all1 = list.files(path=paste(getwd(),"/",f1,sep = ""), "170412-.*.txt", full = TRUE)
data1 = ldply(all1, import_data, matr=matrix1) #plyr library
datacf1<-blancCorrectratio(data1)
datacf1$cond2<-as.character(datacf1$cond2)
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

#to lose the last two bars:
datacfs1<-datacfs1[!(datacfs1$cond2 %in% c("1000","500","250")),]

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
rowIndicesUser <- function(OD, userDefined = 0.2, numberOfTimePoints = 5){
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
  
  if (max(OD)>0.2) {
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


valuesBYwt<-values[values$Background=="BY" & values$Genetested=="gwT",]
valuesBYwt<-valuesBYwt[order(as.numeric(as.character(valuesBYwt$Version))),]
valuesBYwt<-valuesBYwt[valuesBYwt$Version!="1000",]

vRGBini=c(124,231,229)
vRGB=t(vRGBini%*%(t((10:0)+10)/20))/255
colestra=rgb(vRGB)

#now ploting (barplots)
pl<-barplot(valuesBYwt$GFPuseR,col="#34ab83ff")
well<-c()
color<-c()
xval<-c()
time<-c()
OD<-c()
GFP<-c()
mCH<-c()
for (i in 1:nrow(valuesBYwt)) {
  ODv<-datacfs1$SmOD[datacfs1$uniquePos==valuesBYwt$wells[i]]
  well<-c(well,rep(as.character(valuesBYwt$wells[i]),5))
  color<-c(color,rep(colestra[i],5))
  xval<-c(xval,rep(pl[i],5))
  time<-c(time,datacfs1$time[datacfs1$uniquePos==valuesBYwt$wells[i]][rowIndicesUser(ODv)])
  OD<-c(OD,datacfs1$value[datacfs1$uniquePos==valuesBYwt$wells[i]][rowIndicesUser(ODv)])
  GFP<-c(GFP,datacfs1$GFP[datacfs1$uniquePos==valuesBYwt$wells[i]][rowIndicesUser(ODv)])
  mCH<-c(mCH,datacfs1$mCH[datacfs1$uniquePos==valuesBYwt$wells[i]][rowIndicesUser(ODv)])
  plot(datacfs1$time[datacfs1$uniquePos==valuesBYwt$wells[i]],datacfs1$value[datacfs1$uniquePos==valuesBYwt$wells[i]])
  points(datacfs1$time[datacfs1$uniquePos==valuesBYwt$wells[i]][rowIndicesUser(ODv)],datacfs1$value[datacfs1$uniquePos==valuesBYwt$wells[i]][rowIndicesUser(ODv)],pch=16,col=colestra[i])
}
pointvalue<-data.frame(well,color,xval,time,OD,GFP,mCH)
pointvalue$color<-as.character(pointvalue$color)
pointvalue$GFPr<-pointvalue$GFP/pointvalue$OD
pointvalue$mCHr<-pointvalue$mCH/pointvalue$OD


svg(filename = "BYestraGFP.svg",width = 7,height = 3.7)
barplot(valuesBYwt$GFPuseR,col="#34ab83ff",ylim=c(0,80000))
points(pointvalue$xval,pointvalue$GFPr,pch=16)
dev.off()

svg(filename = "BYestramCH.svg",width = 7,height = 3.7)
barplot(valuesBYwt$mCHuseR,col="red",ylim=c(0,2500))
points(pointvalue$xval,pointvalue$mCHr,pch=16)
dev.off()

valuesRMwt<-values[values$Background=="RM" & values$Genetested=="gwT",]
valuesRMwt<-valuesRMwt[order(as.numeric(as.character(valuesRMwt$Version))),]
valuesRMwt<-valuesRMwt[valuesRMwt$Version!="1000",]

vRGBini=c(124,231,229)
vRGB=t(vRGBini%*%(t((10:0)+10)/20))/255
colestra=rgb(vRGB)

pl<-barplot(valuesRMwt$GFPuseR,col="#34ab83ff")
well<-c()
color<-c()
xval<-c()
time<-c()
OD<-c()
GFP<-c()
mCH<-c()
for (i in 1:nrow(valuesRMwt)) {
  ODv<-datacfs1$SmOD[datacfs1$uniquePos==valuesRMwt$wells[i]]
  well<-c(well,rep(as.character(valuesRMwt$wells[i]),5))
  color<-c(color,rep(colestra[i],5))
  xval<-c(xval,rep(pl[i],5))
  time<-c(time,datacfs1$time[datacfs1$uniquePos==valuesRMwt$wells[i]][rowIndicesUser(ODv)])
  OD<-c(OD,datacfs1$value[datacfs1$uniquePos==valuesRMwt$wells[i]][rowIndicesUser(ODv)])
  GFP<-c(GFP,datacfs1$GFP[datacfs1$uniquePos==valuesRMwt$wells[i]][rowIndicesUser(ODv)])
  mCH<-c(mCH,datacfs1$mCH[datacfs1$uniquePos==valuesRMwt$wells[i]][rowIndicesUser(ODv)])
  plot(datacfs1$time[datacfs1$uniquePos==valuesRMwt$wells[i]],datacfs1$value[datacfs1$uniquePos==valuesRMwt$wells[i]])
  points(datacfs1$time[datacfs1$uniquePos==valuesRMwt$wells[i]][rowIndicesUser(ODv)],datacfs1$value[datacfs1$uniquePos==valuesRMwt$wells[i]][rowIndicesUser(ODv)],pch=16,col=colestra[i])
}
pointvalueRM<-data.frame(well,color,xval,time,OD,GFP,mCH)
pointvalueRM$color<-as.character(pointvalueRM$color)
pointvalueRM$GFPr<-pointvalueRM$GFP/pointvalueRM$OD
pointvalueRM$mCHr<-pointvalueRM$mCH/pointvalueRM$OD

svg(filename = "RMestraGFP.svg",width = 7,height = 3.7)
barplot(valuesRMwt$GFPuseR,col="#34ab83ff",ylim=c(0,80000))
points(pointvalueRM$xval,pointvalueRM$GFPr,pch=16)
dev.off()

svg(filename = "RMestramCH.svg",width = 7,height = 3.7)
barplot(valuesRMwt$mCHuseR,col="red",ylim=c(0,2500))
points(pointvalueRM$xval,pointvalueRM$mCHr,pch=16)
dev.off()
