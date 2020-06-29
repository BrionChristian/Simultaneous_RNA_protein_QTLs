#this R script reads plate reader data and process growth curve and fluorescence ratio for tagged strains (GFP and mCherry)
#CRISPR swap speficicaly of tileA wit BY GPD1-GFP DNA (data from 2019-06-19)

library(reshape2)
library(plyr)
par(mfrow = c(1, 1))

#functions used for reading plate reader data
timeC <- function(Strtime) { #trasnform time value into numeric (hours)
  val<-as.numeric(strsplit(Strtime,":")[[1]]) %*% c(3600,60,1)
  val
}

import_data <- function(file,matr){ #import read data and merge them with matrix information
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
  BlancOD=mean(data$value[data$strain=="Blanc" & data$type=="OD"])
  BlancGFP=min(data$value[data$strain!="EMPTY" & data$type=="GFP"])
  BlancmCH=mean(data$value[data$strain=="Blanc" & data$type=="mCH"])
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
  matrix <- cbind(matrix,t(data.frame(lapply(strsplit(matrix$value,"[-,_]"), function(u) transform(c(u,rep(NA,5-length(u))))))))
  matrix <- matrix[,-3]
  colnames(matrix) <- c("row","col","strain","const","cond","cond2","clone")
  matrix1<-matrix[order(matrix$row),]
}

smouth_func <- function(datacf) { #smooth curves to reduce point-to-point noise
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
f1<-"FOLDER_NAME"

#import data, blank correct and smoothing
mat1<-"/190619Mat.txt"
matrix1 <- matrix_import(f1,mat1)
all1 = list.files(path=paste(getwd(),"/",f1,sep = ""), "190619-.*.txt", full = TRUE)
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


datacfs<-datacfs1[datacfs1$strain=="BY" & datacfs1$cond != "NA",]
datacfs1<-datacfs1[datacfs1$strain=="BY" & datacfs1$cond != "NA",]


#generate PDF to compare curves
pdf(file = "190619-ChB-FM-GPD1-FM10aBYRMGG2.pdf",width = 10,height = 9)

genes<-c("FM10a")

for (j in 1:length(genes)) {
  datacfsT<-datacfs[datacfs$cond==genes[j],]
  levels(as.factor(datacfsT$uniquePos))
  
  plot(0,0,xlim=c(0,18),ylim=c(-4,1),type='n',main=genes[j],xlab="time (h)", ylab="OD (log)")
  for (i in levels(as.factor(datacfsT$uniquePos))) {
    if (length(grep(pattern = "_GG_",x = i))>0) {tempcol<-"grey"}
    if (length(grep(pattern = "_RM_",x = i))>0) {tempcol<-"red"}
    if (length(grep(pattern = "_BY_",x = i))>0) {tempcol<-"blue"}
    points(datacfsT$time[datacfsT$uniquePos==i]/3600,log(datacfsT$SmOD[datacfsT$uniquePos==i]),type="l",lwd=1.5,col=tempcol)
  }
  
  plot(0,0,xlim=c(0,18),ylim=c(0,10000),type='n',main=genes[j],xlab="time (h)", ylab="OD (log)")
  for (i in levels(as.factor(datacfsT$uniquePos))) {
    if (length(grep(pattern = "_GG_",x = i))>0) {tempcol<-"grey"}
    if (length(grep(pattern = "_RM_",x = i))>0) {tempcol<-"red"}
    if (length(grep(pattern = "_BY_",x = i))>0) {tempcol<-"blue"}
    points(datacfsT$time[datacfsT$uniquePos==i]/3600,datacfsT$rSmGFP[datacfsT$uniquePos==i],type="l",lwd=1.5,col=tempcol)
  }
  
  typ<-c("RM","BY","GG")
  colt<-c(rgb(1,0,0,1),rgb(0,0,1,1),rgb(0.4,0.4,0.4,1))
  colt2<-c(rgb(1,0,0,0.2),rgb(0,0,1,0.2),rgb(0.4,0.4,0.4,0.2))
  nSE=3
  plot(0,0,xlim=c(0,18),ylim=c(-4,1),type='n',main=genes[j],xlab="time (h)", ylab="GFP (ratio)")
  for (i in 1:3) {
    datacfsT<-datacfs[datacfs$cond==genes[j] & datacfs$cond2==typ[i],]
    tempcol<-colt[i]
    tempcol2<-colt2[i]
    x<-datacfsT$time/3600
    y<-log(datacfsT$SmOD)
    model <- loess(y ~ x,span=0.2)
    fitted <- predict(model, se = T)
    xx <- c(x, rev(x))
    yy <- c(fitted$fit+nSE*fitted$se.fit, rev(fitted$fit-nSE*fitted$se.fit))
    polygon(xx, yy, border = NA,col=tempcol2)
    lines(x, fitted$fit+nSE*fitted$se.fit, lty = "dotted",col=tempcol)
    lines(x, fitted$fit-nSE*fitted$se.fit, lty = "dotted",col=tempcol)
    lines(x, fitted$fit,col=tempcol)
  }
  
  
  plot(0,0,xlim=c(0,18),ylim=c(0,10000),type='n',main=genes[j],xlab="time (h)", ylab="GFP (ratio)")
  for (i in 1:3) {
    datacfsT<-datacfs[datacfs$cond==genes[j] & datacfs$cond2==typ[i],]
    tempcol<-colt[i]
    tempcol2<-colt2[i]
    x<-datacfsT$time/3600
    y<-datacfsT$rSmGFP
    model <- loess(y ~ x,span=0.2)
    fitted <- predict(model, se = T)
    xx <- c(x, rev(x))
    yy <- c(fitted$fit+nSE*fitted$se.fit, rev(fitted$fit-nSE*fitted$se.fit))
    polygon(xx, yy, border = NA,col=tempcol2)
    lines(x, fitted$fit+nSE*fitted$se.fit, lty = "dotted",col=tempcol)
    lines(x, fitted$fit-nSE*fitted$se.fit, lty = "dotted",col=tempcol)
    lines(x, fitted$fit,col=tempcol)
  }
}

##################################
#Determine Time Points at Mid Log

#functions for getting ratio at specific phases of the growth
rowIndicesMidlog <- function(OD, numberOfTimePoints = 5){
  flank <- floor((numberOfTimePoints-1)/2)   # Translate to the amount flanking on each side (prevents even numbers)
  log_OD <- log(OD)
  log_OD[1:2] <- NA # Remove first two point for jumpy data
  log_OD[log_OD < -5] <- NA # Want to remove very small data (inaccurate)
  max_log <- max(log_OD, na.rm = TRUE)
  min_log <- max(log(0.05),min(log_OD, na.rm = TRUE))
  halfmax <- mean(c(max_log, min_log))
  print(c(exp(max_log),exp(min_log)))
  abs_dist <- abs(log_OD - halfmax)
  max_point <- which(log_OD == max_log)[1]
  abs_dist[max_point:length(OD)] <- NA
  k <- which(abs_dist == min(abs_dist, na.rm = TRUE))[1]
  indices <- (k-flank):(k+flank)
  return(indices)
}
rowIndicesUser <- function(OD, userDefined = 0.4, numberOfTimePoints = 5){
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

#get all fluorescence ratio across the plate at different growth point in a loop
levels(as.factor(datacfs1$uniquePos[datacfs1$strain!="Blanc"]))
Report<-c()
Genetested<-c()
Clone<-c()
Version<-c()
Background<-c()
ODmid<-c()
GFPmidR<-c()
mCHmidR<-c()
GFPuseR<-c()
mCHuseR<-c()
GFPsatR<-c()
mCHsatR<-c()
ODuse<-c()
ODsat<-c()
wells<-levels(as.factor(datacfs1$uniquePos[datacfs1$strain!="Blanc"]))
for (i in 1:length(wells)) {
  OD<-datacfs1$SmOD[datacfs1$uniquePos==wells[i]]
  ODm<-mean(datacfs1$SmOD[datacfs1$uniquePos==wells[i]][rowIndicesMidlog(OD)])
  GFPm<-mean(datacfs1$SmGFP[datacfs1$uniquePos==wells[i]][rowIndicesMidlog(OD)])
  mCHm<-mean(datacfs1$SmmCH[datacfs1$uniquePos==wells[i]][rowIndicesMidlog(OD)])
  ODu<-mean(datacfs1$SmOD[datacfs1$uniquePos==wells[i]][rowIndicesUser(OD)])
  GFPu<-mean(datacfs1$SmGFP[datacfs1$uniquePos==wells[i]][rowIndicesUser(OD)])
  mCHu<-mean(datacfs1$SmmCH[datacfs1$uniquePos==wells[i]][rowIndicesUser(OD)])
  ODs<-mean(datacfs1$SmOD[datacfs1$uniquePos==wells[i]][rowIndicesSat(OD)])
  GFPs<-mean(datacfs1$SmGFP[datacfs1$uniquePos==wells[i]][rowIndicesSat(OD)])
  mCHs<-mean(datacfs1$SmmCH[datacfs1$uniquePos==wells[i]][rowIndicesSat(OD)])

  Report[i]<-as.character(datacfs1$const[datacfs1$uniquePos==wells[i]][1])
  Genetested[i]<-as.character(datacfs1$cond[datacfs1$uniquePos==wells[i]][1])
  Version[i]<-as.character(datacfs1$cond2[datacfs1$uniquePos==wells[i]][1])
  Clone[i]<-as.character(datacfs1$clone[datacfs1$uniquePos==wells[i]][1])
  Background[i]<-as.character(datacfs1$strain[datacfs1$uniquePos==wells[i]][1])
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
values<-data.frame(Report,Genetested,Version,Clone,Background,ODmid,GFPmidR,mCHmidR,GFPuseR,mCHuseR,GFPsatR,mCHsatR,ODuse,ODsat)
values$test<-as.factor(paste(values$Genetested,values$Version,sep="_"))
values2<-values[order(values$test),]
values2$Clone<-as.numeric(as.character(values2$Clone))%%10
values2$test2<-paste(values2$test,values2$Clone,sep="-")
values3<-values2[values2$Genetested!="NA",]
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)

genes<-c("FM10a")

limHIGH <- 10000
plot(values2$GFPmidR~values2$test,ylim = c(0,limHIGH),pch=16,cex=0.6,col="green")
text(as.numeric(values2$test)+as.numeric(as.character(values2$Clone))*0.07-0.4,values2$GFPmidR,labels = values2$Clone,cex=0.7, col='black')

plot(values2$GFPuseR~values2$test,pch=16,cex=0.6,col="green") #,ylim = c(0,limHIGH)
text(as.numeric(values2$test)+as.numeric(as.character(values2$Clone))*0.07-0.4,values2$GFPuseR,labels = values2$Clone,cex=0.7, col='black')

plot(values2$GFPsatR~values2$test,ylim = c(0,limHIGH),pch=16,cex=0.6,col="green")
text(as.numeric(values2$test)+as.numeric(as.character(values2$Clone))*0.07-0.4,values2$GFPsatR,labels = values2$Clone,cex=0.7, col='black')

dev.off()

#now ploting (boxplot)
values3dc<-dcast(values2,formula = Version~Clone,fun.aggregate = mean,value.var = "GFPuseR")
values3<-values2
values3$test<-factor(as.character(values3$test),levels = c("FM10a_BY","FM10a_RM","FM10a_GG"))
values3$colplot<-rep("#0000ffff",nrow(values3))
values3$colplot[values3$Version=="RM"]<-"#e3d905ff"
values3$colplot[values3$Version=="GG"]<-"#2f2f7dff"

limHIGH<-7000
svg(filename = "chrX_a_GG_swaptest1.svg",width = 3.65,height = 4.5)
boxplot(t(values3dc[c(1,3,2),-1]),pch=1,cex=0.4,col="#00CC80",ylim = c(3000,limHIGH))
for(i in 1:8) {
  for (j in 1:3) {
    values3temp<-values3[values3$test==levels(values3$test)[j] & values3$Clone==i,]
    if (nrow(values3temp)==2) {
      points(rep(as.numeric(values3temp$test[1])+as.numeric(as.character(values3temp$Clone[1]))*0.07-0.315,2),values3temp$GFPuseR,cex=0.7,pch=16,col=values3temp$colplot[1],typ="l")
    }
  }
}
points(as.numeric(values3$test)+as.numeric(as.character(values3$Clone))*0.07-0.315,values3$GFPuseR,cex=0.7,pch=16,col=values3$colplot)
points(as.numeric(values3$test)+as.numeric(as.character(values3$Clone))*0.07-0.315,values3$GFPuseR,cex=0.7,pch=1,col='black')
dev.off()

region<-levels(values3$Version)[c(1,3,2)]

fuseval<-t(values3dc[c(1,3,2),-1])
fusevallog<-log2(fuseval)
test<-t.test(c(fuseval[,1],fuseval[,2]),fuseval[,3])
print(paste(log2(test$estimate[2]/test$estimate[1]),test$p.value))

#same with log2
test<-t.test(c(fusevallog[,1],fusevallog[,2]),fusevallog[,3])
print(paste(test$estimate[2]-test$estimate[1],test$p.value))



