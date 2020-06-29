#this script draw circulare and linear plasmid map from tabulated map text file and coordonate (for linear draw)
library("PlotRegionHighlighter")

myWD<-"YOUR/WORK/DIRECTORY" #to change
setwd(myWD)

rad<-1
wid<-c(0.1,0.065)
names(wid)<-c("gene","other")

#property across plasmid
plasmidsUsed<-c("pTAG-leu","pdCAS9","pZ3EV","pRS415-mCH","pRS425-mCH","pTRP-CEN-dCAS9","pTRP-2u-dCAS9") #all the tab file of the plasmid map
lineareP1<-c(5727,11046,1862,2,2,2,2)
lineareP2<-c(8741,11044,7019,1,1,1,1)

#input for each plasmid
plval<-1
liP1<-lineareP1[plval] #lineare coordonate start
liP2<-lineareP2[plval] #lineare coordonate end
pl<-plasmidsUsed[plval] #plasmid name

#reformat the plasmid data
plasmidINFO<-read.table(pl,header = T,quote = "",sep = "\t")

plasmidINFOtot<-plasmidINFO[plasmidINFO$Type=="plasmid",]
plasmidINFO2<-plasmidINFO[plasmidINFO$Type!="plasmid",]
plasmidINFO2$type2<-rep("other",nrow(plasmidINFO2))
plasmidINFO2$type2[plasmidINFO2$Type=="gene"]<-"gene"
plasmidINFO2$aStart<- -((plasmidINFO2$Start/plasmidINFOtot$End[1])*(2*pi))+(0.5*pi)
plasmidINFO2$aStop<- -((plasmidINFO2$End/plasmidINFOtot$End[1])*(2*pi))+(0.5*pi)
plasmidINFO2<-plasmidINFO2[order(plasmidINFO2$Start),]
plasmidINFO2$col<-rep(c("grey","darkgrey"),ceiling(nrow(plasmidINFO2)/2))[1:nrow(plasmidINFO2)]
plasmidINFO2$col[plasmidINFO2$Type=="gene"]<-"green"
plasmidINFO2$col[plasmidINFO2$Type=="terminater"]<-"cyan"
plasmidINFO2$col[plasmidINFO2$Type=="promoter"]<-"lightblue"
plasmidINFO2$n<-ceiling(4000*((plasmidINFO2$End-plasmidINFO2$Start)/plasmidINFOtot$End[1]))

#draw circulare plasmid (save as svg file, can be edited on Inkscape)
svg(filename = paste(pl,"_circ.svg",sep=""),width = 6,height = 6.5)
plot(0,0,typ="n",xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),main=paste(plasmidINFOtot$Name[1],plasmidINFOtot$End[1]))
lines(createCircle(c(0,0), r=1, n = 1000, begin = 0, end = 2*pi),lwd=3)
for (i in 1:nrow(plasmidINFO2)) {
  if (plasmidINFO2$direct[i]=="Forward") {
    Cint<-createCircle(c(0,0), r=1-wid[plasmidINFO2$type2[i]], n = plasmidINFO2$n[i], begin = plasmidINFO2$aStart[i], end = plasmidINFO2$aStop[i])
    Cext<-createCircle(c(0,0), r=1+wid[plasmidINFO2$type2[i]], n = plasmidINFO2$n[i], begin = plasmidINFO2$aStart[i], end = plasmidINFO2$aStop[i])
    Carrow<-createCircle(c(0,0), r=1, n = 1, begin = plasmidINFO2$aStop[i], end = plasmidINFO2$aStop[i])
    arroB<-min(90,floor(plasmidINFO2$n[i]/2))
    obj<-rbind(Cint[c(1:(nrow(Cint)-arroB)),],Carrow,Cext[c((nrow(Cext)-floor(arroB*0.9)):1),])
  } else if (plasmidINFO2$direct[i]=="Reverse") {
    Cint<-createCircle(c(0,0), r=1-wid[plasmidINFO2$type2[i]], n = plasmidINFO2$n[i], begin = plasmidINFO2$aStart[i], end = plasmidINFO2$aStop[i])
    Cext<-createCircle(c(0,0), r=1+wid[plasmidINFO2$type2[i]], n = plasmidINFO2$n[i], begin = plasmidINFO2$aStart[i], end = plasmidINFO2$aStop[i])
    Carrow<-createCircle(c(0,0), r=1, n = 1, begin = plasmidINFO2$aStart[i], end = plasmidINFO2$aStart[i])
    arroB<-min(90,floor(plasmidINFO2$n[i]/2))
    obj<-rbind(Cint[c(arroB:(nrow(Cint))),],Cext[c((nrow(Cext)):floor(arroB*0.9)),],Carrow)
  } else {
    Cint<-createCircle(c(0,0), r=1-wid[plasmidINFO2$type2[i]], n = plasmidINFO2$n[i], begin = plasmidINFO2$aStart[i], end = plasmidINFO2$aStop[i])
    Cext<-createCircle(c(0,0), r=1+wid[plasmidINFO2$type2[i]], n = plasmidINFO2$n[i], begin = plasmidINFO2$aStart[i], end = plasmidINFO2$aStop[i])
    obj<-rbind(Cint[c(1:(nrow(Cint))),],Cext[c((nrow(Cext)):1),])
  }
  polygon(obj,col=plasmidINFO2$col[i])
}
dev.off()


#draw linear fragment (save as svg file, can be edited on Inkscape)
linearFrom<-c(liP1,liP2)
plasmidINFO22<-plasmidINFO2
if (linearFrom[2]<linearFrom[1]) {
  plasmidINFO22$End<-plasmidINFO2$End+plasmidINFOtot$End[1]
  plasmidINFO22$Start<-plasmidINFO2$Start+plasmidINFOtot$End[1]
  plasmidINFO22<-rbind(plasmidINFO2,plasmidINFO22)
  linearFrom[2]<-linearFrom[2]+plasmidINFOtot$End[1]
}
plasmidINFO3<-plasmidINFO22[(plasmidINFO22$Start>linearFrom[1] & plasmidINFO22$Start<linearFrom[2]) | (plasmidINFO22$End>linearFrom[1] & plasmidINFO22$End<linearFrom[2]),]
plasmidINFO3$direct[plasmidINFO3$Start<linearFrom[1] & plasmidINFO3$direct=="Reverse"] <- "no"
plasmidINFO3$Start[plasmidINFO3$Start<linearFrom[1]]<-linearFrom[1]
plasmidINFO3$direct[plasmidINFO3$End>linearFrom[2] & plasmidINFO3$direct=="Forward"] <- "no"
plasmidINFO3$End[plasmidINFO3$End>linearFrom[2]]<-linearFrom[2]
plasmidINFO3$Start<-plasmidINFO3$Start-linearFrom[1]
plasmidINFO3$End<-plasmidINFO3$End-linearFrom[1]


svg(filename = paste(pl,"_lin.svg",sep=""),width = 15,height = 8)
plot(0,0,typ="n",xlim=c(0,13000),ylim=c(-1.4,1.4),main=paste(plasmidINFOtot$Name[1],plasmidINFOtot$End[1]))
lines(linearFrom-linearFrom[1],c(0,0),lwd=3, lend="butt")
for (i in 1:nrow(plasmidINFO3)) {
  if (plasmidINFO3$direct[i]=="Forward") {
    arroB<-min(250,floor((plasmidINFO3$End[i]-plasmidINFO3$Start[i])/2))
    x<-c(plasmidINFO3$Start[i],plasmidINFO3$End[i]-arroB,plasmidINFO3$End[i],plasmidINFO3$End[i]-arroB,plasmidINFO3$Start[i])
    y<-c(wid[plasmidINFO3$type2[i]],wid[plasmidINFO3$type2[i]],0,-wid[plasmidINFO3$type2[i]],-wid[plasmidINFO3$type2[i]])
    obj<-data.frame(x,y)
  } else if (plasmidINFO3$direct[i]=="Reverse") {
    arroB<-min(250,floor((plasmidINFO3$End[i]-plasmidINFO3$Start[i])/2))
    x<-c(plasmidINFO3$End[i],plasmidINFO3$Start[i]+arroB,plasmidINFO3$Start[i],plasmidINFO3$Start[i]+arroB,plasmidINFO3$End[i])
    y<-c(wid[plasmidINFO3$type2[i]],wid[plasmidINFO3$type2[i]],0,-wid[plasmidINFO3$type2[i]],-wid[plasmidINFO3$type2[i]])
    obj<-data.frame(x,y)
  } else {
    x<-c(plasmidINFO3$Start[i],plasmidINFO3$End[i],plasmidINFO3$End[i],plasmidINFO3$Start[i])
    y<-c(wid[plasmidINFO3$type2[i]],wid[plasmidINFO3$type2[i]],-wid[plasmidINFO3$type2[i]],-wid[plasmidINFO3$type2[i]])
    obj<-data.frame(x,y)
  }
  polygon(obj,col=plasmidINFO3$col[i])
}
dev.off()

linearFrom-linearFrom[1] #linear size