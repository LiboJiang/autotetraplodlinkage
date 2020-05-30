
chr <- read.csv("Marker_info.csv",header=T)


load("F.RData")




dat <- filtered_data
dat[which(dat==0)] <- 1111
dat[which(dat==1)] <- 1112
dat[which(dat==2)] <- 1122
dat[which(dat==3)] <- 1222
dat[which(dat==4)] <- 2222


#01  02  10  11  12  13  20  21  22 
#314 108 327 275 180  58 138 235 135 



source("auto_sim.R")
source("auto_est.R")
source("auto-debug.R")



#flg <- c()
#for(i in 1:12){
  
#  ii <- which(chr$﻿CHR==paste('CHR',i,sep=" "))
#  d1 <- match(chr$Marker[ii],rownames(filtered_data))
#  dd1 <- na.omit(d1)
#  nd <- length(dd1)
  
#  rc <- 0
#  for(j in 1:(nd-1)){
#    Two_res <- auto_scan(dat=dat,se=c(dd1[j],dd1[j+1]))
#    rc <- c(rc,as.numeric(Two_res[1,9]))
#  }
  
#  rtm <- cbind(rep(i,nd),rownames(filtered_data)[dd1],rc)
#  flg <- rbind(flg,rtm)
#}



#se <- which(paste(filtered_data[,1],filtered_data[,2],sep="")=="21")[6:7]

#Two_res <- auto_scan(dat=dat,se=se)


#write.csv(Two_res,file="res.csv")

i <- 1

ii <- which(chr$﻿CHR==paste('CHR',i,sep=" "))
d1 <- match(chr$Marker[ii],rownames(filtered_data))
dd1 <- na.omit(d1)
nd <- length(dd1)

LG1_RAW <- dat[dd1,]

rc <- c()
for(j in 1:(nd-1)){
  Two_res <- auto_scan(dat=dat,se=c(dd1[j],dd1[j+1]))
  rc <- rbind(rc,Two_res)
}



save(rc,file="LG1.RData")
