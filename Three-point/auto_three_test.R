




source("auto_three_util.R")
source("auto_three_mai.R")



#index <- matrix(c(1,2,3,4,
#                  1,2,3,4,
#                  1,4,3,2,
#                  1,2,3,4,
#                  4,2,3,1,
#                  3,2,1,4),nrow=6,byrow=T)


#pmap <- matrix(c(1,1,2,2,
#                  1,1,2,2,
#                  1,1,2,2,
#                  1,1,2,2,
#                  1,1,2,2,
#                  1,1,2,2),nrow=6,byrow=T)




#M1 <- ndat[1,-c(1,2)]
#M2 <- ndat[2,-c(1,2)]
#M3 <- ndat[3,-c(1,2)]


#rett <- auto_three_scan(ndat=ndat,phase=index)
#rett

load("LG1.RData")

chr <- read.csv("Marker_info.csv",header=T)


load("F.RData")



dat <- filtered_data
dat[which(dat==0)] <- 1111
dat[which(dat==1)] <- 1112
dat[which(dat==2)] <- 1122
dat[which(dat==3)] <- 1222
dat[which(dat==4)] <- 2222



rc1 <- c()
for(i in 1:83){
  ii <- c((2*(i-1)+1):(2*i))
  cn <- unique(c(t(rc[ii,1:2])))
  ndat <- dat[match(cn,rownames(dat)),]
  
  cp1 <- c(t(rc[ii[1],11:14]))
  cp2 <- c(t(rc[ii[2],11:14]))
  
  if(cp1[2]!=cp2[1]){
    
    st12 <- as.numeric(strsplit(cp1[2],"")[[1]])
    st21 <- as.numeric(strsplit(cp2[1],"")[[1]])
    ci1 <- match(st12,st21)
    st21 <- as.numeric(strsplit(cp2[1],"")[[1]])
    
    tmp11 <- as.numeric(strsplit(cp2[1],"")[[1]])[ci1]
    tmp12 <- as.numeric(strsplit(cp2[2],"")[[1]])[ci1]
    tmp111 <- paste(tmp11[1],tmp11[2],tmp11[3],tmp11[4],sep="")
    tmp121 <- paste(tmp12[1],tmp12[2],tmp12[3],tmp12[4],sep="")
  }else{
    tmp111 <- cp2[1]
    tmp121 <- cp2[2]
  }
  
  if(cp1[4]!=cp2[3]){
    
    st14 <- as.numeric(strsplit(cp1[4],"")[[1]])
    st23 <- as.numeric(strsplit(cp2[3],"")[[1]])
    ci2 <- match(st14,st23)
    st24 <- as.numeric(strsplit(cp2[4],"")[[1]])
    
    tmp21 <- as.numeric(strsplit(cp2[3],"")[[1]])[ci2]
    tmp22 <- as.numeric(strsplit(cp2[4],"")[[1]])[ci2]
    tmp211 <- paste(tmp21[1],tmp21[2],tmp21[3],tmp21[4],sep="")
    tmp221 <- paste(tmp22[1],tmp22[2],tmp22[3],tmp22[4],sep="")
  }else{
    tmp211 <- cp2[3]
    tmp221 <- cp2[4]
  }
  ncp2 <- c(tmp111,tmp121,tmp211,tmp221)
  
  
  allp <- c(cp1[1:2],ncp2[2],cp1[3:4],ncp2[4])
  
  Phase <- rbind(as.numeric(strsplit(allp[1],"")[[1]]),
                 as.numeric(strsplit(allp[2],"")[[1]]),
                 as.numeric(strsplit(allp[3],"")[[1]]),
                 as.numeric(strsplit(allp[4],"")[[1]]),
                 as.numeric(strsplit(allp[5],"")[[1]]),
                 as.numeric(strsplit(allp[6],"")[[1]]))
  
  Three_res <- auto_three_scan(ndat=ndat,phase=Phase)
  #Two_res <- auto_scan(dat=dat,)
  rc1 <- rbind(rc1,Three_res)
}

save(rc1,file="Three_LG1.RData")


#allp <- c("1234","4321","2341","1234","1234","1234")

