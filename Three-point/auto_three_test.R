
setwd("E:/Pro1/AutoL1/Auto_three")



source("auto_three_util.R")
source("auto_three_mai.R")
source("phase-debug.R")

index <- matrix(c(1,2,3,4,
                  1,3,2,4,
                  1,2,3,4,
                  4,3,2,1,
                  1,2,3,4,
                  1,4,2,3),nrow=6,byrow=T)


pmap <- matrix(c(1,1,2,2,
                  1,1,2,2,
                  1,1,2,2,
                  1,1,2,2,
                  1,1,2,2,
                  1,1,2,2),nrow=6,byrow=T)


dat <- Frz1_three_sim1(ff1=ff1,n=500,m=10,pmap=pmap,index=index)


i <- 1
ret <- three_rfz_phase(M1=dat$zm1[,i],M2=dat$zm2[,i],M3=dat$zm3[,i],pmap=pmap,index=index)
ret

