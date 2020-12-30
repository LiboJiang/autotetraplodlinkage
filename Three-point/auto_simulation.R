
setwd("E:/Pro1/AutoL1/simulation_R/")



source("../Auto_three/auto_three_util.R")
source("../Auto_three/auto_three_mai.R")
source("../Auto_three/phase-debug.R")

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

load("smf.RData")

apply(smf, 1, fr)
ret11 <- list()
for(i in 1:100){
  dat1 <- Frz1_three_sim1(ff1=smf[1,],n=500,m=1,pmap=pmap,index=index)
  dat2 <- Frz1_three_sim1(ff1=smf[2,],n=500,m=1,pmap=pmap,index=index)
  dat3 <- Frz1_three_sim1(ff1=smf[3,],n=500,m=1,pmap=pmap,index=index)
  dat4 <- Frz1_three_sim1(ff1=smf[4,],n=500,m=1,pmap=pmap,index=index)
  dat5 <- Frz1_three_sim1(ff1=smf[5,],n=500,m=1,pmap=pmap,index=index)
  
  
  
  
  ret1 <- three_rfz_phase(M1=dat1$zm1[,1],M2=dat1$zm2[,1],M3=dat1$zm3[,1],pmap=pmap,index=index)
  ret1
  ret2 <- three_rfz_phase(M1=dat2$zm1[,1],M2=dat2$zm2[,1],M3=dat2$zm3[,1],pmap=pmap,index=index)
  ret2
  ret3 <- three_rfz_phase(M1=dat3$zm1[,1],M2=dat3$zm2[,1],M3=dat3$zm3[,1],pmap=pmap,index=index)
  ret3
  ret4 <- three_rfz_phase(M1=dat4$zm1[,1],M2=dat4$zm2[,1],M3=dat4$zm3[,1],pmap=pmap,index=index)
  ret4
  ret5 <- three_rfz_phase(M1=dat5$zm1[,1],M2=dat5$zm2[,1],M3=dat5$zm3[,1],pmap=pmap,index=index)
  ret5
  ret11[[i]] <- rbind(ret1$res1,ret2$res1,ret3$res1,ret4$res1,ret5$res1)
}




