

setwd("E:/Pro1/AutoL1/Auto_two")


source("auto_sim.R")
source("auto_est.R")
source("auto-debug.R")

#Simulation for a=0.05,b=0.1,r=0.05
#1111=0, 1112=1, 1122=2, 1222=3, 2222=4
# Zygote model exclude HO*HO 
#1x0| 2x0| 0x1| 1x1| 2x1| 0x2| 1x2| 2x2| 1x3|
pmap <- matrix(c(2,1,2,1,
                 1,2,1,2,
                 1,2,1,2,
                 1,2,1,2),nrow=4,byrow=T)

index <- matrix(c(1,2,3,4,
                  3,2,1,4,
                  1,2,3,4,
                  3,2,1,4),nrow=4,byrow=T)

two_z <- Frz1_two_sim1(ff1=f1,n=500,m=50,pmap=pmap,index=index)

index1 <- matrix(c(1,2,3,4,
                   3,2,1,4,
                   3,2,1,4,
                   1,2,3,4),nrow=4,byrow=T)
i <- 3

two_rfz_phase_g2(M1=two_z$zm1[,i],M2=two_z$zm2[,i],pmap=pmap,index=index1)

res <- c()
for(i in 1:50){
  tmp <- two_rfz_phase_g2(M1=two_z$zm1[,i],M2=two_z$zm2[,i],pmap=pmap,index=index1)
  cat("i= ",i,"\n")
  res <- rbind(res,tmp)
}

colMeans(res)

apply(res,2,sd)
