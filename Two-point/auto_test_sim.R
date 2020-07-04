






source("auto_est1.r")
source("auto_sim1.R")
source("auto_util1.R")



pmap <- matrix(c(1,2,1,2,
                 1,2,1,2,
                 1,1,2,2,
                 1,2,2,1),nrow=4,byrow=T)

index <- matrix(c(1,2,3,4,
                  1,2,3,4,
                  1,4,3,2,
                  1,2,3,4),nrow=4,byrow=T)

two_z <- Frz1_two_sim1(ff1=f1,n=200,m=50,pmap=pmap,index=index)

ret <- c()
for(i in 1:50){
  ret <- rbind(ret,two_rfz_phase_g2(M1=two_z$zm1[,i],M2=two_z$zm2[,i],pmap=pmap,index=index))
  cat("number=",i,"\n")
}

colMeans(ret)
apply(ret,2, sd)


two_z <- Frz1_two_sim1(ff1=f1,n=500,m=50,pmap=pmap,index=index)

ret <- c()
for(i in 1:50){
  ret <- rbind(ret,two_rfz_phase_g2(M1=two_z$zm1[,i],M2=two_z$zm2[,i],pmap=pmap,index=index))
  cat("number=",i,"\n")
}

colMeans(ret)
apply(ret,2, sd)
