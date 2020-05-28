




source("auto_two_sim.R")
source("auto_two_est.R")
source("auto_two_debug.R")

#1
#Simulation for a=0.05,b=0.1,r=0.05
#1111=0, 1112=1, 1122=2, 1222=3, 2222=4
# Game model exclude 0 and 4 
#M1 <- c(1,2,3,4)
#M2 <- c(1,2,3,4)
#index <- matrix(c(1,2,3,4,
#                  1,2,4,3),nrow=2,byrow=T)
#two_g <- Fr1_two_sim(f=f1,n=500,m=20,M1=M1,M2=M2,index=index)


#index1 <- matrix(c(1,2,3,4,
#                  1,2,4,3),nrow=2,byrow=T)
#pa_two_g <- c()
#for( i in 1:20){
#  pa_two_g <- rbind(pa_two_g,two_rf1(M=two_g[,i],M1=M1,M2=M2,index1=index1))
#}
#colMeans(pa_two_g)
#apply(pa_two_g, 2,sd)


#Simulation for a=0.05,b=0.1,r=0.05
#1111=0, 1112=1, 1122=2, 1222=3, 2222=4
# Zygote model exclude HO*HO 
#1x0| 2x0| 0x1| 1x1| 2x1| 0x2| 1x2| 2x2| 1x3|
pmap <- matrix(c(1,1,1,2,
                 1,1,1,2,
                 1,1,2,2,
                 1,1,2,2),nrow=4,byrow=T)

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

#two_rfz1(M1=two_z$zm1[,i],M2=two_z$zm2[,i],pmap=pmap,index=index1)

two_rfz_phase_g2(M1=two_z$zm1[,i],M2=two_z$zm2[,i],pmap=pmap,index=index1)



ret <- two_full_cal(M1=two_z$zm1[,i],M2=two_z$zm2[,i],pm=pmap)




