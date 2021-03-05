# Autotetraploidlinkage
A General Framework for Statistical Linkage Analysis in Autotetraplod

1.	Introduction
At present, this package is used to linkage analysis in the full-sib of autotetraploid based on the two-point and three-point model. This guide gives some brief instructions on how to perform the tasks of linkage analysis by this package. The outline of this guide is as follows: 
2.	Data format

ID                        P1   P2  IND_1  IND_10  IND_100  IND_101

solcap_snp_c2_23780       2    0     9       0       1        1

solcap_snp_c2_23781       0    2     2       1       1        1

solcap_snp_c2_23803       1    1     2       1       2        1

solcap_snp_c2_23804       0    2     1       1       1        1

solcap_snp_c2_238045      0    2     1       1       1        1

solcap_snp_c2_7632        0    2     9       1       1        1

solcap_snp_c2_23643       1    0     0       1       0        1

solcap_snp_c2_23669       3    2     9       2       2        3

solcap_snp_c2_23678       1    2     2       2       2        1

solcap_snp_c2_23717       1    1     2       1       2        1

Five genotypes (aaaa=0, Aaaa=1,AAaa=2, AAAa=3, AAAA=4) and missing data (coded as 9 ) are valid marker values. 

3.	Computer simulation

#Two-point model

#load functions

source("auto_sim.R") 

source("auto_est.R")

source("auto-debug.R")


#Simulation for a=0.05,b=0.1,r=0.05

#1111=0, 1112=1, 1122=2, 1222=3, 2222=4

#set the marker type

pmap <- matrix(c(2,1,2,1,1,2,1,2,1,2,1,2,1,2,1,2),nrow=4,byrow=T)
                 
#set the linkage phase

index <- matrix(c(1,2,3,4，3,2,1,4,1,2,3,4,3,2,1,4),nrow=4,byrow=T)
                  
#data simulation

two_z <- Frz1_two_sim1(ff1=f1,n=500,m=50,pmap=pmap,index=index)

#f1 indicates the f frequency, n indicates the sample size, m indicates the number of markers.

#estimate recombination fraction and DR

ret_two <-two_rfz_phase_g2(M1=two_z$zm1[,i],M2=two_z$zm2[,i],pmap=pmap,index=index1)

#Three-point model

#load functions

source("auto_three_util.R")

source("auto_three_mai.R")

source("phase-debug.R")

#set the linkage phase

index <- matrix(c(1,2,3,4,1,3,2,4,1,2,3,4,4,3,2,1,1,2,3,4,
               1,4,2,3),nrow=6,byrow=T)

#set the marker type
pmap <- matrix(c(1,1,2,2,
              1,1,2,2,
              1,1,2,2,
              1,1,2,2,
              1,1,2,2,
              1,1,2,2),nrow=6,byrow=T)

#data simulation by three-point model
dat <- Frz1_three_sim1(ff1=ff1,n=500,m=10,pmap=pmap,index=index)
#ff1 indicates the g frequency
# estimate recombination fraction, cocand DR
ret_three<- three_rfz_phase(M1=dat$zm1[,i],M2=dat$zm2[,i],M3=dat$zm3[,i],
pmap=pmap,index=index)

2.	Work example
#chr5 and chr9
#read genotype file
g5 <- read.csv("../data/geno5.csv")
g51 <- as.matrix(g5[,-1])
g9 <- read.csv("../data/geno9.csv")
g91 <- as.matrix(g9[,-1])

#two-point analysis for the chromosome 5 and 9
ex_chr5 <- work_test(geno=g51)
ex_chr9 <- work_test(geno=g91)

#three-point analysis for the chromosome 5 and 9
three_chr5 <- three_scan1(geno=g51)
three_chr9 <- three_scan1(geno=g91)
