



setwd("E:/Pro1/AutoL1/Auto_two")


source("auto_sim.R")
source("auto_est.R")
source("auto-debug.R")



#chr5
g5 <- read.csv("../data/geno5.csv")

g51 <- as.matrix(g5[,-1])


#ex_test5 <- work_test(geno=g51[1:3,])
ex_chr5 <- work_test(geno=g51)





#chr9
g9 <- read.csv("../data/geno9.csv")

g91 <- as.matrix(g9[,-1])


#ex_test9 <- work_test(geno=g91[1:3,])
ex_chr9 <- work_test(geno=g91)



#library(onemap)

write.csv(ex_chr5$para,file="ch5_Para.csv")
write.csv(ex_chr5$ph1,file="ch5_ph1.csv")
write.csv(ex_chr5$ph2,file="ch5_ph2.csv")


write.csv(ex_chr9$para,file="ch9_Para.csv")
write.csv(ex_chr9$ph1,file="ch9_ph1.csv")
write.csv(ex_chr9$ph2,file="ch9_ph2.csv")