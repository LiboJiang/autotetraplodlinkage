
setwd("E:/Pro1/AutoL1/Auto_three")



source("auto_three_util.R")
source("auto_three_mai.R")
source("phase-debug.R")

#chr5
g5 <- read.csv("../data/geno5.csv")

g51 <- as.matrix(g5[,-1])

ret5 <- three_scan1(geno=g51)
  


#chr5
g9 <- read.csv("../data/geno9.csv")

g91 <- as.matrix(g9[,-1])

ret9 <- three_scan(geno=g91)

  