
source("auto_est1.r")
source("auto_sim1.R")
source("auto_util1.R")

dat <- read.csv("example.csv")[,-1]



rc <- c()
for(j in 1:10){
  Two_res <- auto_scan_mc(dat=dat,se=c(j,j+1))
  cat("Number=",j,"\n")
  rc <- rbind(rc,Two_res)
}









