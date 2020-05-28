

###Game
#a=0.05 b=0.1 r=0.05

f1 <- c(0.0401535,0.0012,0.0082,0.0004465,0.0532,0.0054465,0.88,0.001185908,0.01016759)


f2 <- c(0.04071,0.00130,0.0082,0.00446,0.00353,0.0431,0.88221,0.00736,0.00245)


FM <- function(){matrix(c(1,2,2,2,5,5,5,6,6,6,
               2,1,2,2,5,6,6,5,5,6,
               2,2,1,2,6,5,6,5,6,5,
               2,2,2,1,6,6,5,6,5,5,
               3,3,4,4,7,8,8,8,8,9,
               3,4,3,4,8,7,8,8,9,8,
               3,4,4,3,8,8,7,9,8,8,
               4,3,3,4,8,8,9,7,8,8,
               4,3,4,3,8,9,8,8,7,8,
               4,4,3,3,9,8,8,8,8,7),nrow=10,ncol=10,byrow=T)}



FG <- function(M1,M2){
  M1 <- sort(M1)
  m1_l <- cbind(matrix(rep(M1,2),nrow=2,byrow =T),matrix(M1[combn(4,2)],2))
  m1_l1 <- sort_g(m1_l)
  m1_g <- m1_l1[1,]*10+m1_l1[2,]
  M2 <- sort(M2)
  m2_l <- cbind(matrix(rep(M2,2),nrow=2,byrow =T),matrix(M2[combn(4,2)],2))
  m2_l1 <- sort_g(m2_l)
  m2_g <- m2_l1[1,]*10+m2_l1[2,]

  
  #game <- matrix(rep(m1_g,10),10,10,byrow = F)+matrix(rep(m2_g,10)*100,10,10,byrow = T)
  
  m11 <- unique(m1_g);m22 <- unique(m2_g)
  game1 <- matrix(NA,length(m11),length(m22))
  for(i in 1:length(m11)){
    for(j in 1:length(m22)){
      game1[i,j] <- m11[i]*100+m22[j]
    }
  }
  
  game2 <- matrix(NA,length(m1_g),length(m2_g))
  for(i in 1:length(m1_g)){
    for(j in 1:length(m2_g)){
      game2[i,j] <- m1_g[i]*100+m2_g[j]
    }
  }
  return(list(m1=m1_g,m2=m2_g,game1=game1,game2=game2))
  
}





Fr1_two_sim <- function(f,n,m,M1=c(1,2,3,4),M2=c(1,2,3,4),index){
  
  P1f <- two_phase(M1=index[1,],M2=index[2,])
  fm <- HPP(f=f,P1f)
  
  
  
  tfg <- FG(M1=M1,M2=M2)
  gamec <- tfg$game1
  
  DL <- TM(g1=tfg$m1,g2=tfg$m2)
  
  fm1 <- DL$dL%*%fm%*%t(DL$dR)
  
  gamecv <- c(t(gamec))
  
  gam_m <- matrix(NA,n,m)
  for(i in 1:m){
    gam_m[,i] <- sample(gamecv,n,replace=T,prob=c(t(fm1)))#gamecv[xx.ord]
  }
  
  
  return(gam_m)
}


sort_g <-function(gg){
  
  pn <- dim(gg)[2]
  
  for (i in 1:pn){
    t1 <- gg[1,i]
    t2 <- gg[2,i]
    if(t2<t1){
      gg[1,i] <- t2
      gg[2,i] <- t1
    }
  }
  return(gg)
}


####ZYGOTE


FG1 <- function(P1M1,P1M2,P2M1,P2M2){
  
  if(all(P1M1==1)){
    tfz1 <- FG(M1=P2M1,M2=P1M1)
  }else{
    tfz1 <- FG(M1=P1M1,M2=P2M1)
  }
  
  zoc1 <- c(t(tfz1$game2))
  
  if(all(P1M1==1)){
    tfz1 <- FG(M1=P2M2,M2=P1M2)
  }else{
    tfz2 <- FG(M1=P1M2,M2=P2M2)
  }
  
  zoc2 <- c(t(tfz2$game2))
  
  zo11 <- unique(zoc1);zo22 <- unique(zoc2)
  zygote1 <- matrix(NA,length(zo11),length(zo22))
  for(i in 1:length(zo11)){
    for(j in 1:length(zo22)){
      zygote1[i,j] <- zo11[i]*10000+zo22[j]
    }
  }
  return(list(zoc1=zoc1,zoc2=zoc2,zygote1=zygote1))
}


TM_zy <- function(g1,g2){
  
  g11 <- unique(g1)
  g22 <- unique(g2)
  
  n1 <- length(g11)
  n2 <- length(g22)
  
  Tr1 <- matrix(0,n1,100)
  Tr2 <- matrix(0,n2,100)
  
  for(i in 1:n1){
    Tr1[i,which(g11[i]==g1)] <- 1
  }
  
  for(i in 1:n2){
    Tr2[i,which(g22[i]==g2)] <- 1
  }
  return(list(dL=Tr1,dR=Tr2))
}







Frz1_two_sim <- function(ff1,n,m,P1M1=P1M1,P1M2=P1M2,P2M1=P2M1,P2M2=P2M2){
  
  
  
  fm <- HP(f=ff1)
  fmz <- kronecker(fm,fm)
  
  zygotec <- FG1(P1M1=P1M1,P1M2=P1M2,P2M1=P2M1,P2M2=P2M2)
  
  
  DL <- TM_zy(g1=zygotec$zoc1,g2=zygotec$zoc2)
  
  fmz1 <- DL$dL%*%fmz%*%t(DL$dR)
  
  zygotecv <- c(t(zygotec$zygote1))
  
  zyg_m1 <- zyg_m2 <- matrix(NA,n,m)
  for(i in 1:m){
    tmp <- sample(zygotecv,n,replace=T,prob=c(t(fmz1)))#gamecv[xx.ord]
    zyg_m1[,i] <- signif(tmp/10000,4)
    zyg_m2[,i] <- tmp%%10000
  }
  zo <- list(zm1=zyg_m1,zm2=zyg_m2)
  return(zo)
}

Frz1_two_sim1 <- function(ff1,n,m,pmap,index){
  
  
  
  P1f <- two_phase(M1=index[1,],M2=index[2,])
  P2f <- two_phase(M1=index[3,],M2=index[4,])
  
  fmz <- HP2(f=ff1,P1f=P1f,P2f=P2f)$HP
  #fmz <- kronecker(fm$HP,fm$HP)
  
  zygotec <-  FG1(P1M1=pmap[1,],P1M2=pmap[2,],P2M1=pmap[3,],P2M2=pmap[4,])
  
  
  DL <- TM_zy(g1=zygotec$zoc1,g2=zygotec$zoc2)
  
  fmz1 <- DL$dL%*%fmz%*%t(DL$dR)
  
  zygotecv <- c(t(zygotec$zygote1))
  
  zyg_m1 <- zyg_m2 <- matrix(NA,n,m)
  for(i in 1:m){
    tmp <- sample(zygotecv,n,replace=T,prob=c(t(fmz1)))#gamecv[xx.ord]
    zyg_m1[,i] <- signif(tmp/10000,4)
    zyg_m2[,i] <- tmp%%10000
  }
  zo <- list(zm1=zyg_m1,zm2=zyg_m2)
  return(zo)
}



PFF <- function(g1,g2,P1f){
  
  ind <- P1f
  g11 <- unique(g1)
  g22 <- unique(g2)
  
  n1 <- length(g11)
  n2 <- length(g22)
  
  z1 <- matrix("",nrow=n1,ncol=10)
  
  for(i in 1:n1){
    ti1 <- which(g11[i]==g1)
    if(length(ti1)==1){
      z1[i,] <- as.character(ind[ti1,])
    }else{
      ncom <- c()
      for(j in ti1){
        ncom <- paste(ncom,ind[j,],sep="")
        z1[i,] <- ncom
      }
    }
  }
  
  #z <- matrix("",nrow=n1,ncol=n2)
  for(i in 1:n2){
    ti2 <- which(g22[i]==g2)
    if(length(ti2)==1){
      z1[,i] <- as.character(z1[,ti2])
    }else{
      ncom <- c()
      for(j in ti2){
        ncom <- paste(ncom,z1[,j],sep="")
        z1[,i] <- ncom
      }
    }
  }
  z <- matrix(z1[1:n1,1:n2],n1,n2)
  return(z)
}




HPP <-function(f,P1f){
  
  f1 <- f[1]/4;f2 <- f[2]/12;f3 <- f[3]/12;f4 <- f[4]/12;f5 <- f[5]/12;
  f6 <- f[6]/12;f7 <- f[7]/6;f8 <- f[8]/24;f9 <- f[9]/6
  
  ind <- P1f
  ind[which(ind==1)] <- f1
  ind[which(ind==2)] <- f2
  ind[which(ind==3)] <- f3
  ind[which(ind==4)] <- f4
  ind[which(ind==5)] <- f5
  ind[which(ind==6)] <- f6
  ind[which(ind==7)] <- f7
  ind[which(ind==8)] <- f8
  ind[which(ind==9)] <- f9
  
  
  return(ind)
  
  
  
}
