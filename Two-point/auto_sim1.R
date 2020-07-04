


f1 <- c(0.0401535,0.0012,0.0082,0.0004465,0.0532,0.0054465,0.88,0.001185908,0.01016759)




Frz1_two_sim1 <-function(ff1,n,m,pmap,index){
  
  
  
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



FG <- function(M1,M2){
  M1 <- M1
  m1_l <- cbind(matrix(rep(M1,2),nrow=2,byrow =T),matrix(M1[combn(4,2)],2))
  m1_l1 <- m1_l
  m1_g <- m1_l1[1,]*10+m1_l1[2,]
  M2 <- M2
  m2_l <- cbind(matrix(rep(M2,2),nrow=2,byrow =T),matrix(M2[combn(4,2)],2))
  m2_l1 <- m2_l
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


FG1 <-function(P1M1,P1M2,P2M1,P2M2){
  
  #if(all(P1M1==1)){
  #  tfz1 <- FG(M1=P2M1,M2=P1M1)
  #}else{
    tfz1 <- FG(M1=P1M1,M2=P2M1)
  #}
  
  zoc1 <- c(t(tfz1$game2))
  
  #if(all(P1M2==1)){
  #  tfz2 <- FG(M1=P2M2,M2=P1M2)
  #}else{
    tfz2 <- FG(M1=P1M2,M2=P2M2)
  #}
  
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




HP2 <- function(f,P1f,P2f){
  
  ind <- P1f
  t1 <- which(ind==1)
  ind[which(ind==1)] <- f[1]/length(t1)
  t2 <- which(ind==2)
  ind[which(ind==2)] <- f[2]/length(t2)
  t3 <- which(ind==3)
  ind[which(ind==3)] <- f[3]/length(t3)
  t4 <- which(ind==4)
  ind[which(ind==4)] <- f[4]/length(t4)
  t5 <- which(ind==5)
  ind[which(ind==5)] <- f[5]/length(t5)
  t6 <- which(ind==6)
  ind[which(ind==6)] <- f[6]/length(t6)
  t7 <- which(ind==7)
  ind[which(ind==7)] <- f[7]/length(t7)
  t8 <- which(ind==8)
  ind[which(ind==8)] <- f[8]/length(t8)
  t9 <- which(ind==9)
  ind[which(ind==9)] <- f[9]/length(t9)
  
  ind1 <- ind
  ff1 <- f/c(length(t1),length(t2),length(t3),length(t4),length(t5),length(t6),length(t7),length(t8),length(t9))
  
  ind <- P2f
  t1 <- which(ind==1)
  ind[which(ind==1)] <- f[1]/length(t1)
  t2 <- which(ind==2)
  ind[which(ind==2)] <- f[2]/length(t2)
  t3 <- which(ind==3)
  ind[which(ind==3)] <- f[3]/length(t3)
  t4 <- which(ind==4)
  ind[which(ind==4)] <- f[4]/length(t4)
  t5 <- which(ind==5)
  ind[which(ind==5)] <- f[5]/length(t5)
  t6 <- which(ind==6)
  ind[which(ind==6)] <- f[6]/length(t6)
  t7 <- which(ind==7)
  ind[which(ind==7)] <- f[7]/length(t7)
  t8 <- which(ind==8)
  ind[which(ind==8)] <- f[8]/length(t8)
  t9 <- which(ind==9)
  ind[which(ind==9)] <- f[9]/length(t9)
  
  ff2 <- f/c(length(t1),length(t2),length(t3),length(t4),length(t5),length(t6),length(t7),length(t8),length(t9))
  
  f2 <- c(t(kronecker(ff1,t(ff2))))
  
  tmp2 <- kronecker(ind1,ind)
  return(list(HP=tmp2,f2=f2,fp1=ff1,fp2=ff2))
  
  
  
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