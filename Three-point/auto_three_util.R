



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








FG_three <- function(M1,M2){
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
      tmp1 <- sort(c(m11[i]%/%10,m11[i]%%10,m22[j]%/%10,m22[j]%%10))
      game1[i,j] <- sum(tmp1*c(1000,100,10,1))
    }
  }
  
  game2 <- matrix(NA,length(m1_g),length(m2_g))
  for(i in 1:length(m1_g)){
    for(j in 1:length(m2_g)){
      tmp2 <- sort(c( m1_g[i]%/%10, m1_g[i]%%10,m2_g[j]%/%10,m2_g[j]%%10))
      game2[i,j] <- sum(tmp2*c(1000,100,10,1))
    }
  }
  return(list(m1=m1_g,m2=m2_g,game1=game1,game2=game2))
  
}




FG1_three <-function(P1M1,P1M2,P1M3,P2M1,P2M2,P2M3){
  
  tfz1 <- FG_three(M1=P1M1,M2=P2M1)
  zoc1 <- c(t(tfz1$game2))
  
  tfz2 <- FG_three(M1=P1M2,M2=P2M2)
  zoc2 <- c(t(tfz2$game2))
  
  tfz3 <- FG_three(M1=P1M3,M2=P2M3)
  zoc3 <- c(t(tfz3$game2))
  
  zo11 <- unique(zoc1);zo22 <- unique(zoc2);zo33 <- unique(zoc3);
  zygote1 <- array(NA,dim=c(length(zo11),length(zo22),length(zo33)))
  for(i in 1:length(zo11)){
    for(j in 1:length(zo22)){
      for(k in 1:length(zo33)){
        zygote1[i,j,k] <- zo11[i]*10000*10000+zo22[j]*10000+zo33[k]
      }
    }
  }
  zoc13 <- c()
  for(k in 1:length(zoc3)){
    for(i in 1:length(zoc1)){
      zoc13 <- c(zoc13,zoc1[i]*10000+zoc3[k])
    }
  }
  return(list(zoc1=zoc1,zoc2=zoc2,zoc3=zoc3,zoc13=zoc13,zygote1=melt_m(zygote1)))
}




TM_zy_three <- function(g1,g2){
  
  g11 <- unique(g1)
  g22 <- unique(g2)
  
  n1 <- length(g11)
  n2 <- length(g22)
  
  Tr1 <- matrix(0,n1,10000)
  Tr2 <- matrix(0,n2,100)
  
  for(i in 1:n1){
    Tr1[i,which(g11[i]==g1)] <- 1
  }
  
  for(i in 1:n2){
    Tr2[i,which(g22[i]==g2)] <- 1
  }
  return(list(dL=Tr1,dR=Tr2))
}



findex <-function(){
  
  nM <- matrix(NA,59,59)
  ni <- 1:59
  for(i in 1:59){
    for(j in 1:59){
      nM[i,j] <- paste(ni[i],ni[j],sep="-")
    }
  }
  return(nM)
}




frt <- function(fii){
  
  pp <- c()
  for(i in 1:length(fii)){
    pp <- rbind(pp,strsplit(as.character(fii[i]),"-")[[1]])
  }
  
  fpm <- c()
  for(ii in 1:59){
    
    fpm <- cbind(fpm,((pp[,1]==ii)+(pp[,2]==ii))/2)
    
  }
  
  return(fpm)
  
}


melt_m <- function(Array){
  
  ii <- dim(Array)[3]
  
  rc <- c()
  for(i in 1:ii){
    rc <- rbind(rc,Array[,,i])
  }
  
  rc
}
FM_three <- function(Pf1,Pf2){
  
  
  tmp1 <- matrix(NA,10000,100)
  
  for(i in 1:100){
    for(j in 1:10){
      tmp1[((i-1)*100+1):(i*100),((j-1)*10+1):(j*10)] <- paste(Pf1[i,j],"-",Pf2,sep="")
    }
  }
  
  return(tmp1)
}


PF_zy_thr <-function(g1,g2,P1f,P2f,fii){
  
  ind <- FM_three(P1f,P2f)
  g11 <- unique(g1)
  g22 <- unique(g2)
  
  n1 <- length(g11)
  n2 <- length(g22)
  
  z1 <- matrix("",nrow=n1,ncol=100)
  
  for(i in 1:n1){
    ti1 <- which(g11[i]==g1)
    if(length(ti1)==1){
      z1[i,] <- as.character(ind[ti1,])
    }else{
      ncom <- c()
      for(j in ti1){
        ncom <- paste(ncom,ind[j,],sep="_")
      }
      z1[i,] <- ncom
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
        ncom <- paste(ncom,z1[,j],sep="_")
      }
      z1[,i] <- ncom
    }
  }
  z <- matrix(z1[1:n1,1:n2],n1,n2)
  nk <- c()
  nc <- c()
  for(i in 1:n1){
    for(j in 1:n2){
      ii <- strsplit(z[i,j],"_")[[1]]
      ii1 <- (ii[which(ii!="")])
      nc <- c(nc,ii1)
      nk <- c(nk,length(ii1))
    }
  }
  nkk <- c()
  for(i in 1:length(nk)){
    nkk <- c(nkk,rep(i,nk[i]))
  }
  ncc <- nc
  for(i in 1:length(fii)){
    ncc[which(ncc==fii[i])] <- i
    #cat(i,"\n")
  }
  trm <- matrix(0,length(fii),length(nc))
  for(i in 1:length(fii)){
    trm[i,which(nc==fii[i])] <- 1
    #cat(i,"\n")
  }
  
  return(list(ncc=ncc,nkk=nkk,trm=trm))
}








HP_thr <- function(f,P1f,P2f){
  
  t1s <- sort(unique(c(P1f)))
  
  ff1 <- f[t1s]
  ff1 <- ff1/sum(ff1)
  ind <- P1f
  st1 <- c()
  k <- 1
  for(i in t1s){
    tt <- which(ind==i)
    st1 <- c(st1,length(tt))
    ind[which(ind==i)] <- ff1[k]/length(tt)
    k <- k +1
  }
  
  ind1 <- ind
  ff11 <- f
  ff11[t1s] <- ff1/st1
  
  t2s <- sort(unique(c(P2f)))
  
  ff2 <- f[t2s]
  ff2 <- ff2/sum(ff2)
  ind <- P2f
  st2 <- c()
  k <-1
  for(i in t2s){
    
    tt <- which(ind==i)
    st2 <- c(st2,length(tt))
    ind[which(ind==i)] <- ff2[k]/length(tt)
    k <- k + 1
  }
  ff22 <- f
  ff22[t2s] <- ff2/st2
  
  f2 <- c(t(kronecker(ff11,t(ff22))))
  
  tmp2 <- kronecker(ind1,ind)
  return(list(HP=tmp2,f2=f2,fp1=ff1,fp2=ff2))
  
  
  
}


sapp_L <-function(n){
  log(sqrt(2*pi*n))+n*log(n/exp(1))
}






ff1 <- c(0,rep(1/59,59))


Frz1_three_sim1 <- function(ff1,n,m,pmap,index){
  
  
  P1f <- three_phase(M1=index[1,],M2=index[2,],M3=index[3,])
  P2f <- three_phase(M1=index[4,],M2=index[5,],M3=index[6,])
  
  
  fmz <- HP_thr(f=ff1,P1f=P1f,P2f=P2f)$HP
  
  zygotec <-  FG1_three(P1M1=pmap[1,],P1M2=pmap[2,],P1M3=pmap[3,],P2M1=pmap[4,],P2M2=pmap[5,],P2M3=pmap[6,])
  
  
  DL <- TM_zy_three(g1=zygotec$zoc13,g2=zygotec$zoc2)
  
  fmz1 <- DL$dL%*%fmz%*%t(DL$dR)
  
  zygotecv <- c(t(zygotec$zygote1))
  
  zyg_m1 <- zyg_m2 <- zyg_m3 <- matrix(NA,n,m)
  for(i in 1:m){
    tmp <- sample(zygotecv,n,replace=T,prob=c(t(fmz1)))#gamecv[xx.ord]
    tmp1 <- t(apply(as.matrix(tmp),1,function(x){as.numeric(strsplit(as.character(x),"")[[1]])}))
    zyg_m1[,i] <- tmp1[,1]*1000+tmp1[,2]*100+tmp1[,3]*10+tmp1[,4]
    zyg_m2[,i] <- tmp1[,5]*1000+tmp1[,6]*100+tmp1[,7]*10+tmp1[,8]
    zyg_m3[,i] <- tmp1[,9]*1000+tmp1[,10]*100+tmp1[,11]*10+tmp1[,12]
  }
  zo <- list(zm1=zyg_m1,zm2=zyg_m2,zm3=zyg_m3)
  return(zo)
}