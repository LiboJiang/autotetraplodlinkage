

two_phase1 <- function(M1,M2){
  
  
  mc <- rbind(M1,M2)
  #mc <- sort_g(mc)
  mc1 <- c(mc[1,]*10+mc[2,])
  
  M11 <- sort(M1)
  m1_l <- cbind(matrix(rep(M11,2),nrow=2,byrow =T),matrix(M11[combn(4,2)],2))
  m1_l1 <- sort_g(m1_l)
  m1_g <- m1_l1[1,]*10+m1_l1[2,]
  
  M22 <- sort(M2)
  m2_l <- cbind(matrix(rep(M22,2),nrow=2,byrow =T),matrix(M22[combn(4,2)],2))
  m2_l1 <- sort_g(m2_l)
  m2_g <- m2_l1[1,]*10+m2_l1[2,]
  
  #m12_1 <- m1_l1[1,]*10+ m2_l1[1,]
  #m12_2 <- m1_l1[2,]*10+ m2_l1[2,]
  
  nm <- matrix(0,nrow=length(m1_g),ncol=length(m2_g))
  for(i in 1:length(m1_g)){
    for(j in 1:length(m2_g)){
      tmp1 <- c(m1_l1[1,i],m2_l1[1,j])
      tmp2 <- c(m1_l1[2,i],m2_l1[2,j])
      kl1 <- m1_l1[1,i]*10+m2_l1[1,j]
      kl2 <- m1_l1[2,i]*10+m2_l1[2,j]
      kl3 <- signif(kl1/10,1)*10+kl2%%10
      kl4 <- signif(kl2/10,1)*10+kl1%%10
      #1
      if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(length(which(kl1==mc1))>0)&&(length(which(kl2==mc1))>0)){
        nm[j,i]=1#nm[j,i]+1
      }
      #2
      if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(length(which(kl1==mc1))==0)&&(length(which(kl1==mc1))==0)){
        nm[j,i]=2#nm[j,i]+1
      }
      #3
      if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(((length(which(kl1==mc1))==0)&&(length(which(kl2==mc1))>0))
                                                              ||(((length(which(kl1==mc1))>0)&&(length(which(kl2==mc1))==0))))){
        nm[j,i]=3#nm[j,i]+1
      }
      #4
      if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(length(which(kl1==mc1))==0)&&(length(which(kl2==mc1))==0)){
        nm[j,i]=4#nm[j,i]+1
      }
      #5
      if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(((length(which(kl1==mc1))==0)&&(length(which(kl2==mc1))>0))
                                                              ||(((length(which(kl1==mc1))>0)&&any(length(which(kl2==mc1))==0))))){
        nm[j,i]=5#nm[j,i]+1
      }
      #6
      if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(length(which(kl1==mc1))==0)&&(length(which(kl2==mc1))==0)){
        nm[j,i]=6#nm[j,i]+1
      }
      #7
      if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&((((length(which(kl1==mc1))>0)&&(length(which(kl2==mc1))>0))&&
                                                               (length(which(kl3==mc1))==0)&&(length(which(kl4==mc1))==0))||(((length(which(kl1==mc1))==0)&&(length(which(kl2==mc1))==0))&&
                                                                                                                             (length(which(kl3==mc1))>0)&&(length(which(kl4==mc1))>0)))){
        nm[j,i]=7#nm[j,i]+1
      }
      #8
      if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(((((length(which(kl1==mc1))==0)&&(length(which(kl2==mc1))>0))||((length(which(kl1==mc1))>0)&&(length(which(kl2==mc1))==0)))&&
                                                               ((length(which(kl3==mc1))==0)&&(length(which(kl4==mc1))==0)))||((((length(which(kl3==mc1))==0)&&(length(which(kl4==mc1))>0))||
                                                                                                                                ((length(which(kl3==mc1))>0)&&(length(which(kl4==mc1))==0)))&&((length(which(kl1==mc1))==0)&&(length(which(kl2==mc1))==0)))))
      {
        nm[j,i]=8#nm[j,i]+1
      }
      #9
      if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&((length(which(kl1==mc1))==0)&&(length(which(kl2==mc1))==0))&&
         ((length(which(kl3==mc1))==0)&&(length(which(kl4==mc1))==0))){
        nm[j,i]= 9#nm[j,i]+1
      }
    }
  }
  return(nm)
}




two_phase <- function(M1,M2){
  
  mc <- expand.grid(1:4,1:4,1:4,1:4)
  
  mc1 <- apply(mc,1,function(x){
    if(max(table(x))==1)
      return(1)
    else{
      return(0)
    }
  })
  
  mc2 <- as.matrix(mc[as.logical(mc1),])
  
  
  mc <- rbind(M1,M2)
  #mc <- sort_g(mc)
  mc12 <- c(mc[1,]*10+mc[2,])
  
  
  
  M11 <- sort(M1)
  m1_l <- cbind(matrix(rep(M11,2),nrow=2,byrow =T),matrix(M11[combn(4,2)],2))
  m1_l1 <- sort_g(m1_l)
  m1_g <- m1_l1[1,]*10+m1_l1[2,]
  
  M22 <- sort(M2)
  m2_l <- cbind(matrix(rep(M22,2),nrow=2,byrow =T),matrix(M22[combn(4,2)],2))
  m2_l1 <- sort_g(m2_l)
  m2_g <- m2_l1[1,]*10+m2_l1[2,]
  
  nm <- matrix(0,nrow=length(m1_g),ncol=length(m2_g))
  
  for(ii in 1:dim(mc2)[1]){
    mh1 <- M1[mc2[ii,]]
    mh2 <- M2[mc2[ii,]]
    
    for(i in 1:length(m1_g)){
      for(j in 1:length(m2_g)){
        
        kl112 <- m1_l1[1,i]*10+m2_l1[1,j]
        kl212 <- m1_l1[2,i]*10+m2_l1[2,j]
        
        if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(is.element(kl112,mc12))&&
           (is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1]))){
          nm[j,i]=1
        }
        
        if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(!is.element(kl112,mc12))&&
           (!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2]))){
          nm[j,i]=2
        }
        
        if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&
           (((is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])))||
            ((!is.element(kl112,mc12))&&(is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1]))))){
          nm[j,i]=3
        }
        if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&
           (((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])))||
            ((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2]))))){
          nm[j,i]=4
        }
        if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&
           (((is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])))||
            ((!is.element(kl112,mc12))&&(is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1]))))){
          nm[j,i]=5
        }
        if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&
           (((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])))||
            ((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1]))))){
          nm[j,i]=6
        }
        if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&
           (((is.element(kl112,mc12))&&(is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])))||
            ((is.element(kl112,mc12))&&(is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1]))))){
          nm[j,i]=7
        }
        if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&
           (((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])))||
            ((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2]))))){
          nm[j,i]=7
        }
        if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&
           (((is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])))||
            ((!is.element(kl112,mc12))&&(is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1]))))){
          nm[j,i]=8
        }
        if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&
           (((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])))||
            ((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3]))))){
          nm[j,i]=8
        }
        if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&
           (((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])))||
            ((!is.element(kl112,mc12))&&(!is.element(kl212,mc12))&&((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3]))))){
          nm[j,i]=9
        }
      }
    }
  }
  
  
  
  
  
  return(nm)
}



HP2<-function(f,P1f,P2f){
  
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



FM111 <- function(Pf1,Pf2){
  
  
  tmp1 <- matrix(NA,100,100)
  
  for(i in 1:10){
    for(j in 1:10){
      tmp1[((i-1)*10+1):(i*10),((j-1)*10+1):(j*10)] <- Pf1[i,j]*10+Pf2
    }
  }
  
  return(tmp1)
}



PF_zy1 <- function(g1,g2,P1f,P2f){
  
  ind <- FM111(P1f,P2f)
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
        ncom <- paste(ncom,ind[j,],sep="-")
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
        ncom <- paste(ncom,z1[,j],sep="-")
        z1[,i] <- ncom
      }
    }
  }
  z <- matrix(z1[1:n1,1:n2],n1,n2)
  return(z)
}




PF_zy1_1 <- function(g1,g2,P1f,P2f,fii){
  
  ind <- FM111(P1f,P2f)
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
        ncom <- paste(ncom,ind[j,],sep="-")
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
        ncom <- paste(ncom,z1[,j],sep="-")
        z1[,i] <- ncom
      }
    }
  }
  z <- matrix(z1[1:n1,1:n2],n1,n2)
  nk <- c()
  nc <- c()
  for(i in 1:n1){
    for(j in 1:n2){
      ii <- strsplit(z[i,j],"-")[[1]]
      ii1 <- as.numeric(ii[which(ii!="")])
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
  }
  trm <- matrix(0,length(fii),length(nc))
  for(i in 1:length(fii)){
    trm[i,which(nc==fii[i])] <- 1
  }

  return(list(ncc=ncc,nkk=nkk,trm=trm))
}


two_rfz_phase_g2 <- function(M1,M2,pmap,index){
    
  P1f <- two_phase(M1=index[1,],M2=index[2,])
  P2f <- two_phase(M1=index[3,],M2=index[4,])
  
  M <- M1*10000+M2
  GMZ <- FG1(P1M1=pmap[1,],P1M2=pmap[2,],P2M1=pmap[3,],P2M2=pmap[4,])
  #GMZ <- FG1(P1M1=P1M1,P1M2=P1M2,P2M1=P2M1,P2M2=P2M1)
  DLRZ <- TM_zy(g1=GMZ$zoc1,g2=GMZ$zoc2)
  
  
  
  n1 <- length(unique(GMZ$zoc1))
  n2 <- length(unique(GMZ$zoc2))
  count <- matrix(0,n1,n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      count[i,j] <- length(which(GMZ$zygote1[i,j]==M))
    }
  }
  #cat("Count=",sum(count),"\n")
  #if(sum(count)==0){
  #  return(c(NA,NA,NA,NA))
  #}
  #tf <- sample(1:9,100,replace = T)
  #f_1 <- as.numeric(table(tf)/sum(tf))
  f_1 <- c(0.0401535,0.0012,0.0082,0.0004465,0.0532,0.0054465,0.88,0.001185908,0.01016759)
  
  fii <- c(t(findex()))
  frt1 <- frt(fii)
  
  PFCZ1 <- PF_zy1_1(g1=GMZ$zoc1,g2=GMZ$zoc2,P1f,P2f,fii)
  
  iter <- 1
  while(1){
    
    ft <- HP2(f=f_1,P1f=P1f,P2f=P2f)
    fc <- DLRZ$dL%*%ft$HP%*%t(DLRZ$dR)
    f_2 <- ft$f2
    
    EM1 <- f_2[PFCZ1$ncc]*c(t(count))[PFCZ1$nkk]/c(t(fc))[PFCZ1$nkk]

    
    EM2 <- PFCZ1$trm%*%as.matrix(EM1)
    
    E1 <- t(EM2)%*%frt1
    
    f_3 <- E1/sum(count)
    
    #E1 <- D(r=r,g1=GM$m1,g2=GM$m2)
    
    #r1 <- sum(E1*count)/(2*sum(count))
    
    if(all(abs(f_3-f_1)<1e-5)||iter>130)
      break
    
    f_1 <- f_3
    #cat("iter = "," ",iter,"",f_1,"\n")
    iter <- iter +1
  }
  

  iter <- 1
  r <- 0.5
  while(1){
    v1 <- r^2/(9-18*r+10*r^2)
    v2 <- r/(3-2*r)
    
    r1 <- (f_3[3]+f_3[5]+2*(f_3[2]+f_3[4]+f_3[6]+f_3[9])+2*v1*f_3[7]+(1+v2)*f_3[8])/2
    
    #r1 <- sum(E1*count)/(2*sum(count))
    
    if(abs(r1-r)<1e-5||iter>50)
      break
    
    r <- r1
    iter <- iter +1
    
  }
    
    a <- sum(f_3[1:4])
    b <- sum(f_3[c(1,2,5,6)])
    ii <- f_3!=0
    #
    LL <- sum(E1[ii]*log(f_3[ii]),nan.rm=T)
    res <- c(a,b,r,LL)
    res
  }



sapp_L <- function(n){
  log(sqrt(2*pi*n))+n*log(n/exp(1))
}




work_test <- function(geno){
  
  n <- dim(geno)[1]
  geno[which(geno==0)] <- 1111
  geno[which(geno==1)] <- 1112
  geno[which(geno==2)] <- 1122
  geno[which(geno==3)] <- 1222
  geno[which(geno==4)] <- 2222
  #1111=0, 1112=1, 1122=2, 1222=3, 2222=4
  
  
  
  ep <- as.matrix(expand.grid(1:4,1:4,1:4,1:4))
  ep1 <- ep[which(apply(ep,1,function(x){length(unique(x))})==4),]
  phf1 <- phf2 <- matrix(NA,nrow=n,ncol=4)
  
  phf1[1,] <- c(1,2,3,4)
  phf2[1,] <- c(1,2,3,4)
  abr <- c()
  #(n-1)
  for(i in 1:(n-1)){
    
    index1 <- as.numeric(phf1[i,])
    index3 <- as.numeric(phf2[i,])
    
    m1 <- as.numeric(geno[i,])
    m2 <- as.numeric(geno[i+1,])
    m11 <- m1[-c(1:2)]
    m22 <- m2[-c(1:2)]
    
    p1_1 <-  sort(as.numeric(strsplit(as.character(m1[1]),"")[[1]]))
    p1_2 <-  sort(as.numeric(strsplit(as.character(m2[1]),"")[[1]]))
    p2_1 <-  sort(as.numeric(strsplit(as.character(m1[2]),"")[[1]]))
    p2_2 <-  sort(as.numeric(strsplit(as.character(m2[2]),"")[[1]]))
    
    
    pmap <- rbind(p1_1,p1_2,p2_1,p2_2)
    
    if(sum(p1_1)==4||sum(p1_1)==8||sum(p1_2)==4||sum(p1_2)==8){
      index22 <- matrix(c(1,2,3,4),ncol=4,byrow = T)
    }else{
      index22 <- ep1
    }
    
    if(sum(p2_1)==4||sum(p2_1)==8||sum(p2_2)==4||sum(p2_2)==8){
      index44 <- matrix(c(1,2,3,4),ncol=4,byrow = T)
    }else{
      index44 <- ep1
    }
    
    
    n1 <- dim(index22)[1]
    n2 <- dim(index44)[1]
    tt1 <- c()
    for(ii in 1:n1){
      for(jj in 1:n2){
        index <- rbind(index1,as.numeric(index22[ii,]),index3,as.numeric(index44[jj,]))
        tt <- two_rfz_phase_g2(M1=m11,M2=m22,pmap=pmap,index=index)
        if(tt[3]>0.5&&tt[3]<0.99999){
          tt[3] <- 1-tt[3]
        }
        tt1 <- rbind(tt1,c(ii,jj,tt))
      }
    }
    tt1 <- matrix(tt1[which(tt1[,5]>0),],ncol=6,byrow = F)
    tt2 <- matrix(tt1[which(tt1[,5]==min(tt1[,5])),],ncol=6,byrow = F)
    tt3 <- matrix(tt2[which(tt2[,6]==max(tt2[,6])),],ncol=6,byrow = F)[1,]
    
    phf1[i+1,] <- index22[tt3[1],]
    phf2[i+1,] <- index44[tt3[2],]
    
    abr <- rbind(abr,tt3[-c(1:2)])
    
    cat("m= ",i,"p1= ",as.numeric(index22[tt3[1],]),"p2= ",as.numeric(index44[tt3[2],]),"para= ",tt3[-c(1:2)],"\n")
    
  }
  
  retl <- list(ph1=phf1,ph2=phf2,para=abr)
  return(retl)
}