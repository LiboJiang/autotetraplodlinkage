

allele_com <- function(M){
  
  mc <- expand.grid(1:4,1:4,1:4,1:4)
  
  mc1 <- apply(mc,1,function(x){
    if(max(table(x))==1)
      return(1)
    else{
      return(0)
    }
  })
  
  mc2 <- as.matrix(mc[as.logical(mc1),])
  
  mc2[which(mc2==1)] <- M[1]
  mc2[which(mc2==2)] <- M[2]
  mc2[which(mc2==3)] <- M[3]
  mc2[which(mc2==4)] <- M[4]
  
  tm1 <- unique(paste(mc2[,1],mc2[,2],mc2[,3],mc2[,4],sep=""))
  return(as.numeric(tm1))
}


allele_com_z <- function(phase){
  
  n <- dim(phase)[1]
  
  gam1s <- list()
  for(i in 1:n){
    gam1s[[i]] <- allele_com(M=phase[i,])
  }
  
  phm1 <- expand.grid(gam1s[[1]],gam1s[[2]])
  
  phm2 <- expand.grid(gam1s[[3]],gam1s[[4]])
  
  
  CV1 <- c()
  for(i in 1:dim(phm1)[1]){
    v1 <- as.numeric(strsplit(as.character(phm1[i,1]),"")[[1]])
    v2 <- as.numeric(strsplit(as.character(phm1[i,2]),"")[[1]])
    CV1 <- rbind(CV1,phase_tran(v1=v1,v2=v2))
  }
  
  CV2 <- c()
  for(i in 1:dim(phm2)[1]){
    v1 <- as.numeric(strsplit(as.character(phm2[i,1]),"")[[1]])
    v2 <- as.numeric(strsplit(as.character(phm2[i,2]),"")[[1]])
    CV2 <- rbind(CV2,phase_tran(v1=v1,v2=v2))
  }
  
  
  cp11 <- paste(CV1[,1],CV1[,2],CV1[,3],CV1[,4],sep="")
  cp12 <- paste(CV1[,5],CV1[,6],CV1[,7],CV1[,8],sep="")
  cp21 <- paste(CV2[,1],CV2[,2],CV2[,3],CV2[,4],sep="")
  cp22 <- paste(CV2[,5],CV2[,6],CV2[,7],CV2[,8],sep="")
  
  si1 <- !duplicated(paste(cp11,cp12,sep=""))
  si2 <- !duplicated(paste(cp21,cp22,sep=""))
  
  res11 <- phm1[si1,]; res12 <- phm2[si2,]
  res21 <- cbind(cp11,cp12)[si1,]; res22 <- as.matrix(cbind(cp21,cp22))[si2,];
  if(is.null(dim(res11)[1])){
    res11 <- t(as.matrix(res11))
  }
  if(is.null(dim(res12)[1])){
    res12 <- t(as.matrix(res12))
  }
  if(is.null(dim(res21)[1])){
    res21 <- t(as.matrix(res21))
  }
  if(is.null(dim(res22)[1])){
    res22 <- t(as.matrix(res22))
  }
  tp11 <- paste(res11[,1],res11[,2],sep="")
  tp12 <- paste(res12[,1],res12[,2],sep="")
  
  tp21 <- paste(res21[,1],res21[,2],sep="")
  tp22 <- paste(res22[,1],res22[,2],sep="")
  
  res1 <- expand.grid(tp11,tp12)
  res2 <- expand.grid(tp21,tp22)
  
  
  
  fres1 <- cbind(signif(as.numeric(as.character(res1[,1]))/10000,4),
                 (as.numeric(as.character(res1[,1]))%%10000),
                 signif(as.numeric(as.character(res1[,2]))/10000,4),
                 (as.numeric(as.character(res1[,2]))%%10000))
  
  fres2 <- cbind(signif(as.numeric(as.character(res2[,1]))/10000,4),
                 (as.numeric(as.character(res2[,1]))%%10000),
                 signif(as.numeric(as.character(res2[,2]))/10000,4),
                 (as.numeric(as.character(res2[,2]))%%10000))
  return(list(fres1=fres1,fres2=fres2))
}


phase_tran <- function(v1,v2){
  
  stad <- c(1,2,3,4)
  index <- v1==v2
  if(sum(index)==4||sum(index)==3){
    v1 <- stad
    v2 <- stad
  }
  if(sum(index)==2){
    v1[index] <- stad[index]
    v2[index] <- stad[index] 
    v1[!index] <- stad[!index]
    v2[!index] <- rev(stad[!index])
  }
  if(sum(index)==1){
    v1[index] <- stad[index]
    v2[index] <- stad[index] 
    v1[!index] <- stad[!index][c(1,2,3)]
    v2[!index] <- stad[!index][c(3,1,2)]
  }
  if(sum(index)==0){
    
    v1 <- c(1,2,3,4)
    v2 <- c(2,1,4,3)
  }
  if(is.null(dim(v1))){
    v <- cbind(t(as.matrix(v1)),t(as.matrix(v2)))
  }else{
    v <- cbind(v1,v2)
  }
  
  return(v)
}



phase_tran1 <- function(v1,v2){
  
  stad <- c(1,2,3,4)
  index <- v1==v2
  if(sum(index)==4||sum(index)==3){
    v1 <- stad
    v2 <- stad
  }
  if(sum(index)==2){
    v1[index] <- stad[index]
    v2[index] <- stad[index] 
    v1[!index] <- stad[!index]
    v2[!index] <- rev(stad[!index])
  }
  if(sum(index)==1){
    re1 <- paire(stad[!index])
    v11 <- c();v12 <- c()
    for(i in 1:dim(re1)[1]){
      v1[index] <- stad[index]
      v2[index] <- stad[index]
      v1[!index] <- stad[!index]
      v2[!index] <- re1[i,]
      v11 <- rbind(v11,v1);v12 <- rbind(v12,v2)
    }
    v1 <- v11
    v2 <- v12
  }
  if(sum(index)==0){
    re1 <- paire(stad[!index])
    v11 <- c();v12 <- c()
    for(i in 1:dim(re1)[1]){
      v1[index] <- stad[index]
      v2[index] <- stad[index]
      v1[!index] <- stad[!index]
      v2[!index] <- re1[i,]
      v11 <- rbind(v11,v1);v12 <- rbind(v12,v2)
    }
    v1 <- v11
    v2 <- v12
  }
  if(is.null(dim(v1))){
    v <- cbind(t(as.matrix(v1)),t(as.matrix(v2)))
  }else{
    v <- cbind(v1,v2)
  }
  
  return(v)
}




paire <- function(V){
  if(length(V)==4){
    vt <- as.matrix(expand.grid(V,V,V,V))
  }
  if(length(V)==3){
    vt <- as.matrix(expand.grid(V,V,V))
  }
  
  vi1 <- apply(vt,1,function(x){
    if(max(table(x))==1)
      return(1)
    else{
      return(0)
    }
  })
  vi2 <- as.matrix(vt[as.logical(vi1),])
  
  vi3 <- apply(vi2,1,function(x,y){
    sum(x!=y)
  },y=V)
  vi4 <- vi2[which(vi3==length(V)),]
  return(vi4)
}

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
  
  #m11 <- matrix(as.numeric(unlist(strsplit(as.character(M1),""))),ncol=4,byrow=T)
  #m22 <- matrix(as.numeric(unlist(strsplit(as.character(M2),""))),ncol=4,byrow=T)
  #M <-(m11[,1]*1000+m11[,4]*100+m11[,2]*10+m11[,3])*10000+(m22[,1]*1000+m22[,4]*100+m22[,2]*10+m22[,3])
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

Haldan <- function(rf){
  
  
  -0.5 * log(1 - 2 * rf)
  
  
}


two_full_cal <- function(M1,M2,pm){
  
 
  tph <- allele_com_z(phase=pm)
  allph1 <- tph$fres1
  allph2 <- tph$fres2

  
  nph <- dim(allph1)[1]
  
  cat("Phase="," ",nph,"\n")
  
  ret1 <- c()
  for(i in 1:nph){
    

    
    ph21 <- as.numeric(strsplit(as.character(allph2[i,1]),"")[[1]])
    ph22 <- as.numeric(strsplit(as.character(allph2[i,2]),"")[[1]])
    ph23 <- as.numeric(strsplit(as.character(allph2[i,3]),"")[[1]])
    ph24 <- as.numeric(strsplit(as.character(allph2[i,4]),"")[[1]])
    
    ph2 <- rbind(ph21,ph22,ph23,ph24)
    
    #cat("Phase="," ",i,"\n")
    ret <- two_rfz_phase_g2(M1=M1,M2=M2,pmap=pm,index=ph2)
    
    #cat(allph1[i,1],"",allph1[i,2],"",allph1[i,3],"",allph1[i,4],ret,"\n")
    
    ret1 <- rbind(ret1,c(as.numeric(allph1[i,1:4]),ret,as.numeric(allph2[i,1:4])))
  }
  
  scr <- ret1[,7]
  lz <- which(scr<0)
  if(length(lz)>0){
    ret1 <- ret1[-lz,]
  }else{
    ret1 <- ret1
  }
  if(is.null(dim(ret1))){
    ret1 <- t(as.matrix(ret1))
  }
  ret11 <- ret1[which(ret1[,7]==min(ret1[,7])),]
  if(is.matrix(ret11)){
    ret22 <- ret11[which(ret11[,8]==max(ret11[,8]))[1],]
  }else{
    ret22 <- ret11
  }       
  
  if(ret22[7] >0.75)
    ret22[7] <- 1-ret22[7]
  return(ret22)
  
}



auto_scan <- function(dat,se=c(1,100)){
  
  if(!is.null(se)){
    
    ndat <- dat[se,]
  }else{
    ndat <- dat
  }
  
  n <- dim(ndat)[1]
  mname <- rownames(ndat)
  if(n==2){
    nc <- as.matrix(c(1,2))
  }else{
    nc <- combn(n,2)
  }

  fret <- c()
  for(i in 1:dim(nc)[2]){
    
    m1 <- ndat[nc[1,i],]
    m2 <- ndat[nc[2,i],]
    m11 <- m1[-c(1,2)]
    m22 <- m2[-c(1,2)]
    
    P1M1 <- as.numeric(strsplit(as.character(m1[1]),"")[[1]])
    P1M2 <- as.numeric(strsplit(as.character(m2[1]),"")[[1]])
    P2M1 <- as.numeric(strsplit(as.character(m1[2]),"")[[1]])
    P2M2 <- as.numeric(strsplit(as.character(m2[2]),"")[[1]])
    pm <- matrix(c(P1M1,P1M2,P2M1,P2M2),nrow=4,byrow = T)
    
    cat("M1: ",mname[nc[1,i]]," M2: ",mname[nc[2,i]],"\n")
    ret <- two_full_cal(M1=m11,M2=m22,pm=pm)
    ret1 <- c(mname[nc[1,i]],mname[nc[2,i]],ret)
    cat("",ret,"\n")
    fret <- rbind(fret,ret1)
  }
  
  return(fret)
}





auto_scan_ph <- function(dat,se=c(1,100)){
  
  if(!is.null(se)){
    
    ndat <- dat[se,]
  }else{
    ndat <- dat
  }
  
  n <- dim(ndat)[1]
  mname <- rownames(ndat)
  if(n==2){
    nc <- as.matrix(c(1,2))
  }else{
    nc <- combn(n,2)
  }
  
  fret <- c()
  for(i in 1:dim(nc)[2]){
    
    m1 <- ndat[nc[1,i],]
    m2 <- ndat[nc[2,i],]
    m11 <- m1[-c(1,2)]
    m22 <- m2[-c(1,2)]
    
    P1M1 <- as.numeric(strsplit(as.character(m1[1]),"")[[1]])
    P1M2 <- as.numeric(strsplit(as.character(m2[1]),"")[[1]])
    P2M1 <- as.numeric(strsplit(as.character(m1[2]),"")[[1]])
    P2M2 <- as.numeric(strsplit(as.character(m2[2]),"")[[1]])
    pm <- matrix(c(P1M1,P1M2,P2M1,P2M2),nrow=4,byrow = T)
    
    if(all(P2M1==1)){
      pm[1,] <- P2M1
      pm[3,] <- P1M1
    }
    
    if(all(P2M2==1)){
      pm[2,] <- P2M2
      pm[4,] <- P1M2
    }
    
    
    cat("M1: ",mname[nc[1,i]]," M2: ",mname[nc[2,i]],"\n")
    ret <- two_full_cal_ph(M1=m11,M2=m22,pm=pm)
    #ret1 <- c(mname[nc[1,i]],mname[nc[2,i]],ret)
    #cat("",ret,"\n")
    fret <- rbind(fret,ret)
  }
  
  return(fret)
}




