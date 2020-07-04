
Haldan <- function(rf){
  
  
  -0.5 * log(1 - 2 * rf)*100
  
  
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




PF_zy1_1 <-function(g1,g2,P1f,P2f,fii){
  
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

frt <-function(fii){
  
  pp <- c()
  for(i in 1:length(fii)){
    pp <- rbind(pp,strsplit(as.character(fii[i]),"")[[1]])
  }
  
  fpm <- c()
  for(ii in 1:9){
    
    fpm <- cbind(fpm,((pp[,1]==ii)+(pp[,2]==ii))/2)
    
  }
  
  return(fpm)
  
}

findex <- function(){
  
  nM <- matrix(NA,9,9)
  for(i in 1:9){
    for(j in 1:9){
      nM[i,j] <- i*10+j
    }
  }
  return(nM)
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

allele_com <-function(M){
  
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


allele_com_z <-function(phase){
  
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