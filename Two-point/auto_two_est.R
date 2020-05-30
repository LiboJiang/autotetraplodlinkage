



D <- function(r,g1,g2){
  
  mode1 <- c(11,22,33,44,12,13,14,23,24,34)
  
  v1 <- 2*r^2/(9-18*r+10*r^2)
  v2 <- 1+r/(3-2*r)
  
  dm <- matrix(c(0,2,2,2,1,1,1,2,2,2,
               2,0,2,2,1,2,2,1,1,2,
               2,2,0,2,2,1,2,1,2,1,
               2,2,2,0,2,2,1,2,1,1,
               1,1,2,2,v1,v2,v2,v2,v2,2,
               1,2,1,2,v2,v1,v2,v2,2,v2,
               1,2,2,1,v2,v2,v1,2,v2,v2,
               2,1,1,2,v2,v2,2,v1,v2,v2,
               2,1,2,1,v2,2,v2,v2,v1,v2,
               2,2,1,1,2,v2,v2,v2,v2,v1),nrow=10,byrow=T)
 
  
  g11 <- unique(g1)
  g22 <- unique(g2)
  dm1 <- dm[match(g11,mode1),match(g22,mode1)]
  
  
  return(dm1)
}


HP <- function(f){
  
  f1 <- f[1]/4;f2 <- f[2]/12;f3 <- f[3]/12;f4 <- f[4]/12;f5 <- f[5]/12;
  f6 <- f[6]/12;f7 <- f[7]/6;f8 <- f[8]/24;f9 <- f[9]/6
  
  ind <- FM()
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

TM <- function(g1,g2){
  
  g11 <- unique(g1)
  g22 <- unique(g2)
  
  n1 <- length(g11)
  n2 <- length(g22)
  
  Tr1 <- matrix(0,n1,10)
  Tr2 <- matrix(0,n2,10)
  
  for(i in 1:n1){
    Tr1[i,which(g11[i]==g1)] <- 1
  }
  
  for(i in 1:n2){
    Tr2[i,which(g22[i]==g2)] <- 1
  }
  return(list(dL=Tr1,dR=Tr2))
}


PF <- function(g1,g2){
  
  ind <- FM()
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



two_rf1 <- function(M,M1,M2){
  
  GM <- FG(M1=M1,M2=M2)
  DLR <- TM(g1=GM$m1,g2=GM$m2)
  PFC <- PF(g1=GM$m1,g2=GM$m2)
  n1 <- length(unique(GM$m1))
  n2 <- length(unique(GM$m2))
  count <- matrix(0,n1,n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      count[i,j] <- length(which(GM$game1[i,j]==M))
    }
  }
  
  
  #ind <- FM()
  
  #N1 <- c()
  
  #for(i in 1:9){
  #  N1 <- c(N1,sum(count[which(ind==i)]))
  #}
  f_1 <- rep(1/9,9)
  f_2 <- f_1*c(1/4,1/12,1/12,1/12,1/12,1/12,1/6,1/24,1/6)
  #r <- 0.05
  iter <- 1
  while(1){
    
    f_2 <- f_1*c(1/4,1/12,1/12,1/12,1/12,1/12,1/6,1/24,1/6)
    ft <- HP(f_1)
    fc <- DLR$dL%*%ft%*%t(DLR$dR)
    fe <-rep(0,9)
    for(k in 1:9){
      tmp <- 0
      for(i in 1:n1){
        for(j in 1:n2){
          if(fc[i,j]==0){
            next
          }else{
            tmp <- c(tmp,sum(as.numeric(strsplit(PFC[i,j],"")[[1]])==k)*f_2[k]*count[i,j]/fc[i,j])
          }
        }
      }
      fe[k] <- sum(tmp)
    }
    
    
    f_3 <- fe/sum(count)
    
    #E1 <- D(r=r,g1=GM$m1,g2=GM$m2)
    
    #r1 <- sum(E1*count)/(2*sum(count))
    
    if(all(abs(f_3-f_1)<1e-5))
      break
    
    f_1 <- f_3
    cat("iter = "," ",iter,"\n")
    iter <- iter +1
  }
  
  r <- 0.05
  while(1){
    v1 <- r^2/(9-18*r+10*r^2)
    v2 <- r/(3-2*r)
    
    r1 <- (f_3[3]+f_3[5]+2*(f_3[2]+f_3[4]+f_3[6]+f_3[9])+2*v1*f_3[7]+(1+v2)*f_3[8])/2
    
    #r1 <- sum(E1*count)/(2*sum(count))
    
    if(abs(r1-r)<1e-5)
      break
    
    r <- r1
    
  }
  
  a <- sum(f_1[1:4])
  b <- sum(f_1[c(1,2,5,6)])
  #a
  #b
  
  
  #r
  res <- c(a,b,r)
  return(res)
}

####zygote########






HP1 <- function(f){
  
  f1 <- f[1]/4;f2 <- f[2]/12;f3 <- f[3]/12;f4 <- f[4]/12;f5 <- f[5]/12;
  f6 <- f[6]/12;f7 <- f[7]/6;f8 <- f[8]/24;f9 <- f[9]/6
  
  ind <- FM()
  ind[which(ind==1)] <- f1
  ind[which(ind==2)] <- f2
  ind[which(ind==3)] <- f3
  ind[which(ind==4)] <- f4
  ind[which(ind==5)] <- f5
  ind[which(ind==6)] <- f6
  ind[which(ind==7)] <- f7
  ind[which(ind==8)] <- f8
  ind[which(ind==9)] <- f9
  
  tmr1 <- matrix(1,nrow=10,ncol=10)
  
  tmp2 <- kronecker(ind,ind)
  return(tmp2)
  
  
  
}




PF_zy <- function(g1,g2){
  
  ind <- FM11()
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








findex <- function(){
  
  nM <- matrix(NA,9,9)
  for(i in 1:9){
    for(j in 1:9){
      nM[i,j] <- i*10+j
    }
  }
  return(nM)
}


frt <- function(fii){
  
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









FM1 <- function(){
  
  tmp <- matrix(c(1,2,2,2,5,5,5,6,6,6,
                  2,1,2,2,5,6,6,5,5,6,
                  2,2,1,2,6,5,6,5,6,5,
                  2,2,2,1,6,6,5,6,5,5,
                  3,3,4,4,7,8,8,8,8,9,
                  3,4,3,4,8,7,8,8,9,8,
                  3,4,4,3,8,8,7,9,8,8,
                  4,3,3,4,8,8,9,7,8,8,
                  4,3,4,3,8,9,8,8,7,8,
                  4,4,3,3,9,8,8,8,8,7),nrow=10,ncol=10,byrow=T)
  
  tmr <- matrix(1,nrow=10,ncol=10)
  
  tmp1 <- kronecker(tmp,tmr)
  
  return(tmp1)
}



FM11 <- function(){
  
  tmp <- matrix(c(1,2,2,2,5,5,5,6,6,6,
                  2,1,2,2,5,6,6,5,5,6,
                  2,2,1,2,6,5,6,5,6,5,
                  2,2,2,1,6,6,5,6,5,5,
                  3,3,4,4,7,8,8,8,8,9,
                  3,4,3,4,8,7,8,8,9,8,
                  3,4,4,3,8,8,7,9,8,8,
                  4,3,3,4,8,8,9,7,8,8,
                  4,3,4,3,8,9,8,8,7,8,
                  4,4,3,3,9,8,8,8,8,7),nrow=10,ncol=10,byrow=T)
  tmp1 <- matrix(NA,100,100)
  
  for(i in 1:10){
    for(j in 1:10){
      tmp1[((i-1)*10+1):(i*10),((j-1)*10+1):(j*10)] <- tmp[i,j]*10+tmp
    }
  }
  
  return(tmp1)
}




two_rfz <- function(M1,M2,P1M1,P1M2,P2M1,P2M2){
  
  M <- M1*10000+M2
  GMZ <- FG1(P1M1=P1M1,P1M2=P1M2,P2M1=P2M1,P2M2=P2M2)
  DLRZ <- TM_zy(g1=GMZ$zoc1,g2=GMZ$zoc2)
  PFCZ <- PF_zy(g1=GMZ$zoc1,g2=GMZ$zoc2)
  
  
  n1 <- length(unique(GMZ$zoc1))
  n2 <- length(unique(GMZ$zoc2))
  count <- matrix(0,n1,n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      count[i,j] <- length(which(GMZ$zygote1[i,j]==M))
    }
  }
  
  f_1 <- rep(1/9,9)
  f_2 <- f_1*c(1/4,1/12,1/12,1/12,1/12,1/12,1/6,1/24,1/6)
  fii <- c(t(findex()))
  frt1 <- frt(fii)
  #r <- 0.05
  iter <- 1
  while(1){
    
    f_11 <- as.matrix(f_1*c(1/4,1/12,1/12,1/12,1/12,1/12,1/6,1/24,1/6))
    f_2 <- c(t(kronecker(f_11,t(f_11))))
    ft <- HP1(f_1)
    fc <- DLRZ$dL%*%ft%*%t(DLRZ$dR)
    fe <- rep(0,length(f_2))
    for(k in 1:length(f_2)){
      tmp <- 0
      for(i in 1:n1){
        for(j in 1:n2){
          if(fc[i,j]==0){
            next
          }else{
            tmp <- c(tmp,sum(as.numeric(strsplit(PFCZ[i,j],"-")[[1]])==fii[k],na.rm=T)*f_2[k]*count[i,j]/fc[i,j])
          }
        }
      }
      fe[k] <- sum(tmp)
    }
    
    
    f_3 <- t(as.matrix(fe))%*%frt1/sum(count)
    
    #E1 <- D(r=r,g1=GM$m1,g2=GM$m2)
    
    #r1 <- sum(E1*count)/(2*sum(count))
    
    if(all(abs(f_3-f_1)<1e-5))
      break
    
    f_1 <- f_3
    cat("iter = "," ",iter,"",f_1,"\n")
    iter <- iter +1
  }
  
  
  
  r <- 0.05
  while(1){
    v1 <- r^2/(9-18*r+10*r^2)
    v2 <- r/(3-2*r)
    
    r1 <- (f_3[3]+f_3[5]+2*(f_3[2]+f_3[4]+f_3[6]+f_3[9])+2*v1*f_3[7]+(1+v2)*f_3[8])/2
    
    #r1 <- sum(E1*count)/(2*sum(count))
    
    if(abs(r1-r)<1e-5)
      break
    
    r <- r1
    
  }
  
  a <- sum(f_1[1:4])
  b <- sum(f_1[c(1,2,5,6)])
  res <- c(a,b,r)
  
  
  return(res)
}





