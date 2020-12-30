
three_rfz_phase <- function(M1,M2,M3,pmap,index){
  
  P1f <- three_phase(M1=index[1,],M2=index[2,],M3=index[3,])
  P2f <- three_phase(M1=index[4,],M2=index[5,],M3=index[6,])
  
  
  M <- M1*10000*10000+M2*10000+M3
  
  GMZ <- FG1_three(P1M1=pmap[1,],P1M2=pmap[2,],P1M3=pmap[3,],P2M1=pmap[4,],P2M2=pmap[5,],P2M3=pmap[6,])

  DLRZ <- TM_zy_three(g1=GMZ$zoc13,g2=GMZ$zoc2)
  
  
  
  n1 <- length(unique(GMZ$zoc13))
  n2 <- length(unique(GMZ$zoc2))

  count <- matrix(0,n1,n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      count[i,j] <- length(which(GMZ$zygote1[i,j]==M))
    }
  }
  
  cat("count= ",sum(count),"\n")

  f_1 <- rep(0.0001,59)
  
  fii <- c(t(findex()))
  frt1 <- frt(fii)
  
  PFCZ1 <- PF_zy_thr(g1=GMZ$zoc13,g2=GMZ$zoc2,P1f,P2f,fii)
  
  iter <- 1
  while(1){
    
    ft <- HP_thr(f=f_1,P1f=P1f,P2f=P2f)
    fc <- DLRZ$dL%*%ft$HP%*%t(DLRZ$dR)
    f_2 <- ft$f2
    
    EM1 <- f_2[as.numeric(PFCZ1$ncc)]*c(t(count))[PFCZ1$nkk]/c(t(fc))[PFCZ1$nkk]
    
    #cat("EM1= ",sum(EM1),"\n")
 
    EM2 <- PFCZ1$trm%*%as.matrix(EM1)
    
    E1 <- t(EM2)%*%frt1
    
    f_3 <- E1/sum(count)

    if(all(abs(f_3-f_1)<1e-5)||iter>230)
      break
    f_1 <- f_3
    #cat("iter = "," ",iter,"",f_1,"\n")
    iter <- iter +1
  }
  
  r <- c(0,0,0)
  iter <- 1
  while(1){
    v11 <- r[1]^2/(9-18*r[1]+10*r[1]^2)
    v21 <- r[2]^2/(9-18*r[2]+10*r[2]^2)
    v12 <- r[1]/(3-2*r[1])
    v22 <- r[2]/(3-2*r[2])
    
    r1 <- sum(f_3[c(2,4,6,7,10,11,13,14,17,18,20,21,24,25,28,32,33,35,37,41,44,47,48,54,56,59)])+
      sum(f_3[c(3,8,9,15,16,22,23,26,27,34,40,45,46,55)])/2+sum(f_3[c(29,38,42,49,57)])*v11+
      sum(f_3[c(30,31,36,39,43,50,51,52,53,58)])*(1+v12)/2
    
    
    r2 <- sum(f_3[c(2,4,5,7,9,11,14,18,19,21,23,27,28,31,32,34,35,38,39,41,44,46,48,53,55,57)])+
      sum(f_3[c(3,8,10,12,13,20,29,30,36,37,40,45,47,56)])/2+sum(f_3[c(15,24,42,50,59)])*v21+
      sum(f_3[c(16,17,22,25,43,49,50,51,52,54,58)])*(1+v22)/2
    
    r3 <- sum(f_3[c(5:11)])+sum(f_3[c(19:25)])+sum(f_3[c(33:39)])+sum(f_3[c(55:59,52)])+
      sum(f_3[c(12:18)])+sum(f_3[c(26:32)])+sum(f_3[c(40,41,44,49,50,51)])/2+3*sum(f_3[c(45,46,47,48,53,54)])/4+
      sum(f_3[c(42,42,49)])*v11/2+sum(f_3[c(42,42,50)])*v12/2-2*f_3[42]*v11*v12+
      sum(f_3[c(43,43,49,51)])*v21/2-f_3[52]*v21/2-f_3[50]*v12*v21+
      sum(f_3[c(43,43,49,51)])*v22/2-f_3[52]*v22/2-f_3[49]*v11*v22-sum(f_3[c(43,43,51)])*v21*v22+f_3[52]*v21*v22
    cr <- c(r1,r2,r3)
    if(max(abs(cr-r))<1e-5||max(cr)>0.9999999)
      break
    
    r <- cr
    iter <- iter +1
    
  }
  coc <- (r[1]+r[2]-r[3])/(2*r[1]*r[2])
  a <- sum(f_3[1:25])
  c <- sum(f_3[c(1,2,5:7,12:14,19:21,26:28,33:35,40,41,45:48,55,56)])
  b <- sum(f_3[c(1:11,26:39)])
  
  rnc <- rbind(r,(1-r))
  r <- apply(rnc, 2, min)
  #coc <- (r[1]+r[2]-r[3])/(2*r[1]*r[2])
  
  ii <- f_1!=0
  LL <- sapp_L(sum(E1))-sum(sapp_L(E1[ii]))+sum(E1[ii]*log(f_1[ii]),nan.rm=T)
  res1 <- c(a,b,c,r,coc,LL)
  rm(P1f,P2f,GMZ,DLRZ,PFCZ1,EM1,EM2)
  gc()
  res <- list(res1=res1,res2=f_1)
  return(res)
}



three_scan <- function(geno,ph1,ph2){
  
  
  
  n <- dim(geno)[1]
  
  geno[which(geno==0)] <- 1111
  geno[which(geno==1)] <- 1112
  geno[which(geno==2)] <- 1122
  geno[which(geno==3)] <- 1222
  geno[which(geno==4)] <- 2222
  
  para  <- c(); fl <- c()
  for(i in 1:300){
    
    il <- (2*i-1):(2*i+1)
    if(max(il)>n){
      break
    }
    gg <- geno[il,]
    
    m1 <- as.numeric(gg[1,])
    m2 <- as.numeric(gg[2,])
    m3 <- as.numeric(gg[3,])
    
    
    m11 <- m1[-c(1:2)]
    m22 <- m2[-c(1:2)]
    m33 <- m3[-c(1:2)]
    
    p1_1 <-  sort(as.numeric(strsplit(as.character(m1[1]),"")[[1]]))
    p1_2 <-  sort(as.numeric(strsplit(as.character(m2[1]),"")[[1]]))
    p1_3 <-  sort(as.numeric(strsplit(as.character(m3[1]),"")[[1]]))
    p2_1 <-  sort(as.numeric(strsplit(as.character(m1[2]),"")[[1]]))
    p2_2 <-  sort(as.numeric(strsplit(as.character(m2[2]),"")[[1]]))
    p2_3 <-  sort(as.numeric(strsplit(as.character(m3[2]),"")[[1]]))
    
    pmap <- rbind(p1_1,p1_2,p1_3,p2_1,p2_2,p2_3)
    
    phase1 <- ph1[il,];phase2 <- ph2[il,]
    
    index <- rbind(phase1,phase2)
    
    ret <- three_rfz_phase(M1=m11,M2=m22,M3=m33,pmap=pmap,index=index)
    para <- rbind(para,ret$res1)
    cat("para= ",ret$res1,"\n")
    fl <- rbind(fl,ret$res2)
  }
  
  retl <- list(para=para,fl=fl)
  return(retl)
}


three_scan1 <- function(geno){
  
  
  
  n <- dim(geno)[1]
  
  geno[which(geno==0)] <- 1111
  geno[which(geno==1)] <- 1112
  geno[which(geno==2)] <- 1122
  geno[which(geno==3)] <- 1222
  geno[which(geno==4)] <- 2222
  
  ep <- as.matrix(expand.grid(1:4,1:4,1:4,1:4))
  ep1 <- ep[which(apply(ep,1,function(x){length(unique(x))})==4),]
  phf1 <- phf2 <- matrix(NA,nrow=n,ncol=4)
  
  
  
  para  <- c(); 
  for(i in 1:300){
    
    il <- (2*i-1):(2*i+1)
    if(max(il)>n){
      break
    }
    gg <- geno[il,]
    
    m1 <- as.numeric(gg[1,])
    m2 <- as.numeric(gg[2,])
    m3 <- as.numeric(gg[3,])
    
    
    m11 <- m1[-c(1:2)]
    m22 <- m2[-c(1:2)]
    m33 <- m3[-c(1:2)]
    
    p1_1 <-  sort(as.numeric(strsplit(as.character(m1[1]),"")[[1]]))
    p1_2 <-  sort(as.numeric(strsplit(as.character(m2[1]),"")[[1]]))
    p1_3 <-  sort(as.numeric(strsplit(as.character(m3[1]),"")[[1]]))
    p2_1 <-  sort(as.numeric(strsplit(as.character(m1[2]),"")[[1]]))
    p2_2 <-  sort(as.numeric(strsplit(as.character(m2[2]),"")[[1]]))
    p2_3 <-  sort(as.numeric(strsplit(as.character(m3[2]),"")[[1]]))
    
    pmap <- rbind(p1_1,p1_2,p1_3,p2_1,p2_2,p2_3)
    
    index1 <- as.numeric(phf1[il[1],])
    index4 <- as.numeric(phf2[il[1],])
    
    
    if(sum(p1_1)==4||sum(p1_1)==8||sum(p1_2)==4||sum(p1_2)==8){
      index22 <- matrix(c(1,2,3,4),ncol=4,byrow = T)
    }else{
      index22 <- ep1
    }
    
    if(sum(p1_2)==4||sum(p1_2)==8||sum(p1_3)==4||sum(p1_3)==8){
      index33 <- matrix(c(1,2,3,4),ncol=4,byrow = T)
    }else{
      index33 <- ep1
    }
    
    
    if(sum(p2_1)==4||sum(p2_1)==8||sum(p2_2)==4||sum(p2_2)==8){
      index55 <- matrix(c(1,2,3,4),ncol=4,byrow = T)
    }else{
      index55 <- ep1
    }
    
    if(sum(p2_2)==4||sum(p2_2)==8||sum(p2_3)==4||sum(p2_3)==8){
      index66 <- matrix(c(1,2,3,4),ncol=4,byrow = T)
    }else{
      index66 <- ep1
    }
    
    
    n1 <- dim(index22)[1]
    n2 <- dim(index33)[1]
    n3 <- dim(index55)[1]
    n4 <- dim(index66)[1]
    tt1 <- c()
    for(ii in 1:n1){
      for(jj in 1:n2){
        for(kk in 1:n3){
          for(hh in 1:n4){
            index <- rbind(index1,as.numeric(index22[ii,]),as.numeric(index33[ii,]),
                           index3,as.numeric(index55[jj,]),as.numeric(index66[jj,]))
            tt <- three_rfz_phase(M1=m11,M2=m22,M3=m33,pmap=pmap,index=index)
            tt1 <- rbind(tt1,c(ii,jj,kk,hh,tt$res1))
          }
        }
      }
    }
    
    rt <- tt1[,8:10]
    rts <- rowSums(rt)
    tt3 <- tt1[which(rts==min(rts)),]
    
    
    para <- rbind(para,tt3[-c(1:4)])
    
    cat("m= ",i,"para= ",tt3[-c(1:4)],"\n")
    
    phf1[il[2],] <- index22[tt3[1],]
    phf1[il[3],] <- index33[tt3[2],]
    
    phf2[il[2],] <- index55[tt3[1],]
    phf2[il[3],] <- index66[tt3[2],]
    
  }
  return(para)
}

three_scan2 <- function(geno,ph1,ph2,iii){
  
  
  
  n <- dim(geno)[1]
  
  geno[which(geno==0)] <- 1111
  geno[which(geno==1)] <- 1112
  geno[which(geno==2)] <- 1122
  geno[which(geno==3)] <- 1222
  geno[which(geno==4)] <- 2222
  
  para  <- c(); fl <- c()
  for(i in iii){
    
    il <- (2*i-1):(2*i+1)
    if(max(il)>n){
      break
    }
    gg <- geno[il,]
    
    m1 <- as.numeric(gg[1,])
    m2 <- as.numeric(gg[2,])
    m3 <- as.numeric(gg[3,])
    
    
    m11 <- m1[-c(1:2)]
    m22 <- m2[-c(1:2)]
    m33 <- m3[-c(1:2)]
    
    p1_1 <-  sort(as.numeric(strsplit(as.character(m1[1]),"")[[1]]))
    p1_2 <-  sort(as.numeric(strsplit(as.character(m2[1]),"")[[1]]))
    p1_3 <-  sort(as.numeric(strsplit(as.character(m3[1]),"")[[1]]))
    p2_1 <-  sort(as.numeric(strsplit(as.character(m1[2]),"")[[1]]))
    p2_2 <-  sort(as.numeric(strsplit(as.character(m2[2]),"")[[1]]))
    p2_3 <-  sort(as.numeric(strsplit(as.character(m3[2]),"")[[1]]))
    
    pmap <- rbind(p1_1,p1_2,p1_3,p2_1,p2_2,p2_3)
    
    phase1 <- ph1[il,];phase2 <- ph2[il,]
    
    index <- rbind(phase1,phase2)
    
    ret <- three_rfz_phase(M1=m11,M2=m22,M3=m33,pmap=pmap,index=index)
    para <- rbind(para,ret$res1)
    cat("i= ",i,"para= ",ret$res1,"\n")
    fl <- rbind(fl,ret$res2)
  }
  
  retl <- list(para=para,fl=fl)
} 

get_con_param<-function(parm.id)
{
  for (e in commandArgs())
  {
    ta = strsplit(e,"=", fixed=TRUE);
    if(! is.na( ta[[1]][2]))
    {
      temp = ta[[1]][2];
      if( ta[[1]][1] == parm.id) {
        return (as.character(temp));
      }
    }
  }
  
  return(NA);
}





fr <- function(f_3){
  
  r <- c(0,0,0)
  iter <- 1
  while(1){
    v11 <- r[1]^2/(9-18*r[1]+10*r[1]^2)
    v21 <- r[2]^2/(9-18*r[2]+10*r[2]^2)
    v12 <- r[1]/(3-2*r[1])
    v22 <- r[2]/(3-2*r[2])
    
    r1 <- sum(f_3[c(2,4,6,7,10,11,13,14,17,18,20,21,24,25,28,32,33,35,37,41,44,47,48,54,56,59)])+
      sum(f_3[c(3,8,9,15,16,22,23,26,27,34,40,45,46,55)])/2+sum(f_3[c(29,38,42,49,57)])*v11+
      sum(f_3[c(30,31,36,39,43,50,51,52,53,58)])*(1+v12)/2
    
    
    r2 <- sum(f_3[c(2,4,5,7,9,11,14,18,19,21,23,27,28,31,32,34,35,38,39,41,44,46,48,53,55,57)])+
      sum(f_3[c(3,8,10,12,13,20,29,30,36,37,40,45,47,56)])/2+sum(f_3[c(15,24,42,50,59)])*v21+
      sum(f_3[c(16,17,22,25,43,49,50,51,52,54,58)])*(1+v22)/2
    
    r3 <- sum(f_3[c(5:11)])+sum(f_3[c(19:25)])+sum(f_3[c(33:39)])+sum(f_3[c(55:59,52)])+
      sum(f_3[c(12:18)])+sum(f_3[c(26:32)])+sum(f_3[c(40,41,44,49,50,51)])/2+3*sum(f_3[c(45,46,47,48,53,54)])/4+
      sum(f_3[c(42,42,49)])*v11/2+sum(f_3[c(42,42,50)])*v12/2-2*f_3[42]*v11*v12+
      sum(f_3[c(43,43,49,51)])*v21/2-f_3[52]*v21/2-f_3[50]*v12*v21+
      sum(f_3[c(43,43,49,51)])*v22/2-f_3[52]*v22/2-f_3[49]*v11*v22-sum(f_3[c(43,43,51)])*v21*v22+f_3[52]*v21*v22
    cr <- c(r1,r2,r3)
    if(max(abs(cr-r))<1e-5||max(cr)>0.9999999)
      break
    
    r <- cr
    iter <- iter +1
    
  }
  coc <- (r[1]+r[2]-r[3])/(2*r[1]*r[2])
  a <- sum(f_3[1:25])
  c <- sum(f_3[c(1,2,5:7,12:14,19:21,26:28,33:35,40,41,45:48,55,56)])
  b <- sum(f_3[c(1:11,26:39)])
  
  rnc <- rbind(r,(1-r))
  r <- apply(rnc, 2, min)
  
  #ii <- f_3!=0
  #LL <- sapp_L(sum(E1))-sum(sapp_L(E1[ii]))+sum(E1[ii]*log(f_1[ii]),nan.rm=T)
  res1 <- c(a,b,c,r,coc)
  
  res1
}
