two_rfz_phase_g2 <- function(M1,M2,pmap,index){
  
  P1f <- two_phase(M1=index[1,],M2=index[2,])
  P2f <- two_phase(M1=index[3,],M2=index[4,])
  
  M <- M1*10000+M2
  
  GMZ <- FG1(P1M1=pmap[1,],P1M2=pmap[2,],P2M1=pmap[3,],P2M2=pmap[4,])
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
    
    if(all(abs(f_3-f_1)<1e-5)||iter>130)
      break
    
    f_1 <- f_3
    iter <- iter +1
  }
  
  
  iter <- 1
  r <- 0.01
  while(1){
    v1 <- r^2/(9-18*r+10*r^2)
    v2 <- r/(3-2*r)
    
    r1 <- (f_3[3]+f_3[5]+2*(f_3[2]+f_3[4]+f_3[6]+f_3[9])+2*v1*f_3[7]+(1+v2)*f_3[8])/2

    if(abs(r1-r)<1e-5||iter>150)
      break
    
    if(r1>0.25)
      r1 <- abs(0.5-r1)
    r <- r1
    iter <- iter +1
    

  }
  
  a <- sum(f_3[1:4])
  b <- sum(f_3[c(1,2,5,6)])
  ii <- f_3!=0
  #
  LL <- sum(E1[ii]*log(f_3[ii]),nan.rm=T)
  res <- c(a,b,r*0.77,LL)
  res
}




two_full_cal <- function(M1,M2,pm){
  
  
  
  #npm <- matrix(c(1,1,2,2,
  #                1,1,2,2,
  #                1,1,2,2,
  #                1,1,2,2), nrow = 4, byrow=T)
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
  if(ret22[7] >0.25)
    ret22[7] <- abs(0.5-ret22[7])
  return(ret22)
  
}




auto_scan_mc <- function(dat,se=c(1,100)){
  
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
    
    ret1 <- two_full_cal(M1=m11,M2=m22,pm=pm)
    cat("M1: ",mname[nc[1,i]]," M2: ",mname[nc[2,i]],"\n")
    cat("",ret1,"\n")
    
    fret <- rbind(fret,ret1)
  }
  
  return(fret)
}