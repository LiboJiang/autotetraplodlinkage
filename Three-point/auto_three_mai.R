
three_rfz_phase <- function(M1,M2,M3,pmap,index){
  
  P1f <- three_phase(M1=index[1,],M2=index[2,],M3=index[3,])
  P2f <- three_phase(M1=index[4,],M2=index[5,],M3=index[6,])
  
  
  #m11 <- matrix(as.numeric(unlist(strsplit(as.character(M1),""))),ncol=4,byrow=T)
  #m22 <- matrix(as.numeric(unlist(strsplit(as.character(M2),""))),ncol=4,byrow=T)
  #m33 <- matrix(as.numeric(unlist(strsplit(as.character(M3),""))),ncol=4,byrow=T)
  
  
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
  
 
  f_1 <- c(0,rep(1/59,59))
  
  fii <- c(t(findex()))
  frt1 <- frt(fii)
  
  PFCZ1 <- PF_zy_thr(g1=GMZ$zoc13,g2=GMZ$zoc2,P1f,P2f,fii)
  
  iter <- 1
  while(1){
    
    ft <- HP_thr(f=f_1,P1f=P1f,P2f=P2f)
    fc <- DLRZ$dL%*%ft$HP%*%t(DLRZ$dR)
    f_2 <- ft$f2
    
    EM1 <- f_2[as.numeric(PFCZ1$ncc)]*c(t(count))[PFCZ1$nkk]/c(t(fc))[PFCZ1$nkk]
    
    #EM1[which(is.na(EM1))]  <- 0
    #EM1[which(is.nan(EM1))]  <- 0
    #EM1[which(is.infinite(EM1))] <- 0
    
    
    EM2 <- PFCZ1$trm%*%as.matrix(EM1)
    
    E1 <- t(EM2)%*%frt1
    
    f_3 <- E1/sum(count)
    
    #E1 <- D(r=r,g1=GM$m1,g2=GM$m2)
    
    #r1 <- sum(E1*count)/(2*sum(count))
    
    if(all(abs(f_3-f_1)<1e-5)||iter>130)
      break
    
    f_1 <- f_3
    cat("iter = "," ",iter,"",f_1,"\n")
    iter <- iter +1
  }
  
  f_3 <- f_3[-1]
  r <- c(0.05,0.05,0.05)
  iter <- 1
  while(1){
    v11 <- r[1]^2/(9-18*r[1]+10*r[1]^2)
    v21 <- r[2]^2/(9-18*r[2]+10*r[2]^2)
    v12 <- r[1]/(3-2*r[1])
    v22 <- r[2]/(3-2*r[2])
    
    #r1 <- (f_3[3]+f_3[5]+2*(f_3[2]+f_3[4]+f_3[6]+f_3[9])+2*v1*f_3[7]+(1+v2)*f_3[8])/2
    
    r1 <- sum(f_3[c(2,4,6,7,10,11,13,14,17,18,20,21,24,25,28,32,33,35,37,41,44,47,48,54,56,59)])+
      sum(f_3[c(3,8,9,15,16,22,23,26,27,34,40,45,46,55)])/2+sum(f_3[c(29,38,42,49,57)])*v11+
      sum(f_3[c(30,31,36,39,43,50,51,52,53,58)])*(1+v12)/2
    
    
    r2 <- sum(f_3[c(2,4,5,7,9,11,14,18,19,21,23,27,28,31,32,34,35,38,39,41,44,46,48,53,55,57)])+
      sum(f_3[c(3,8,10,12,13,20,29,30,36,37,40,45,47,56)])/2+sum(f_3[c(15,24,42,50,59)])*v21+
      sum(f_3[c(16,17,22,25,43,49,50,51,52,54,58)])*(1+v22)/2
    #if(r2>0.375){
    #  r2 <- 0.75-r2
    #}
    r3 <- sum(f_3[c(5:11)])+sum(f_3[c(19:25)])+sum(f_3[c(33:39)])+sum(f_3[c(55:59,52)])+
      sum(f_3[c(12:18)])+sum(f_3[c(26:32)])+sum(f_3[c(40,41,44,49,50,51)])/2+3*sum(f_3[c(45,46,47,48,53,54)])/4+
      sum(f_3[c(42,42,49)])*v11/2+sum(f_3[c(42,42,50)])*v12/2-2*f_3[42]*v11*v12+
                                           sum(f_3[c(43,43,49,51)])*v21/2-f_3[52]*v21/2-f_3[50]*v12*v21+
                                           sum(f_3[c(43,43,49,51)])*v22/2-f_3[52]*v22/2-f_3[49]*v11*v22-sum(f_3[c(43,43,51)])*v21*v22+f_3[52]*v21*v22
    #if(r3>0.375){
    #  r3 <- 0.75-r3
    #}
    cr <- c(r1,r2,r3)
    if(max(abs(cr-r))<1e-5||iter>50)
      break
    
    r <- cr
    iter <- iter +1
    
  }
  
  a <- sum(f_3[1:25])
  b <- sum(f_3[c(1,2,5:7,12:14,19:21,26:28,33:35,40,41,45:48,55,56)])
  c <- sum(f_3[c(1:11,26:39)])
  
  rnc <- rbind(r,abs(1-r))
  r <- apply(rnc, 2, min)
  coc <- (r[1]+r[2]-r[3])/(2*r[1]*r[2])
  
  ii <- f_1!=0
  LL <- sapp_L(sum(E1))-sum(sapp_L(E1[ii]))+sum(E1[ii]*log(f_1[ii]),nan.rm=T)
  res <- c(a,b,c,r,coc,LL)
  rm(P1f,P2f,GMZ,DLRZ,PFCZ1,EM1,EM2)
  gc()
  res
}








auto_three_scan <- function(ndat,phase){
  
  
  
  m1 <- ndat[1,]
  m2 <- ndat[2,]
  m3 <- ndat[3,]
  m11 <- m1[-c(1,2)]  
  m22 <- m2[-c(1,2)]
  m33 <- m2[-c(1,2)]
  mname <- rownames(ndat)
  P1M1 <- as.numeric(strsplit(as.character(m1[1]),"")[[1]])
  P1M2 <- as.numeric(strsplit(as.character(m2[1]),"")[[1]])
  P1M3 <- as.numeric(strsplit(as.character(m3[1]),"")[[1]])
  P2M1 <- as.numeric(strsplit(as.character(m1[2]),"")[[1]])
  P2M2 <- as.numeric(strsplit(as.character(m2[2]),"")[[1]])
  P2M3 <- as.numeric(strsplit(as.character(m3[2]),"")[[1]])
  pm <- matrix(c(P1M1,P1M2,P1M3,P2M1,P2M2,P2M3),nrow=6,byrow = T)
    
  if(all(P2M1==1)){
    pm[1,] <- P2M1
    pm[4,] <- P1M1
  }
    
  if(all(P2M2==1)){
    pm[2,] <- P2M2
    pm[5,] <- P1M2
  }
  
  if(all(P2M3==1)){
    pm[3,] <- P2M3
    pm[6,] <- P1M3
  }
  
    
  cat("M1: ",mname[1]," M2: ",mname[2]," M3: ",mname[3],"\n")
  ret <- three_rfz_phase(M1=m11,M2=m22,M3=m33,pmap=pm,index=phase)
  cat("",ret,"\n")
  
  return(ret)
}
