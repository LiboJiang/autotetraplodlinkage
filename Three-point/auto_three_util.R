



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



three_phase <- function(M1,M2,M3){
  
  
  mc <- expand.grid(1:4,1:4,1:4,1:4)
  
  mc1 <- apply(mc,1,function(x){
    if(max(table(x))==1)
      return(1)
    else{
      return(0)
    }
  })
  
  mc2 <- as.matrix(mc[as.logical(mc1),])
  
  
  mc <- rbind(M1,M2,M3)
  #mc <- sort_g(mc)
  mc12 <- c(mc[1,]*10+mc[2,])
  mc23 <- c(mc[2,]*10+mc[3,])
  mc13 <- c(mc[1,]*10+mc[3,])
  
  
  M11 <- sort(M1)
  m1_l <- cbind(matrix(rep(M11,2),nrow=2,byrow =T),matrix(M11[combn(4,2)],2))
  m1_l1 <- sort_g(m1_l)
  m1_g <- m1_l1[1,]*10+m1_l1[2,]
  
  M22 <- sort(M2)
  m2_l <- cbind(matrix(rep(M22,2),nrow=2,byrow =T),matrix(M22[combn(4,2)],2))
  m2_l1 <- sort_g(m2_l)
  m2_g <- m2_l1[1,]*10+m2_l1[2,]
  
  M33 <- sort(M3)
  m3_l <- cbind(matrix(rep(M33,2),nrow=2,byrow =T),matrix(M33[combn(4,2)],2))
  m3_l1 <- sort_g(m3_l)
  m3_g <- m3_l1[1,]*10+m3_l1[2,]
  
  
  
  nm <- array(0,dim=c(length(m1_g),length(m2_g),length(m3_g)))
  
  
  for(ii in 1:dim(mc2)[1]){
    mh1 <- M1[mc2[ii,]]
    mh2 <- M2[mc2[ii,]]
    mh3 <- M3[mc2[ii,]]
    #cat(ii,"\n")
    for(i in 1:length(m1_g)){
      for(j in 1:length(m2_g)){
        for(k in 1:length(m3_g)){
          
          tmp <- c(m1_l1[1,i],m2_l1[1,j],m3_l1[1,k],m1_l1[2,i],m2_l1[2,j],m3_l1[2,k])
          kl112 <- m1_l1[1,i]*10+m2_l1[1,j]
          kl123 <- m2_l1[1,j]*10+m3_l1[1,k]
          kl113 <- m1_l1[1,i]*10+m3_l1[1,k]
          
          kl212 <- m1_l1[2,i]*10+m2_l1[2,j]
          kl223 <- m2_l1[2,j]*10+m3_l1[2,k]
          kl213 <- m1_l1[2,i]*10+m3_l1[2,k]
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
             (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
             ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
              &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))){
            nm[j,i,k]=1
            
          }
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
             (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
             ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
              &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))){
            nm[j,i,k]=2
            
          }
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=3
            
          }
          
          
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
             (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
             (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
               &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))||
              ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
               &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))){
            nm[j,i,k]=4
            
          }
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
             (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
             ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
              &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2]))){
            nm[j,i,k]=5
          }
          
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
             (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
             ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
              &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2]))){
            nm[j,i,k]=6
            
          }
          
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
             (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
             ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
              &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2]))){
            nm[j,i,k]=7
            
          }
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=8
            
          }
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=9
          }
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=10
            
          }
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
             (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
             (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
               &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2]))||
              ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
               &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[2])))){
            nm[j,i,k]=11
            
          }
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=12
            
          }
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=13
            
          }
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=14
            
          }
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
               ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))||
              
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=15
          }
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))||
              
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=16
            
          }
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))||
              
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=17
            
          }
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=18
          }
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
              (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
              (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3]))||
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=19
            
          }
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=20
          }
          
          
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
              (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
              (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3]))||
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=21
            
          }
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2]))))||
              
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3])))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=22
            #cat("",i,j,k," ",m1_l1[1,i],m2_l1[1,j],m3_l1[1,k],"-",m1_l1[2,i],m2_l1[2,j],m3_l1[2,k],"\n")
          }
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=23
            
          }
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                  &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3]))||
                 ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                  &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2])))))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3]))||
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2])))))){
            nm[j,i,k]=24
            
          }
          
          
          if((m1_l1[1,i]==m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                  &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3]))||
                 ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                  &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2])))))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[3])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=25
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=26
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=27
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            
            nm[j,i,k]=28
            
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=29
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=30
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=31
            
          }
          
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))
              ||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            
            nm[j,i,k]=32
            
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))))){
            
            nm[j,i,k]=33
            
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13)))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            
            nm[j,i,k]=34
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             (((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))))){
            
            nm[j,i,k]=35
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            nm[j,i,k]=36
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))){
            
            nm[j,i,k]=37
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                (((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                  &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))||
                 ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                  &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))))){
            nm[j,i,k]=38
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]==m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[3])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[3])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[1])))))){
            nm[j,i,k]=39
            
          }
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1])))||
              
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=40
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                  &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))||
                 ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                  &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1])))))
              ||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))))){
            nm[j,i,k]=41
            #cat("",i,j,k," ",m1_l1[1,i],m2_l1[1,j],m3_l1[1,k],"-",m1_l1[2,i],m2_l1[2,j],m3_l1[2,k],"\n")
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                  &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))||
                 ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                  &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1])))))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))||
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))||
              
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))))||
              
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))))){
            nm[j,i,k]=42
            
          }
          
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))||
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))))||
              
              
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=43
            
            
          }
          
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                  &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))||
                 ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                  &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1])))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13)))&&
              (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[2])&&(m3_l1[2,k]==mh3[1]))||
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[2]))))){
            nm[j,i,k]=44
            
          }
          
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))){
            nm[j,i,k]=45
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))){
            nm[j,i,k]=46
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))){
            nm[j,i,k]=47
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))))){
            nm[j,i,k]=48
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))){
            nm[j,i,k]=49
            
            #cat("",i,j,k," ",m1_l1[1,i],m2_l1[1,j],m3_l1[1,k],"-",m1_l1[2,i],m2_l1[2,j],m3_l1[2,k],"\n")
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
              
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1]))||
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))))){
            nm[j,i,k]=50
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[4])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))){
            nm[j,i,k]=51
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13)))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
              
              ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[2])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
              ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))){
            nm[j,i,k]=52
            
          }
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (
               (((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                 (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                 ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[4])
                  &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
               
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
               ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
               
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
               
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
               ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))
               
             )){
            nm[j,i,k]=53
            
            
          }
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (
               (((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                 (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                 ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                  &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3]))))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
               
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
               ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
               
               ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))||
               
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[1])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[1])&&(m3_l1[2,k]==mh3[3])))
               
             )){
            nm[j,i,k]=54
            
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13)))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[1])
                &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))
             )){
            nm[j,i,k]=55
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]==m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             ((((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13)))&&
               ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
                &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4])))||
              ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
               (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
               (m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[3])
               &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3]))
             )){
            nm[j,i,k]=56
            
          }
          
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (
               (((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                 (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                 ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                  &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4]))))||
               ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))||
               
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[2])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[2])
                 &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))
               
             )){
            
            nm[j,i,k]=57
            
          }
          
          
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (
               (((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                 (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                 ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                  &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4]))))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))||
               
               ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))||
               
               ((is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))||
               ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4])))||
               
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[1])
                 &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[1])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4])))
               
             )){
            
            nm[j,i,k]=58
            
          }
          if((m1_l1[1,i]!=m1_l1[2,i])&&(m2_l1[1,j]!=m2_l1[2,j])&&(m3_l1[1,k]!=m3_l1[2,k])&&
             (
               (((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                 (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13)))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4])))||
               ((!is.element(kl112,mc12))&&(is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))||
               
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[1])&&(m1_l1[2,i]==mh1[2])&&(m2_l1[1,j]==mh2[4])&&(m2_l1[2,j]==mh2[3])
                 &&(m3_l1[1,k]==mh3[3])&&(m3_l1[2,k]==mh3[4])))||
               ((!is.element(kl112,mc12))&&(!is.element(kl123,mc23))&&(!is.element(kl113,mc13))&&
                (!is.element(kl212,mc12))&&(!is.element(kl223,mc23))&&(!is.element(kl213,mc13))&&
                ((m1_l1[1,i]==mh1[2])&&(m1_l1[2,i]==mh1[1])&&(m2_l1[1,j]==mh2[3])&&(m2_l1[2,j]==mh2[4])
                 &&(m3_l1[1,k]==mh3[4])&&(m3_l1[2,k]==mh3[3])))
               
             )){
            nm[j,i,k] <- 59
            
          }
          
          #if(nm[j,i,k]==0){
          #  cat(j,i,k,"\n")
          #  break
          #}
        }
      }
    }
  }
  
  return(melt_m(nm))
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
  
  nM <- matrix(NA,60,60)
  ni <- 0:59
  for(i in 1:60){
    for(j in 1:60){
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
  for(ii in 0:59){
    
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
  
  ff1 <- f[t1s+1]
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
  ff11[t1s+1] <- ff1/st1
  ff11[-(t1s+1)] <- 0
  
  t2s <- sort(unique(c(P2f)))
  
  ff2 <- f[t2s+1]
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
  ff22[t2s+1] <- ff2/st2
  ff22[-(t2s+1)] <- 0
  
  f2 <- c(t(kronecker(ff11,t(ff22))))
  
  tmp2 <- kronecker(ind1,ind)
  return(list(HP=tmp2,f2=f2,fp1=ff1,fp2=ff2))
  
  
  
}


sapp_L <-function(n){
  log(sqrt(2*pi*n))+n*log(n/exp(1))
}

Haldan <- function(rf){
  
  
  -0.5 * log(1 - 2 * rf)
  
  
}