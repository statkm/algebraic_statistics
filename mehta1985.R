##networkcmh: c_sを計算する
##exactcmh: networkcmhを利用して正確検定


networkcmh<-function(xxx){
  #Mehta1985 (1,1)和を検定統計量とした場合
  K<-dim(xxx)[3]
  #(1,1)セルの上限下限
  yy <-apply(cbind(rep(0,K),xxx[1,1,]-xxx[2,2,]),1,max)
  aa <- cumsum(yy)
  zz <- apply(cbind(xxx[1,1,]+xxx[1,2,],xxx[1,1,]+xxx[2,1,]),1,min)
  bb <- cumsum(zz)
  
  #周辺和
  tt <- apply(xxx, c(2,3),sum)[1,]
  mm <- apply(xxx, c(1,3),sum)[1,]
  nn <- apply(xxx, c(1,3),sum)[2,]
  
  cc <- vector("list",K+1)
  cc[[1]] <- 0
  
  
  for( i in 0:(K-1) ){
    if(i==0){
      for(j in 1:(zz[i+1]-yy[i+1]+1) )
        cc[[i+2]][j]<- log(choose(mm[i+1], yy[i+1]+j-1)*choose(nn[i+1], tt[i+1]-(yy[i+1]+j-1)))
    }else{
      cc[[i+2]]<-rep(-Inf,bb[i+1]-aa[i+1]+1); 
      for(k in 1:( bb[i]-aa[i]+1 )){
        for(j in 1:( zz[i+1]-yy[i+1]+1 )){
          delta <- max( cc[[i+2]][k+j-1], 
                        cc[[i+1]][k]+log(choose(mm[i+1], yy[i+1]+j-1)*choose(nn[i+1], tt[i+1]-(yy[i+1]+j-1))) )
          cc[[i+2]][k+j-1]<- 
            log(exp(cc[[i+2]][k+j-1]-delta) + 
                  exp(cc[[i+1]][k]+log(choose(mm[i+1], yy[i+1]+j-1)*choose(nn[i+1], tt[i+1]-(yy[i+1]+j-1)))-delta) ) + delta
        }
      }
    }
  } #ccの計算
  
  return( list(logged=cc[[K+1]],aa=aa,bb=bb ) )
}






exactcmh <- function(xxx,alternative="two.sided", correct =T){
  K<-dim(xxx)[3];
  obs <- sum(xxx[1,1,]);
  
  yy <-apply(cbind(rep(0,K),xxx[1,1,]-xxx[2,2,]),1,max)
  aa <- cumsum(yy)
  zz <- apply(cbind(xxx[1,1,]+xxx[1,2,],xxx[1,1,]+xxx[2,1,]),1,min)
  bb <- cumsum(zz)
  
#  logged<- networkcmh(xxx)$logged
#  delta <-max( networkcmh(xxx)$logged )
#  logden <- log(sum(exp(logged-delta)))+delta #分母
#  den <- exp(logden)
#  cc<-exp( logged-logden )
  
  cc <- exp(networkcmh(xxx)$logged)
  den <- sum(cc)
  
  cri<-obs-aa[K]+1 #添字としてのobsの位置
  switch(alternative,
         "greater" = {
           num <- sum( tail(cc,bb[K]-aa[K]+1-cri) ) + cc[cri]*(1-correct/2)
         }, 
         "lower" = {
           num <- sum( head(cc,cri-1) ) + cc[cri]*(1-correct/2)
         }, 
         { #両側
           num <- sum( tail(cc,bb[K]-aa[K]+1-cri) ) + cc[cri]*(1-correct/2)
           if(num<den/2) num <- 2*num
           else num <- 2*(den-num+cc[cri]*(1-correct) )
         }
  )
  return(list(p.value=min(num/den,1), obs=obs,data=xxx, dteststat=cc/den,obsmin=aa[K],obsmax=bb[K]))
}