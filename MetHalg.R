#Metropolis-Hastings Algorithm
#return a sequence of pvalues


#20171203 burn-in added
#20171210 両側検定, mcmcのmethodの選択
#20171217 半数補正
#20171219 棄却域修正
#20171222 棄却域修正
#20180101 MCMCにskipを導入
#20180110 採択確率の計算を対数を用いて計算
#20180206 correctを消去
#20180212 mid p-value実装
#20180214 mid p-value修正


#非心超幾何確率函数比例部分
pnoncenthyp<-function(xxx,odds){ #xxxは2x2分割表
  return( choose(sum(xxx[1,]),xxx[1,1]) * choose(sum(xxx[2,]),xxx[2,1]) * odds^(xxx[1,1]) )
}

# logged
plognoncenthyp<-function(xxx,odds){ #xxxは2x2分割表
  return( log(choose(sum(xxx[1,]),xxx[1,1])) + log(choose(sum(xxx[2,]),xxx[2,1])) + xxx[1,1]*log(odds) )
}



Methalg<-function(N,xxx, odds,nburnin, alternative="two.sided", midp=F,
                  method="random", mcmcstep=1){
  K<-dim(xxx)[3]
  obs<-sum(xxx[1,1,]);
  min_str <-apply(cbind(rep(0,K),xxx[1,1,]-xxx[2,2,]),1,max)
  obs_min <- sum(min_str)
  max_str <- apply(cbind(xxx[1,1,]+xxx[1,2,],xxx[1,1,]+xxx[2,1,]),1,min)
  obs_max <- sum(max_str)
  cnt<-0;
  sig<-0;
  sig2<-0; #棄却域の端
  pval<-numeric(N);
  teststat<-numeric(N*mcmcstep+nburnin)
  xcur<-xxx
  
  #reject region
  if(alternative=="lower"){
    rej_max <-Inf; rej_min <- obs
  }else{
    rej_max <-obs; rej_min <- -Inf
  }
  #rej reg fin
  
  
  while(cnt<nburnin+N*mcmcstep){ #MCMCloop
    #推移しうる状態
    xtemp<-xcur
    #実際に推移する状態
    xnext<-xtemp
    
    if(method=="each"){
      eps<-2*rbinom(n=K,size=1,p=0.5)-1
      for(rmb in 1:K){
        xtemp[1,1,rmb]<-xtemp[1,1,rmb]+eps[rmb]
        xtemp[1,2,rmb]<-xtemp[1,2,rmb]-eps[rmb]
        xtemp[2,1,rmb]<-xtemp[2,1,rmb]-eps[rmb]
        xtemp[2,2,rmb]<-xtemp[2,2,rmb]+eps[rmb]
        #MH algorithm
        if(xtemp[1,1,rmb]<0|xtemp[1,2,rmb]<0|xtemp[2,1,rmb]<0|xtemp[2,2,rmb]<0){
          xnext[,,rmb]<-xcur[,,rmb]
        }else{
          if( runif(1) < exp( plognoncenthyp(xtemp[,,rmb],odds)-plognoncenthyp(xcur[,,rmb],odds) ) ) 
            xnext[,,rmb]<-xtemp[,,rmb]
          else
            xnext[,,rmb]<-xcur[,,rmb]
        }
      }
    } #each
    
    
    if(method=="random"){
      eps<-2*rbinom(1,1,p=0.5)-1
      mb<-sample(1:K,1)
        xtemp[1,1,mb]<-xtemp[1,1,mb]+eps; xtemp[1,2,mb]<-xtemp[1,2,mb]-eps
        xtemp[2,1,mb]<-xtemp[2,1,mb]-eps; xtemp[2,2,mb]<-xtemp[2,2,mb]+eps

        #MH algorithm
        if( (xtemp[1,1,mb]<0)||(xtemp[1,2,mb]<0)||(xtemp[2,1,mb]<0)||(xtemp[2,2,mb]<0) ){
          xnext[,,mb]<-xcur[,,mb]
        }else{
          if( runif(1) < exp(plognoncenthyp(xtemp[,,mb],odds)-plognoncenthyp(xcur[,,mb],odds) ) ) 
            xnext[,,mb]<-xtemp[,,mb]
          else  
            xnext[,,mb]<- xcur[,,mb]
        }
    } #random 
    
    
    
    if(method=="diag"){
#      eps<-c(0,0)
#      while(!sum(eps^2)) eps<-sample(1:3,2,replace = T)-2
      eps<-sample(1:3,2,replace = T)-2
      mb<-sample(1:K,2,replace = F)
      xtemp[,,mb[1]]<-xtemp[,,mb[1]]+matrix(c(eps[1],-eps[1],-eps[1],eps[1]),ncol=2)
      xtemp[,,mb[2]]<-xtemp[,,mb[2]]+matrix(c(eps[2],-eps[2],-eps[2],eps[2]),ncol=2)      
      
      #MH algorithm
        if( (xtemp[1,1,mb[1]]<0)||(xtemp[1,2,mb[1]]<0)||(xtemp[2,1,mb[1]]<0)||(xtemp[2,2,mb[1]]<0)
            ||(xtemp[1,1,mb[2]]<0)||(xtemp[1,2,mb[2]]<0)||(xtemp[2,1,mb[2]]<0)||(xtemp[2,2,mb[2]]<0)  ){
          xnext[,,mb[1]]<-xcur[,,mb[1]]; xnext[,,mb[2]]<-xcur[,,mb[2]]
        }else{
          if( runif(1) < exp(plognoncenthyp(xtemp[,,mb[1]],odds) + plognoncenthyp(xtemp[,,mb[2]],odds)
                             -plognoncenthyp(xcur[,,mb[1]],odds) - plognoncenthyp(xcur[,,mb[2]],odds) ) 
          ){
            xnext[,,mb[1]]<-xtemp[,,mb[1]]; xnext[,,mb[2]]<-xtemp[,,mb[2]]
          } 
          else{
            xnext[,,mb[1]]<- xcur[,,mb[1]]; xnext[,,mb[2]]<- xcur[,,mb[2]]
          }  
        }
    } #diag 
    
    
    
    
    ##############transition###############
    
    xcur<-xnext

    if(cnt>=nburnin){ #burnin終了
      
      if( (cnt-nburnin)%%mcmcstep==0 ){
        if( ( sum(xnext[1,1,]) >= rej_max) || ( sum(xnext[1,1,]) <= rej_min ) ){ #rej. region
          if( ( sum(xnext[1,1,]) == rej_max ) || ( sum(xnext[1,1,]) == rej_min ) ){ #boundary
            sig2<-sig2+rbinom(1,1,prob = 1-midp/2) #mid p-value
            sig<-sig+rbinom(1,1,prob = 1-midp/2)
          }else{ #not on boundary
            sig<-sig+1
          }
        } #rej. region fin
        pval[ (cnt-nburnin)%/%mcmcstep +1 ]<-sig/( (cnt-nburnin)%/%mcmcstep +1 )
 
      } #mcmcstep
      
    } #reject check(after burnin)
    cnt<-cnt+1           
    teststat[cnt]<- sum(xnext[1,1,])
  
  }  #MCMCloop
  
  
  pvalest <- pval[N]
  if(alternative=="two.sided"){
     pvalest <- min( 2*(1-pval[N]+sig2/N*(1-midp)),2*pval[N], 1)
  }
  
  
  return( list(p.value=pvalest,pval=pval,teststat=teststat
               , obs=list(obs_min=obs_min,obs_max=obs_max,rej_min=rej_min,rej_max=rej_max) , data=xxx ) )

  
}
