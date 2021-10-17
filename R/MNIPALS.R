#' Caculate the estimator on the MNIPALS method
#'
#' @param data is the orignal data set
#' @param data0 is the missing data set
#' @param real is to judge whether the data set is a real missing data set
#' @param example is to judge whether the data set is a simulation example.


#' 
#' @return XMNIPALS, MSEMNIPALS, MAEMNIPALS, REMNIPALS, GCVMNIPALS,timeMNIPALS
#' @export
#'

#' @examples 
#'  library(MASS)   
#'  etatol=0.9
#'  n=100;p=10;per=0.1
#'  mu=as.matrix(runif(p,0,10))
#'  sigma=as.matrix(runif(p,0,1))
#'  ro=as.matrix(c(runif(p,-1,1)))
#'  RO=ro%*%t(ro);diag(RO)=1
#'  Sigma=sigma%*%t(sigma)*RO
#'  X0=data=mvrnorm(n,mu,Sigma)
#'  m=round(per*n*p,digits=0)
#'  mr=sample(1:(n*p),m,replace=FALSE)
#'  X0[mr]=NA;data0=X0
#'  MNIPALS(data=data,data0=data0,real=FALSE,example=FALSE)


#the MNIPALS method 
MNIPALS=function(data=0,data0,real=TRUE,example=FALSE)
#It defaults that the data set is a real data set
{#1
  if(real||example){#2
    etatol=0.7
  }else{#2
    etatol=0.9
  }#2
  lll=0
  time=system.time(#2
    while(lll==0){#3
      X0=data0
      n=nrow(X0);p=ncol(X0)
      mr=which(is.na(X0)==TRUE)
      m=nrow(as.matrix(mr))
      cm0=colMeans(X0,na.rm=T)
      ina=as.matrix(mr%%n)
      jna=as.matrix(floor((mr+n-1)/n))
      data0[is.na(data0)]=cm0[ceiling(which(is.na(X0))/n)]	
      X=as.matrix(data0)
      Z0=Z=scale(X,center=TRUE,scale=FALSE)
      tol=1e-5;nb=100;niter=0;d=1
      t0=matrix(0,n,1);a0=matrix(0,p,1)
      t1=matrix(Z[,which.max(diag(t(Z)%*%Z))],,1)
      while((d>=tol) & (niter<=nb)){#4
        niter=niter+1
        thethat=matrix(t(ginv(t(t1)%*%t1)%*%(t(t1)%*%Z)),p,1)
        a1hat=matrix(thethat/sqrt(sum(thethat*thethat)),p,1)
        t1hat=matrix(t(ginv(t(a1hat)%*%a1hat)%*%(t(a1hat)%*%t(Z))),n,1)
        d=sqrt(t(t1hat-t1)%*%(t1hat-t1))
        t1=t1hat
      }#4
      for( i in 1:n){#4
        M=is.na(X0[i,])
        job=which(M==FALSE);jna=which(M==TRUE)
        piob=nrow(as.matrix(job));pina=nrow(as.matrix(jna))
        while((piob>0)&(pina>0)){#5
          Qi=matrix(0,p,p)
          for( u in 1:piob){#6
            Qi[job[u],u]=1
          }#6
          for( v in 1:pina){#6
            Qi[jna[v],v+piob]=1
          }#6
          zi=Z0[i,]
          zQi=zi%*%Qi
          ZQi=Z0%*%Qi
          a1Qi=t(t(a1hat)%*%Qi)
          ziob=t(as.matrix(zQi[,1:piob]))
          zina=t(as.matrix(zQi[,piob+(1:pina)]))
          Ziob=matrix(ZQi[,1:piob],n,piob)
          Zina=matrix(ZQi[,piob+(1:pina)],n,pina)
          a1iob=matrix(t(a1Qi)[,1:piob],piob,1)
          a1ina=matrix(t(a1Qi)[,piob+(1:pina)],pina,1)
          zinahat=t1hat[i,]*t(a1ina)
          ZQi[i,piob+(1:pina)]=zinahat
          Zi=ZQi%*%t(Qi)
          Z=Z0=Zi
          pina=0
        }#5
      }#4
      XMNIPALS=Xnew=Z+matrix(1,n,p)%*%diag(cm0)
      lll=1
    }#3
  )#2
  if(real){#2
    MSEMNIPALS= MAEMNIPALS= REMNIPALS='NULL'
  }else{#2
    MSEMNIPALS=(1/m)*t(Xnew[mr]-data[mr])%*%(Xnew[mr]-data[mr])
    MAEMNIPALS=(1/m)*sum(abs(Xnew[mr]-data[mr]))	
    REMNIPALS=(sum(abs(data[mr]-Xnew[mr])))/(sum(data[mr]))
  }#2
  lambdaMNIPALS=svd(cor(XMNIPALS))$d
  lMNIPALS=lambdaMNIPALS/sum(lambdaMNIPALS);J=rep(lMNIPALS,times=p);dim(J)=c(p,p)
  upper.tri(J,diag=T);J[lower.tri(J)]=0;J;dim(J)=c(p,p)
  etaMNIPALS=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
  wwMNIPALS=which(etaMNIPALS>=etatol);kMNIPALS=wwMNIPALS[1] 
  lambdaMNIPALSpk=lambdaMNIPALS[(kMNIPALS+1):p]
  GCVMNIPALS=sum(lambdaMNIPALSpk)*p/(p-kMNIPALS)^2
  return(list(XMNIPALS=XMNIPALS,MSEMNIPALS=MSEMNIPALS,MAEMNIPALS=MAEMNIPALS,REMNIPALS=REMNIPALS,GCVMNIPALS=GCVMNIPALS,timeMNIPALS=time))
}#1

