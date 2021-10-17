#' Caculate the estimator on the MRPCA method
#'
#' @param data is the orignal data set
#' @param data0 is the missing data set
#' @param real is to judge whether the data set is a real missing data set
#' @param example is to judge whether the data set is a simulation example.


#' 
#' @return XMRPCA, MSEMRPCA, MAEMRPCA, REMRPCA, GCVMRPCA,timeMRPCA
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
#'  MRPCA(data=data,data0=data0,real=FALSE,example=FALSE)

#the MRPCA method
MRPCA=function(data=0,data0,real=TRUE,example=FALSE)
#It defaults that the data set is a real data set
{#1
  if (real||example){#2
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
      Z=scale(X,center=TRUE,scale=FALSE)
      niter=0;d=1;tol=1e-5;nb=10
      while((d>=tol) & (niter<=nb)){#4
        niter=niter+1
        Xold=X
        Zold=Z
        R=cor(Z)
        lambda=svd(R)$d
        l=lambda/sum(lambda)
        J=rep(l,times=p);dim(J)=c(p,p)
        upper.tri(J,diag=T);J[lower.tri(J)]=0
        eta=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
        ww=which(eta>=etatol)
        k=ww[1] 
        Lambda=svd(Z)$d
        A=svd(Z)$v
        B=svd(Z)$u
        Lambdak=diag(sqrt(lambda[1:k]),k,k)
        Ak=matrix(A[,1:k],p,k);Bk=matrix(B[,1:k],n,k)
        Lambdapk=diag(sqrt(lambda[(k+1):p]),p-k,p-k)
        sigma2hat=sum(diag(Lambdapk%*%Lambdapk))/(p-k)
        for( i in 1:n){#5
          M=is.na(X0[i,])
          job=which(M==FALSE);jna=which(M==TRUE)
          piob=nrow(as.matrix(job));pina=nrow(as.matrix(jna))
          while((piob>0)&(pina>0)){#6
            Qi=matrix(0,p,p)
            for( u in 1:piob){#7
              Qi[job[u],u]=1
            }#7
            for( v in 1:pina){#7
              Qi[jna[v],v+piob]=1
            }#7
            zi=Z[i,]
            zQi=zi%*%Qi
            ZQi=Z%*%Qi
            AQi=t(t(Ak)%*%Qi)
            ziob=matrix(zQi[,1:piob],1,piob)
            zina=matrix(zQi[,piob+(1:pina)],1,pina)
            Ziob=matrix(ZQi[,1:piob],n,piob,byrow=FALSE)
            Zina=matrix(ZQi[,piob+(1:pina)],n,pina,byrow=FALSE)
            Aiob=matrix(AQi[1:piob,],piob,k,byrow=FALSE)
            Aina=matrix(AQi[piob+(1:pina),],pina,k,byrow=FALSE)
            Cihat=n^(-1/2)*Aina%*%(Lambdak%*%Lambdak-sigma2hat*diag(k))^(1/2)
            tihat=n^(1/2)*solve(Lambdak)%*%(Lambdak%*%Lambdak-sigma2hat*diag(k))^(1/2)%*%Bk[i,]
            zinahat=Cihat%*%tihat
            ZQi[i,piob+(1:pina)]=zinahat
            Zi=ZQi%*%t(Qi)
            Z=Zi
            pina=0
          }#6
        }#5
        ZMRPCA=Znew=Z
        d=sqrt(sum(diag((t(Zold-Znew)%*%(Zold-Znew)))))
      }#4
      XMRPCA=Xnew=Znew+matrix(rep(1,n*p),ncol=p)%*%diag(cm0)
      lll=1
    }#3
  )#2
  if(real){#2
    MSEMRPCA= MAEMRPCA= REMRPCA='NULL'
  }else{#2
    MSEMRPCA=(1/m)*t(Xnew[mr]-data[mr])%*%(Xnew[mr]-data[mr])
    MAEMRPCA=(1/m)*sum(abs(Xnew[mr]-data[mr]))	
    REMRPCA=(sum(abs(data[mr]-Xnew[mr])))/(sum(data[mr]))
  }#2
  lambdaMRPCA=svd(cor(XMRPCA))$d
  lMRPCA=lambdaMRPCA/sum(lambdaMRPCA);J=rep(lMRPCA,times=p);dim(J)=c(p,p)
  upper.tri(J,diag=T);J[lower.tri(J)]=0;dim(J)=c(p,p)
  etaMRPCA=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
  wwMRPCA=which(etaMRPCA>=etatol);kMRPCA=wwMRPCA[1] 
  lambdaMRPCApk=lambdaMRPCA[(kMRPCA+1):p]
  GCVMRPCA=sum(lambdaMRPCApk)*p/(p-kMRPCA)^2
return(list(XMRPCA=XMRPCA,MSEMRPCA=MSEMRPCA,MAEMRPCA=MAEMRPCA,REMRPCA=REMRPCA,GCVMRPCA=GCVMRPCA,timeMRPCA=time))
}#1


