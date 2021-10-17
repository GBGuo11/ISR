#' Caculate the estimator on the MMLPCA method
#'
#' @param data is the orignal data set
#' @param data0 is the missing data set
#' @param real is to judge whether the data set is a real missing data set
#' @param example is to judge whether the data set is a simulation example.


#' 
#' @return XMMLPCA, MSEMMLPCA, MAEMMLPCA, REMMLPCA, GCVMMLPCA,timeMMLPCA
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
#'  MMLPCA(data=data,data0=data0,real=FALSE,example=FALSE)


##the MMLPCA method 
MMLPCA=function(data=0,data0,real=TRUE,example=FALSE)
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
      tol=1e-5;nb=10;niter=0;d=1;SS=1
      R=cor(Z)
      lambda=svd(R)$d;l=lambda/sum(lambda);J=rep(l,times=p);dim(J)=c(p,p)
      upper.tri(J,diag=T);J[lower.tri(J)]=0;eta=matrix(colSums(J),nrow=1,ncol=p,byrow=FALSE)
      ww=which(eta>=etatol);k=ww[1]
      while((SS>=tol)&(niter<=nb)){#4
        niter=niter+1
        Zold=Z
        R=cor(Z)
        A=svd(Z)$v
        Ak=matrix(A[,1:k],p,k)
        for (i in 1:n) {#5
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
            zinahat=t(Aina%*%t(Aiob)%*%t(ziob))
            ZQi[i,piob+(1:pina)]=zinahat
            Zi=ZQi%*%t(Qi)
            Z=Zi
            pina=0
          }#6
        }#5
        Zrow=Znew=Z
        S1=sum((data0[mr]-Zrow[mr])^2)
        B=svd(Z)$u
        Bk=matrix(B[,1:k],n,k)
        for (j in 1:p) {#5
          N=is.na(X0[,j])
          iob=which(N==FALSE);ina=which(N==TRUE)
          njob=nrow(as.matrix(iob));njna=nrow(as.matrix(ina))
          while((njob>0)&(njna>0)){#6
            Qj=matrix(0,n,n)
            for(u in 1:njob){#7
              Qj[u,iob[u]]=1
            }#7
            for(v in 1:njna){#7
              Qj[v+njob,ina[v]]=1
            }#7
            zj=Z[,j]
            zQj=Qj%*%zj
            ZQj=Qj%*%Z
            BQj=t(t(Bk)%*%Qj)
            zjob=matrix(zQj[1:njob,],njob,1)
            zjna=matrix(zQj[njob+(1:njna),],njna,1)
            Zjob=matrix(ZQj[1:njob,],njob,p,byrow=FALSE)
            Zjna=matrix(ZQj[njob+(1:njna),],njna,p,byrow=FALSE)
            Bjob=matrix(BQj[1:njob,],njob,k,byrow=FALSE)
            Bjna=matrix(t(BQj)[,njob+(1:njna)],njna,k,byrow=FALSE)
            zjnahat=Bjna%*%t(Bjob)%*%zjob
            ZQj[njob+(1:njna),j]=zjnahat
            Zj=t(Qj)%*%ZQj
            Z=Zj
            njna=0
          }#6
        }#5
        ZMMLPCA=Zcol=Znew=Z
        S2=sum((data0[mr]-Zcol[mr])^2)
        SS=abs(S2-S1)/S2
      }#4
      XMMLPCA=Xnew=Znew+matrix(rep(1,n*p),ncol=p)%*%diag(cm0)
      lll=1
    }#3
  )#2
  if(real){#2
    MSEMMLPCA= MAEMMLPCA= REMMLPCA='NULL'
  }else{#2
    MSEMMLPCA=(1/m)*t(Xnew[mr]-data[mr])%*%(Xnew[mr]-data[mr])
    MAEMMLPCA=(1/m)*sum(abs(data[mr]-Xnew[mr]))  
    REMMLPCA=(sum(abs(data[mr]-Xnew[mr])))/(sum(data[mr]))
  }#2
  lambdaMMLPCA=svd(cor(XMMLPCA))$d
  lMMLPCA=lambdaMMLPCA/sum(lambdaMMLPCA);J=rep(lMMLPCA,times=p);dim(J)=c(p,p)
  upper.tri(J,diag=T);J[lower.tri(J)]=0;dim(J)=c(p,p)
  etaMMLPCA=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
  wwMMLPCA=which(etaMMLPCA>=etatol);kMMLPCA=wwMMLPCA[1] 
  lambdaMMLPCApk=lambdaMMLPCA[(kMMLPCA+1):p]
  GCVMMLPCA=sum(lambdaMMLPCApk)*p/(p-kMMLPCA)^2
  return(list(XMMLPCA=XMMLPCA,MSEMMLPCA=MSEMMLPCA,MAEMMLPCA=MAEMMLPCA,REMMLPCA=REMMLPCA,GCVMMLPCA=GCVMMLPCA,timeMMLPCA=time))
}#1
