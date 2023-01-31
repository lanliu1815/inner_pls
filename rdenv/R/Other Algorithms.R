#' SIMPLS for envelope estimation
#'
#' Specifically, the algorithm in Section 6.5.1, Cook's Book.
#'
#' @return
#' \item{Gammahat}{Estimator for the envelope subspace.}
#' \item{Wu}{The original \eqn{W_u} matrix.}
SIMPLS_gamma<-function(M, U, u){

  r<-dim(M)[1]
  if(u>r){
    stop("u is larger than r...")
  }

  Wu<-matrix(NA, ncol=u, nrow=r)  # This is just the matrix to output.

  #---k=0----

  Wu[,1]<-eigen(U, symmetric=T, only.values=F)$vectors[,1]  # The eigenvector corresponding to the largest eigenvalue of U

  #---1<=k<=u-1---

  if(u>1){
    for(i in 1:(u-1)){

      QEcal<-diag(r)-Pr((M%*%Wu[,1:i]), rank=i)

      Wu[,(i+1)]<-eigen(crossprod(QEcal, (U%*%QEcal)), symmetric=T, only.values=F)$vectors[,1]  # 'Largest' eigenvector of QEcal_U_QEcal
    }
  }

  Gammahat<-qr.Q(qr(Wu))

  return(list(Gammahat=Gammahat, Wu=Wu))

}


#' Improvement of \code{\link{SIMPLS_gamma}}.
#'
#' Use eigenvectors in calculating \eqn{Q_{Ecal}}. But it seems this does not make
#' a different at all.
SIMPLS_gamma2<-function(M, U, u){

  r<-dim(M)[1]
  if(u>r){
    stop("u is larger than r...")
  }

  Wu<-matrix(NA, ncol=u, nrow=r)  # This is just the matrix to output.

  #---k=0----

  Wu[,1]<-eigen(U, symmetric=T, only.values=F)$vectors[,1]  # The eigenvector corresponding to the largest eigenvalue of U

  #---1<=k<=u-1---

  if(u>1){
    for(i in 1:(u-1)){

      eig.MWuM<-eigen(tcrossprod(M%*%Wu[,1:i],NULL), symmetric=T)
      QEcal<-tcrossprod(eig.MWuM$vectors[,(i+1):r], NULL)

      Wu[,(i+1)]<-eigen(crossprod(QEcal, (U%*%QEcal)), symmetric=T, only.values=F)$vectors[,1]  # 'Largest' eigenvector of QEcal_U_QEcal
    }
  }

  Gammahat<-qr.Q(qr(Wu))

  return(list(Gammahat=Gammahat, Wu=Wu))

}

#' SIMPLS for envelope regression estimation
#'
#' Deals with sample estimation. Invokes \code{\link{SIMPLS_gamma}} to calculate \eqn{\Gamma}.
#'
#' @return
#' \item{beta.PLS}
#' \item{Gamma.PLS}
#' \item{mu.PLS}
SIMPLS_env<-function(X, Y, u){

  n<-nrow(X)
  SX<-stats::cov(X)*(n-1)/n
  SXY<-stats::cov(X, Y)*(n-1)/n
  resu_SIMPLS<-SIMPLS_gamma(M=SX, U=tcrossprod(SXY), u=u)
  Gam.PLS<-resu_SIMPLS$Gammahat
  beta.PLS<-Gam.PLS%*%chol2inv(chol(t(Gam.PLS)%*%SX%*%Gam.PLS))%*%t(Gam.PLS)%*%SXY
  mu.PLS<-colMeans(Y)-crossprod(beta.PLS, colMeans(X))

  return(list(beta.PLS=beta.PLS, Gamma.PLS=Gam.PLS, mu.PLS=mu.PLS,
              Wu.PLS=resu_SIMPLS$Wu))
}



