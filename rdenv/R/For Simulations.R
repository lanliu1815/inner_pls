#' CLC Envelope Simulation
#'
#' Implement Section 4.1 in the paper of non-Grassmann algorithm (A Note on Fast
#' Envelope Estimation, Cook, Forzani and Su, 2016).
#'
#' @return \code{M, U, betaOLS, u, r, Gamma1=Gamma1, beta=beta, X, Y,
#'   epsilon}.

scenarios_nonGrassmann<-function(case){

  if(case==1){

    p=100; r=100; n=500; u=20

    eta<-matrix(runif(u*p, 0, 10), ncol=p, nrow=u)
    Gamma<-rsomat(r, r)
    Gamma1<-Gamma[,1:u]
    Gamma0<-Gamma[,(u+1):r]
    Sigma1<-sweep(Gamma1, MARGIN=2, runif(u, 49, 51), "*")%*%t(Gamma1)
    Sigma0<-sweep(Gamma0, MARGIN=2, runif(r-u, 49, 51), "*")%*%t(Gamma0)
    Sigma<-Sigma1+Sigma0
    beta<-Gamma1%*%eta

    X<-MASS::mvrnorm(n, rep(0, p), 400*diag(p))
    eps<-MASS::mvrnorm(n, rep(0, r), Sigma)
    Y<-tcrossprod(X, beta)+eps

  } else if(case==2){

    p=100; r=100; n=500; u=5

    eta<-matrix(runif(u*p, 0, 10), ncol=p, nrow=u)
    Gamma<-rsomat(r, r)
    Gamma1<-Gamma[,1:u]
    Gamma0<-Gamma[,(u+1):r]
    Sigma1<-sweep(Gamma1, MARGIN=2, rep(1,u), "*")%*%t(Gamma1)
    Sigma0<-sweep(Gamma0, MARGIN=2, rep(100, r-u), "*")%*%t(Gamma0)
    Sigma<-Sigma1+Sigma0
    beta<-Gamma1%*%eta

    X<-MASS::mvrnorm(n, rep(0, p), 25*diag(p))
    eps<-MASS::mvrnorm(n, rep(0, r), Sigma)
    Y<-tcrossprod(X, beta)+eps

  } else if(case==3||case==4){

    p=30; r=30; n=200; u=5

    eta<-matrix(runif(u*p, 0, 10), ncol=p, nrow=u)

    if(case==3){
      Gamma<-diag(r)
      Gamma1<-Gamma[,1:u]
      Gamma0<-Gamma[,(u+1):r]
      Sigma1<-sweep(Gamma1, MARGIN=2, 1.5^(1:u), "*")%*%t(Gamma1)
      Sigma0<-sweep(Gamma0, MARGIN=2, 1.5^((u+1):r), "*")%*%t(Gamma0)
    } else{
      Gamma<-rsomat(r, r)
      Gamma1<-Gamma[,1:u]
      Gamma0<-Gamma[,(u+1):r]
      Sigma1<-sweep(Gamma1, MARGIN=2, 1.05^(1:u), "*")%*%t(Gamma1)
      Sigma0<-sweep(Gamma0, MARGIN=2, 1.05^((u+1):r), "*")%*%t(Gamma0)
    }

    Sigma<-Sigma1+Sigma0
    beta<-Gamma1%*%eta

    X<-MASS::mvrnorm(n, rep(0, p), 100*diag(p))
    eps<-MASS::mvrnorm(n, rep(0, r), Sigma)
    Y<-tcrossprod(X, beta)+eps
  }


  sigY<-stats::cov(Y)*(n-1)/n
  sigYX<-stats::cov(Y, X)*(n-1)/n
  sigX<-stats::cov(X)*(n-1)/n
  invsigX<-chol2inv(chol(sigX))
  betaOLS<-sigYX%*%invsigX
  U<-tcrossprod(betaOLS, sigYX)
  M<-sigY-U

  return(list(M=M, U=U, betaOLS=betaOLS, u=u, r=r, Gamma1=Gamma1, beta=beta,
              Sigma=Sigma, X=X, Y=Y, epsilon=eps))
}


#' Data Generation for Envelope estimation.
#'
#' \code{p<u} is OK.
#'
#' \eqn{\eta} by independent uniform (0,10) variates.
#' @param SAM Whether to generate data of sample version.
#' @param u u<r.
#' @param dia The first \code{u} elements of \code{dia} are taken as the
#'   eigenvalues corresponding to the envelope subspace.
#'
#' @return \code{list(Sigma=Sigma, U=U, beta=beta, Gamma1=Gamma1,
#' Gamma0=Gamma0, eta=eta)} for population version.
#'
#' \code{return(list(beta=beta, betaOLS=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps,
#'  Gamma1=Gamma1, Gamma0=Gamma0, Sigma=Sigma))} for sample version.
data_gen.env<-function(SAM=F, u, p, r, n, x, dia, upper.eta=10){

  eta<-matrix(runif(u*p, 0, upper.eta), ncol=p, nrow=u)

  Gamma<-rsomat(r, r)
  Gamma1<-as.matrix(Gamma[,1:u])
  Gamma0<-as.matrix(Gamma[,(u+1):r])
#  Sigma1<-sweep(Gamma1, MARGIN=2, dia[1:u], "*")%*%t(Gamma1)
#  Sigma0<-sweep(Gamma0, MARGIN=2, dia[(u+1):r], "*")%*%t(Gamma0)
  Sigma<-sweep(Gamma, MARGIN=2, dia[1:r], "*")%*%t(Gamma)
  beta<-Gamma1%*%eta
  U<-tcrossprod(beta)

  if(!SAM){
    return(list(Sigma=Sigma, U=U, beta=beta, Gamma1=Gamma1, Gamma0=Gamma0, eta=eta))
  } else{

    X<-MASS::mvrnorm(n, rep(0, p), x*diag(p))
    eps<-MASS::mvrnorm(n, rep(0, r), Sigma)
    Y<-tcrossprod(X, beta)+eps

    sigY<-stats::cov(Y)*(n-1)/n
    sigYX<-stats::cov(Y, X)*(n-1)/n
    sigX<-stats::cov(X)*(n-1)/n
    invsigX<-chol2inv(chol(sigX))
    betaOLS<-sigYX%*%invsigX
    U<-tcrossprod(betaOLS, sigYX)
    M<-sigY-U

    return(list(beta=beta, betaOLS=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps, Gamma1=Gamma1,
                Gamma0=Gamma0, Sigma=Sigma))
  }

}

#' Data Generation for Inner Envelope Simulation
#'
#' The envelope subspace will always be the full space.
#'
#' \eqn{\eta_1, \eta_2} by \code{runif(.., 0, upper.eta)}.
#'
#' @param d Dimension of the inner envelope subspaces. \eqn{d=p} is allowed.
#' @param dia Eigenvalues of \eqn{M}, with the first \eqn{d} ones corresponding
#' to the \eqn{S_1} subspace.
#' @param x Used as \code{X<-MASS::mvrnorm(n, rep(0, p), x*diag(p))}, to generate
#' sample data matrix X.
#' @param VERBOSE Whether to print some extra warning information.
data_gen.ienv<-function(SAM=F, d, p, r, dia, n, x, upper.eta=10, VERBOSE=F){

  Gam<-rsomat(r, r)
  Gam1<-as.matrix(Gam[,1:d])
  Gam0<-as.matrix(Gam[,(d+1):r])

  if(d==p){
    if(VERBOSE){
      warning("Degenerate to envelope.\n")
    }
    eta1<-matrix(runif(d*p, 0, upper.eta), nrow=p)
    beta<-tcrossprod(Gam1,eta1)
    B<-Gam3<-0
  } else{
    B<-rsomat(p-d, r-d)
    B0<-qr.Q(qr(B), complete=T)[,(p-d+1):(r-d)]
    Gam3<-Gam0%*%B0

    eta1<-matrix(runif(d*p, 0, upper.eta), nrow=p)
    eta2<-matrix(runif(p*(p-d), 0, upper.eta), nrow=p)
    beta<-tcrossprod(Gam1,eta1)+tcrossprod((Gam0%*%B),eta2)
  }

  Sigma<-sweep(Gam, MARGIN=2, dia[1:r], "*")%*%t(Gam)
  U.pop<-tcrossprod(beta)*x

  if(!SAM){
    return(list(Sigma=Sigma, U.pop=U.pop, beta=beta, Gamma1=Gam1, Gamma12=Gam,
                Gamma3=Gam3, Gamma0=Gam0, B=B, eta1=eta1, eta2=eta2))
  } else{

    #X<-matrix(100*rbinom(n*p, 1, 0.5), nrow=n, ncol=p)
    X<-MASS::mvrnorm(n, rep(0, p), x*diag(p))
    eps<-MASS::mvrnorm(n, rep(0, r), Sigma)
    Y<-tcrossprod(X, beta)+eps

    sigY<-stats::cov(Y)*(n-1)/n
    sigYX<-stats::cov(Y, X)*(n-1)/n
    sigX<-stats::cov(X)*(n-1)/n
    invsigX<-chol2inv(chol(sigX))
    betaOLS<-sigYX%*%invsigX
    U<-tcrossprod(betaOLS, sigYX)
    M<-sigY-U

    return(list(beta=beta, betaOLS=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps,
                Gamma1=Gam1, Gamma12=Gam, Gamma3=Gam3, Gamma0=Gam0, B=B,
                Sigma=Sigma, U.pop=U.pop))
  }
}

#' Data Generation for Inner Envelope Simulation
#'
#' Improved version of \code{\link{data_gen.ienv}}. The \eqn{B} matrix can be
#' specified in the input.
#'
#' @param B The \eqn{B} matrix. When \code{d=p} this input is useless.
data_gen.ienv2<-function(SAM=F, d, p, r, dia, n, x, B=NULL, upper.eta=10, VERBOSE=F){

  Gam<-rsomat(r, r)
  Gam1<-as.matrix(Gam[,1:d])
  Gam0<-as.matrix(Gam[,(d+1):r])

  if(d==p){
    if(VERBOSE){
      warning("Degenerate to envelope.\n")
    }
    eta1<-matrix(runif(d*p, 0, upper.eta), nrow=p)
    beta<-tcrossprod(Gam1,eta1)
    B<-Gam3<-0
  } else {

    if(is.null(B)){
      B<-rsomat(p-d, r-d)
    }
    if(dim(B)[1]!=(r-d) | dim(B)[2]!=(p-d)){
      stop("Wrong dimension of B!")
    }
    B0<-qr.Q(qr(B), complete=T)[,(p-d+1):(r-d)]
    Gam3<-Gam0%*%B0

    eta1<-matrix(runif(d*p, 0, upper.eta), nrow=p)
    eta2<-matrix(runif(p*(p-d), 0, upper.eta), nrow=p)
    beta<-tcrossprod(Gam1,eta1)+tcrossprod((Gam0%*%B),eta2)
  }

  Sigma<-sweep(Gam, MARGIN=2, dia[1:r], "*")%*%t(Gam)
  U<-tcrossprod(beta)

  if(!SAM){
    return(list(Sigma=Sigma, U=U, beta=beta, Gamma1=Gam1, Gamma12=Gam, Gamma3=Gam3, Gamma0=Gam0, B=B))
  } else{

    X<-MASS::mvrnorm(n, rep(0, p), x*diag(p))
    eps<-MASS::mvrnorm(n, rep(0, r), Sigma)
    Y<-tcrossprod(X, beta)+eps

    sigY<-stats::cov(Y)*(n-1)/n
    sigYX<-stats::cov(Y, X)*(n-1)/n
    sigX<-stats::cov(X)*(n-1)/n
    invsigX<-chol2inv(chol(sigX))
    betaOLS<-sigYX%*%invsigX
    U<-tcrossprod(betaOLS, sigYX)
    M<-sigY-U

    return(list(beta=beta, betaOLS=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps,
                Gamma1=Gam1, Gamma12=Gam, Gamma3=Gam3, Gamma0=Gam0, B=B, Sigma=Sigma))
  }
}


#' Data Generation for Predictor Envelope Simulation
#'
#' Notice that \eqn{\beta} is the transposed version! That is, the one used in
#' Cook's notation for predictor envelopes.
#'
#' \eqn{\alpha} by \code{runif(.., 0, upper.eta)}.
#'
#' @param u Dimension of the envelope subspace.
#' @param eig.SiX Eigenvalues of \eqn{\Sigma_X}, with the first \eqn{u} ones corresponding
#' to the envelope subspace.
#' @param eig.SiYvX Eigenvalues of \eqn{\Sigma_{Y|X}}. This matrix is generated
#' by \code{A diag(eigSiYvX) A^T} with \eqn{A} being random generated orthogonal
#' matrix.
#' @param VERBOSE Whether to print some extra warning information.
#'
#' @return The reuturned \code{M, U} matrix is the one directly derived from the objective
#' function, on a general envelope sight. By simple envelope algebra it can be shown
#' that \eqn{E_M(U)} equals the predictor envelope subspace.
data_gen.pred_env<-function(SAM=F, u, r, p, eig.SiX, eig.SiYvX, n, upper.alpha=10, VERBOSE=F){

  Gam<-rsomat(p, p)
  Gam1<-as.matrix(Gam[,1:u])
  Gam0<-as.matrix(Gam[,(u+1):p])

  if(r>p){
    stop("Wrong as r>p!")
  }

  alpha<-matrix(runif(u*r, 0, upper.alpha), nrow=u)
    # beta is the transposed version!!!
  beta<-Gam1%*%alpha

  SiX<-sweep(Gam, MARGIN=2, eig.SiX[1:p], "*")%*%t(Gam)
  #  U<-tcrossprod(beta)  # It should hold that span(U)=span(beta^T)!

  A<-rsomat(r, r)
  SiYvX<-sweep(A, MARGIN=2, eig.SiYvX[1:r], "*")%*%t(A)

  if(!SAM){
    return(list(SiX=SiX, SiYvX=SiYvX, beta=beta, Gamma1=Gam1, Gamma=Gam, Gamma0=Gam0,
                alpha=alpha))
  } else{

    X<-MASS::mvrnorm(n, rep(0, p), SiX)
    eps<-MASS::mvrnorm(n, rep(0, r), SiYvX)
    Y<-X%*%beta+eps

    SY<-stats::cov(Y)*(n-1)/n
    SYX<-stats::cov(Y, X)*(n-1)/n
    SX<-stats::cov(X)*(n-1)/n
    SXY<-t(SYX)
    invSX<-chol2inv(chol(SX))
    invSY<-chol2inv(chol(SY))
    betaOLS<-invSX%*%t(SYX)
    U<-crossprod(SYX, invSY)%*%SYX
    M<-SX-U

    return(list(beta=beta, betaOLS=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps,
                Gamma1=Gam1, Gamma=Gam, Gamma0=Gam0, SiX=SiX, SiYvX=SiYvX, alpha=alpha,
                SX=SX, SYX=SYX, SXY=SXY, SY=SY))
  }
}


#' Data Generation for Inner Predictor Envelope Simulation
#'
#' The envelope subspace will always be the full space. Notice that \eqn{\beta}
#' is the transposed version! That is, the one used in Cook's notation for predictor
#' envelopes.
#'
#' If \eqn{d=r}, then \eqn{B,\eta_2} will be \code{NA}.
#'
#' \eqn{\eta_1, \eta_2} by \code{runif(.., 0, upper.eta)} when \code{special_case=F}.
#'
#' @param SAM Can only generate sample data with normal errors.
#' @param d Dimension of the inner envelope subspaces. \eqn{d=r} is allowed. But
#' \eqn{d>0}a is required.
#' @param eig.SiX Eigenvalues of \eqn{\Sigma_X}, with the first \eqn{d} ones corresponding
#' to the \eqn{S_1} subspace.
#' @param eig.SiYvX Eigenvalues of \eqn{\Sigma_{Y|X}}. This matrix is generated
#' by \code{A diag(eigSiYvX) A^T} with \eqn{A} being random generated orthogonal
#' matrix. Moreover, this matrix will only be used when generating sample data.
#' @param VERBOSE Whether to print some extra warning information.
#' @param special_case \eqn{\eta} is generated to satisfy the setting in the special case, but only
#' suitable for the case when SiYvX has same diagonal elements.. Specifically,
#' they are generated by 'cutting' an orthogonal matrix.
#' @param coef.eta2 Coefficient added to \eqn{\eta_2} to make \eqn{\beta} larger
#' in the \strong{special case}.
data_gen.IPE<-function(SAM=F, d, p, r, eig.SiX, eig.SiYvX, n, upper.eta=10, VERBOSE=F, special_case=F, coef.eta2=100){

  Gam<-rsomat(p, p)
  Gam1<-as.matrix(Gam[,1:d])
  Gam0<-as.matrix(Gam[,(d+1):p])

  if(r>p){
    stop("Wrong as r>p!")
  }

  SiX<-sweep(Gam, MARGIN=2, eig.SiX[1:p], "*")%*%t(Gam)
  #  U<-tcrossprod(beta)  # It should hold that span(U)=span(beta^T)!

  A<-rsomat(r, r)
  SiYvX<-sweep(A, MARGIN=2, eig.SiYvX[1:r], "*")%*%t(A)
  msqrSiYvX<-sweep(A, MARGIN=2, 1/sqrt(eig.SiYvX[1:r]), "*")%*%t(A)
  sqrSiYvX<-sweep(A, MARGIN=2, sqrt(eig.SiYvX[1:r]), "*")%*%t(A)

  if(d==r){
    if(special_case){
      stop("d==r.")
    }

    if(VERBOSE){
      warning("Degenerate to envelope.\n")
    }

    eta1<-matrix(runif(d*r, 0, upper.eta), nrow=r)
    beta<-tcrossprod(Gam1, eta1)
    B<-Gam3<-eta2<-NA

  } else{

    B<-rsomat(r-d, p-d)
    B0<-qr.Q(qr(B), complete=T)[,(r-d+1):(p-d)]
    Gam3<-Gam0%*%B0

    if(special_case){
      eta<-rsomat(r,r)
      eta1<-t(eta[,1:d,drop=F])
      eta2<-t(eta[,(d+1):r,drop=F])*coef.eta2
      eta1<-eta1%*%sqrSiYvX
      #    eta2<-eta2%*%sqrSiX
    } else {
      eta1<-matrix(runif(d*r, 0, upper.eta), nrow=d)
      eta2<-matrix(runif(r*(r-d), 0, upper.eta), nrow=r-d)
    }

   # beta is the transposed version!!!
    beta<-Gam1%*%eta1+Gam0%*%B%*%eta2
  }

  ret_pop<-list(SiX=SiX, SiYvX=SiYvX, beta=beta, eig.SiX=eig.SiX, eig.SiYvX=eig.SiYvX,
                Gamma1=Gam1, Gamma3=Gam3, Gamma0=Gam0, C=B, eta1=eta1, eta2=eta2)

  if(!SAM){
    return(ret_pop)
  } else{

    X<-MASS::mvrnorm(n, rep(0, p), SiX)
    eps<-MASS::mvrnorm(n, rep(0, r), SiYvX)
    Y<-X%*%beta+eps

    SY<-stats::cov(Y)*(n-1)/n
    SYX<-stats::cov(Y, X)*(n-1)/n
    SX<-stats::cov(X)*(n-1)/n
    invSX<-chol2inv(chol(SX))
    invSY<-chol2inv(chol(SY))
    betaOLS<-invSX%*%t(SYX)
    U<-crossprod(SYX, invSY)%*%SYX
    M<-SX-U

    ret_sam<-ret_pop
    ret_sam$betaOLS<-betaOLS
    ret_sam$M<-M
    ret_sam$U<-U
    ret_sam$Y<-Y
    ret_sam$X<-X
    ret_sam$epsilon<-eps
    ret_sam$SX<-SX
    ret_sam$SY<-SY
    ret_sam$SYX<-SYX
    return(ret_sam)
  }
}

#' Data generation for speicial case of IPE. OBSOLETE.
#'
#' Use \code{\link{data_gen.IPE}}.
data_gen.IPE2<-function(SAM=F, d, p, r, eig.SiX, eig.SiYvX, n, coef.eta2=100, VERBOSE=F){

  Gam<-rsomat(p, p)
  Gam1<-as.matrix(Gam[,1:d])
  Gam0<-as.matrix(Gam[,(d+1):p])

  if(r>p){
    stop("Wrong as r>p!")
  }

  SiX<-sweep(Gam, MARGIN=2, eig.SiX[1:p], "*")%*%t(Gam)
  #  U<-tcrossprod(beta)  # It should hold that span(U)=span(beta^T)!

  A<-rsomat(r, r)
  SiYvX<-sweep(A, MARGIN=2, eig.SiYvX[1:r], "*")%*%t(A)
  msqrSiYvX<-sweep(A, MARGIN=2, 1/sqrt(eig.SiYvX[1:r]), "*")%*%t(A)
  sqrSiYvX<-sweep(A, MARGIN=2, sqrt(eig.SiYvX[1:r]), "*")%*%t(A)

  if(d==r){
    stop("d==r.")
  } else{
    B<-rsomat(r-d, p-d)
    B0<-qr.Q(qr(B), complete=T)[,(r-d+1):(p-d)]
    Gam3<-Gam0%*%B0

    eta<-rsomat(r,r)
    eta1<-t(eta[,1:d,drop=F])
    eta2<-t(eta[,(d+1):r,drop=F])*coef.eta2
    eta1<-eta1%*%sqrSiYvX
#    eta2<-eta2%*%msqrSiX

    # beta is the transposed version!!!
    beta<-Gam1%*%eta1+Gam0%*%B%*%eta2
  }

  if(!SAM){
    return(list(SiX=SiX, SiYvX=SiYvX, beta=beta, Gamma1=Gam1, Gam=Gam, Gamma3=Gam3,
                Gamma0=Gam0, C=B, eta1=eta1, eta2=eta2, eig.SiX=eig.SiX))
  } else{

    X<-MASS::mvrnorm(n, rep(0, p), SiX)
    eps<-MASS::mvrnorm(n, rep(0, r), SiYvX)
    Y<-X%*%beta+eps

    SY<-stats::cov(Y)*(n-1)/n
    SYX<-stats::cov(Y, X)*(n-1)/n
    SX<-stats::cov(X)*(n-1)/n
    invSX<-chol2inv(chol(SX))
    invSY<-chol2inv(chol(SY))
    betaOLS<-invSX%*%t(SYX)
    U<-crossprod(SYX, invSY)%*%SYX  # Is this necessary? I wonder.
    M<-SX-U

    return(list(beta=beta, betaOLS=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps,
                Gamma1=Gam1, Gamma=Gam, Gamma3=Gam3, Gamma0=Gam0, C=B,
                SiX=SiX, SiYvX=SiYvX, SX=SX, SY=SY, SYX=SYX, eta1=eta1, eta2=eta2))
  }
}


#' Data Generation for Russian Doll Envelope Simulation.
#'
#' Gamma1 corresponds to the first u1 eigenvalues of Sigma, Gamma2 corrsponds to
#' the next u2 ones, and Gamma3 the rest.
#'
#' @param SAM Whether to use sample version.
#' @param u1,u2 u1<=p; u1+u2>=p; u2!=1. Because we assume that beta is of full rank.
#' @param upper.eta The upper bound for the distribution of eta, which is a
#'  uniform one, with the lower bound being 0.
data_gen.rdenv<-function(SAM=F, u1, u2, p, r, n, x, dia, upper.eta=10){

  Gam<-rsomat(r, r)

  Gam1<-as.matrix(Gam[,1:u1])
  Gam2<-as.matrix(Gam[,(u1+1):(u1+u2)])
  Gam12<-as.matrix(Gam[,1:(u1+u2)])
  Gam0<-as.matrix(Gam[,(u1+1):r])
  if(u1+u2==r){
    Gam3<-NULL
  } else{
    Gam3<-as.matrix(Gam[,(u1+u2+1):r])
  }
  B<-rsomat(p-u1, u2)

  Sigma<-sweep(Gam, MARGIN=2, dia[1:r], "*")%*%t(Gam)

  eta1<-matrix(runif(u1*p, min=0, max=upper.eta), nrow=u1)
  eta2<-matrix(runif((p-u1)*p, min=0, max=upper.eta), nrow=p-u1)
  beta<-Gam1%*%eta1+Gam2%*%B%*%eta2
#  s<-min(abs(c(p-u1, u2-1)))
#  eta1<-matrix(runif(u1*p, min=0, max=10), nrow=u1)  # p*u1, rank=u1
#  eta21<-rsomat(s, u2)  # s*u2
#  eta22<-matrix(runif(s*p, min=0, max=10), nrow=s)  # p*s
#  eta2<-eta21%*%eta22  # p*u2, with rank=s

  if(!SAM){
    return(list(Sigma=Sigma, beta=beta, Gamma1=Gam1, Gamma2=Gam2, Gamma12=Gam12,
                Gamma3=Gam3, Gamma0=Gam0, eta1=eta1, B=B, eta2=eta2,
                U=tcrossprod(beta), M=Sigma))
  } else{

    X<-MASS::mvrnorm(n, rep(0, p), x*diag(p))
    eps<-MASS::mvrnorm(n, rep(0, r), Sigma)
    Y<-tcrossprod(X, beta)+eps

    sigY<-stats::cov(Y)*(n-1)/n
    sigYX<-stats::cov(Y, X)*(n-1)/n
    sigX<-stats::cov(X)*(n-1)/n
    invsigX<-chol2inv(chol(sigX))
    betaOLS<-sigYX%*%invsigX
    U<-tcrossprod(betaOLS, sigYX)
    M<-sigY-U

    return(list(beta=beta, betaOLS=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps, Gamma1=Gam1,
                Gamma2=Gam2, Gamma12=Gam12, Gamma3=Gam3, Gamma0=Gam0, Sigma=Sigma, B=B))
  }
}


#' For Inner Envelope Studying.
#'
#' Generate two sets of inner envelope-structured data, the envelope structure of one is
#'  the full space and the other is not, and one set of envelope-structured data.
#'
#' @param SAM Whether to use sample version.
#' @param u2,u u2 columns of Gamma0 is contained within beta.part. u is the dimension
#' of the envelope-structured data.
#' @param upper.eta The upper bound for the distribution of eta, which is a
#'  uniform one, with the lower bound being 0.
data_gen.ioenv<-function(SAM=F, d, u2=p-d+1, u=r-d, p, r, n, x, dia, upper.eta=10){

  # Gammas, B and etas
  Gamma<-rsomat(r, r)
  Gamma1.ie<-as.matrix(Gamma[,1:d])
  Gamma0.ie<-as.matrix(Gamma[,(d+1):r])
  Gamma1.env<-as.matrix(Gamma[,(r-u+1):r])
  Gamma0.env<-as.matrix(Gamma[,1:(r-u)])

  B.ie.full<-rsomat(p-d, r-d)
  B.ie.part<-rsomat(p-d, u2)
  eta1.ie<-matrix(runif(d*p, 0, upper.eta), nrow=p)
  eta2.ie<-matrix(runif(p*(p-d), 0, upper.eta), nrow=p)
  eta.env<-matrix(runif(u*p, 0, upper.eta), ncol=p, nrow=u)

  # Sigma and betas
  Sigma<-sweep(Gamma, MARGIN=2, dia[1:r], "*")%*%t(Gamma)
  beta.env<-Gamma1.env%*%eta.env
  beta.ie.part<-tcrossprod(Gamma1.ie, eta1.ie)+tcrossprod((Gamma0.ie[,1:(p-d+1)]%*%B.ie.part), eta2.ie)
  beta.ie.full<-tcrossprod(Gamma1.ie, eta1.ie)+tcrossprod((Gamma0.ie%*%B.ie.full), eta2.ie)

  # Sample data
  SiX<-x*diag(p)
  X<-MASS::mvrnorm(n, rep(0, p), SiX)
  eps<-MASS::mvrnorm(n, rep(0, r), Sigma)
  Y.env<-tcrossprod(X, beta.env)+eps
  Y.ie.full<-tcrossprod(X, beta.ie.full)+eps
  Y.ie.part<-tcrossprod(X, beta.ie.part)+eps
  sigX<-stats::cov(X)*(n-1)/n
  invsigX<-chol2inv(chol(sigX))

  sigY.env<-stats::cov(Y.env)*(n-1)/n
  sigYX.env<-stats::cov(Y.env, X)*(n-1)/n
  betaOLS.env<-sigYX.env%*%invsigX
  U.env<-tcrossprod(betaOLS.env, sigYX.env)
  M.env<-sigY.env-U.env

  sigY.ie.full<-stats::cov(Y.ie.full)*(n-1)/n
  sigYX.ie.full<-stats::cov(Y.ie.full, X)*(n-1)/n
  betaOLS.ie.full<-sigYX.ie.full%*%invsigX
  U.ie.full<-tcrossprod(betaOLS.ie.full, sigYX.ie.full)
  M.ie.full<-sigY.ie.full-U.ie.full

  sigY.ie.part<-stats::cov(Y.ie.part)*(n-1)/n
  sigYX.ie.part<-stats::cov(Y.ie.part, X)*(n-1)/n
  betaOLS.ie.part<-sigYX.ie.part%*%invsigX
  U.ie.part<-tcrossprod(betaOLS.ie.part, sigYX.ie.part)
  M.ie.part<-sigY.ie.part-U.ie.part

  return(list(beta.env=beta.env, beta.ie.full=beta.ie.full, beta.ie.part=beta.ie.part,
              M.env=M.env, U.env=U.env, M.ie.full=M.ie.full, U.ie.full=U.ie.full,
              M.ie.part=M.ie.part, U.ie.part=U.ie.part, Gamma1.env=Gamma1.env,
              Gamma0.env=Gamma0.env, Gamma1.ie=Gamma1.ie, Gamma0.ie=Gamma0.ie,
              X=X, Y.env=Y.env, Y.ie.full=Y.ie.full, Y.ie.part=Y.ie.part,
              betaOLS.env=betaOLS.env, betaOLS.ie.full=betaOLS.ie.full,
              betaOLS.ie.part=betaOLS.ie.part, Sigma=Sigma))
}


#' To see if the data by \code{\link{data_gen.rdenv}} can be estimated well by the \code{\link[Renvlp]{envMU}} method.
#'
#' Pu by \code{\link[Renvlp]{ginv}}: \code{tol=.Machine$double.eps^0.8}.
#'
#' @param PRINT If true,the process of choosing starting values will be printed.
#' @return \code{list(d12.sv=d12.sv, a12.sv=a12.sv, d1.sv=d1.sv, a1.sv=a1.sv, d12.clc=d12.clc, a12.clc=a12.clc, d1.clc=d1.clc, a1.clc=a1.clc, noise.beta=noise.beta, noise.M=noise.M)}.
data.test<-function(PRINT=F, type="po", u1, u2, p, r, n, x, dia){
  #  u1<-10; u2<-10; p<-15; r<-30; n<-500; type="po"; dia<-4^(r:1); x<-100

  dat<-data_gen.rdenv(type=type, u1, u2, p, r, n, x, dia)
  M<-dat$M
  U<-dat$U
  beta<-dat$beta
  Gamma1.t<-dat$Gamma1
  Gamma2.t<-dat$Gamma2
  #  Gamma3<-dat$Gamma3
  #  Gamma23<-cbind(dat$Gamma2, dat$Gamma3)
  Gamma12.t<-dat$Gamma12
  real_beta<-dat$coef$beta
  real_M<-dat$coef$M
  #  real_U<-dat$coef$U
  noise.beta<-sqrt(svd(crossprod(beta-real_beta))$d[1])  # L2 norm of matrices
  noise.M<-sqrt(svd(crossprod(M-real_M))$d[1])
  #  noise.U<-sqrt(svd(crossprod(U-real_U))$d[1])

  eig.MU<-eigen(M+U)
  invMU<-sweep(eig.MU$vectors, MARGIN=2, 1/eig.MU$values, "*")%*%t(eig.MU$vectors)

  if(PRINT==T){
    sv<-starval2(M, U, beta, invMU, u1, u2, r, Gamma1.t, Gamma12.t)
  } else{
    sv<-starval(M, U, beta, invMU, u1, u2, r)
  }

  Gam1<-sv$Gamma1.12
  Gam2<-sv$Gamma2
  Gam12<-sv$Gamma12
  d1.sv<-frbnorm(Gam1, Gamma1.t)
  a1.sv<-angle(Gam1, Gamma1.t)
  d12.sv<-frbnorm(Gam12, Gamma12.t)
  a12.sv<-angle(Gam12, Gamma12.t)

  # use envMU()
  Pu<-U%*%Renvlp::ginv(crossprod(U), tol=.Machine$double.eps^0.8)%*%t(U)
  resu.inner<-Renvlp::envMU(M, diag(r)-Pu, r-u1)
  d1.clc<-frbnorm(resu.inner$Gamma0hat, Gamma1.t)
  a1.clc<-angle(resu.inner$Gamma0hat, Gamma1.t)

  resu<-Renvlp::envMU(M, U, u1+u2)
  d12.clc<-frbnorm(resu$Gammahat, Gamma12.t)
  a12.clc<-angle(resu$Gammahat, Gamma12.t)

  return(
    list(d12.sv=d12.sv, a12.sv=a12.sv, d1.sv=d1.sv, a1.sv=a1.sv, d12.clc=d12.clc, a12.clc=a12.clc, d1.clc=d1.clc, a1.clc=a1.clc, noise.beta=noise.beta, noise.M=noise.M)
  )
}

#' Alternating Algorithm Simulation
#'
#' Print the result of sv,clc and rd estimators and the estimation of beta by alternating algorithm.
#'
#' @param u1,u2 u1+u2!=r is the only constraint.
#' @param type This input must exist so it is always possible to check if the codes are right by running them on population version.
#' @param DTB When disturbance is used, envMU(initial=starting values.disturbed). While when DTB=F, envMU(initial=Gam1.sv) when estimating Gam1, where Gam1.sv comes from starval().
rd_alt_alg.test<-function(type="pop", u1=10, u2=10, p=15, r=30, n=500, x=400, dia=1.3^(30:1), DTB=F, PRINT=F){
  # type="po"; u1=15; u2=14; p=2; r=30; x=100; n=400;  dia=1.5^(c(2:30, 1)); DTB=F
  #---Data Preliminaries---

  dat<-data_gen.rdenv(type, u1, u2, p, r, n, x, dia)
  M<-dat$M
  U<-dat$U
  beta<-dat$beta
  real_beta<-dat$coef$beta
  Gamma1.t<-dat$Gamma1
  #  Gamma2<-dat$Gamma2
  #  Gamma3<-dat$Gamma3
  #  Gamma23<-cbind(dat$Gamma2, dat$Gamma3)
  Gamma12.t<-dat$Gamma12
  Gamma23.t<-qr.Q(qr(Gamma1.t), complete=T)[,(u1+1):r]

  eig.MU<-eigen(M+U)
  invMU<-sweep(eig.MU$vectors, MARGIN=2, 1/eig.MU$values, "*")%*%t(eig.MU$vectors)

  #---In the most outer function---

  if(DTB==F){
    sv<-starval(M, U, beta, invMU, u1, u2, r)
    Gam1<-Gam1.sv<-sv$Gamma1.12  # Choose this version of Gamma1
    Gam2<-sv$Gamma2
    Gam12<-Gam12.sv<-sv$Gamma12
    Gam23.sv<-qr.Q(qr(Gam1.sv), complete=T)[,(u1+1):r]
  } else{
    #---Disturbance---
    Gam1<-Gamma1.t
    Gam12<-Gamma12.t+MASS::mvrnorm(r, rep(0,u1+u2), 0.001*diag(u1+u2))  # Disturb Gam12
    Gam12<-qr.Q(qr(Gam12))
    Phi<-qr.solve(Gam12, Gam1)
    Phi<-qr.Q(qr(Phi))
    Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
    Gam1<-Gam12%*%Phi  # Now Gam1 is also disturbed
    Gam2<-Gam12%*%Phi0
    Gam12.sv<-Gam12
    Gam1.sv<-Gam1
    Gam23.sv<-qr.Q(qr(Gam1.sv), complete=T)[,(u1+1):r]
  }


  #---Alternating Algorithm---
  GEidx.g12<-Renvlp::GE(Gam12)
  ca<-Gam12%*%solve(Gam12[GEidx.g12[1:(u1+u2)],])  # This is the C_A without the rearrange of rows.
  GEidx.g1<-Renvlp::GE(Gam1)
  ce<-Gam1%*%solve(Gam1[GEidx.g1[1:u1],])
  cb1.0<-qr.solve(ca, ce)  # Render C_A, C_{B_1} in the outer function
  GEidx.cb1<-Renvlp::GE(cb1.0)
  cb1<-cb1.0%*%solve(cb1.0[GEidx.cb1[1:u1],])

  e1<-eigen(t(Gam1)%*%M%*%Gam1)
  e2<-eigen(t(Gam2)%*%M%*%Gam2)
  e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
  obj1<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))  # Value of obj, ok this way???
  #  cat("obj1: ", obj1, "\n")

  if(PRINT==F){
    rs.al<-alt_alg(u1, u2, r, M, invMU, ca, ce, cb1, GEidx.g1, GEidx.g12, GEidx.cb1, obj1)
  } else{
    rs.al<-alt_alg.p(u1, u2, r, M, invMU, ca, ce, cb1, GEidx.g1, GEidx.g12, GEidx.cb1, obj1)
  }
  #--- ---

  d1<-frbnorm(rs.al$Gamma1, Gamma1.t)
  a1<-angle(rs.al$Gamma1, Gamma1.t)
  d12<-frbnorm(rs.al$Gamma12, Gamma12.t)
  a12<-angle(rs.al$Gamma12, Gamma12.t)

  beta.rd<-tcrossprod(rs.al$Gamma12)%*%beta
  n.beta.rd<-sqrt(eigen(crossprod(beta.rd-real_beta))$values[1])

  d1.sv<-frbnorm(Gam1.sv, Gamma1.t)
  a1.sv<-angle(Gam1.sv, Gamma1.t)
  d12.sv<-frbnorm(Gam12.sv, Gamma12.t)
  a12.sv<-angle(Gam12.sv, Gamma12.t)

  #---Compare with envMU()---
  resu<-envMU2(M, U, u1+u2, initial=Gam12.sv)
  d12.clc<-frbnorm(resu$Gammahat, Gamma12.t)
  a12.clc<-angle(resu$Gammahat, Gamma12.t)

  Pu<-U%*%Renvlp::ginv(crossprod(U), tol=.Machine$double.eps^0.75)%*%t(U)
  resu.inner<-envMU2(M, diag(r)-Pu, r-u1, initial=Gam23.sv)
  d1.clc<-frbnorm(resu.inner$Gamma0hat, Gamma1.t)
  a1.clc<-angle(resu.inner$Gamma0hat, Gamma1.t)

  beta.clc<-tcrossprod(resu$Gammahat)%*%beta
  n.beta.clc<-sqrt(eigen(crossprod(beta.clc-real_beta))$values[1])

  beta.env_real<-tcrossprod(Gamma12.t)%*%beta
  n.beta.env_real<-sqrt(eigen(crossprod(beta.env_real-real_beta))$values[1])

  n.beta.ols<-sqrt(eigen(crossprod(beta-real_beta))$values[1])
  return(
    list(d1.sv=d1.sv, a1.sv=a1.sv, d1.clc=d1.clc, a1.clc=a1.clc, d1=d1, a1=a1,
         d12.sv=d12.sv, a12.sv=a12.sv, d12.clc=d12.clc, a12.clc=a12.clc, d12=d12, a12=a12,
         n.beta.clc=n.beta.clc, n.beta.rd=n.beta.rd, n.beta.env_real=n.beta.env_real,
         n.beta.ols=n.beta.ols)
  )
}



#---Some data---
dia.inner.origin<-c(1,5,10,50,100,500,1000,5000,10000,50000)
dia.inner<-c(50000,10000,5000,1000,500,100,50,10,5,1)
dia.inner2<-c(500000, 400000, 300000, 200000, 100000, 70000, 50000, 20000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10, 5, 1)
#---


#' Data Generation for Inner Envelope Simulation
#'
#' Implement Section 4.1 in the paper of inner envelope (Su, Cook, 2011).
#'
#' The inputs are set to conform to that of udfenv_coef(), so some of them are
#' actually useless. And they are preset to be "Hsg4".
#'
#' @param u1 Taking value of 1 or 4 or 7, according to the paper.
#' @return \code{list(M=Sig, U=U, beta=beta, Gamma1=Gam1, Gamma12=Gam, Gamma3=Gam3,
#'  Gamma0=Gam0)} for population version.
#'
#'  And \code{list(beta=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps, Gamma1=Gam1,
#'  Gamma12=Gam, Gamma3=Gam3, Gamma0=Gam0, coef=coef))} for sample version.
incoef_ie<-function(SAM=F, u1, u2="Hsg4", p, r, n=400, x="Hsg4", dia){

  d<-u1
  Gam<-rsomat(r, r)
  Gam1<-Gam[,1:d]
  Gam0<-Gam[,(d+1):r]
  B<-rsomat(p-d, r-d)
  B0<-qr.Q(qr(B), complete=T)[,(p-d+1):(r-d)]
  Gam3<-Gam0%*%B0

  Sig<-sweep(Gam, MARGIN=2, dia[1:r], "*")%*%t(Gam)

  eta1<-matrix(rnorm(d*p), nrow=p)
  eta2<-matrix(rnorm(p*(p-d)), nrow=p)

  beta<-tcrossprod(Gam1,eta1)+tcrossprod((Gam0%*%B),eta2)

  U<-tcrossprod(beta)
  coef<-list(M=Sig, U=U, beta=beta, Gamma1=Gam1, Gamma12=Gam, Gamma3=Gam3, Gamma0=Gam0, B=B)

  if(!SAM){
    return(coef)
  } else{

    X<-matrix(100*rbinom(n*p, 1, 0.5), nrow=n, ncol=p)
    eps<-MASS::mvrnorm(n, rep(0, r), Sig)
    Y<-tcrossprod(X, beta)+eps

    sigY<-stats::cov(Y)*(n-1)/n
    sigYX<-stats::cov(Y, X)*(n-1)/n
    sigX<-stats::cov(X)*(n-1)/n
    invsigX<-chol2inv(chol(sigX))
    betaOLS<-sigYX%*%invsigX
    U<-tcrossprod(betaOLS, sigYX)
    M<-sigY-U

    return(list(beta=betaOLS, M=M, U=U, Y=Y, X=X, epsilon=eps, Gamma1=Gam1, Gamma12=Gam, Gamma3=Gam3, Gamma0=Gam0, B=B, coef=coef))
  }

}

