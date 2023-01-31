#' Inner Predictor Envelope Estimation
#'
#' Overall estimation. Either mere using the starting values or further optimizing
#' by non-Grassmann technique are both good.
#'
#' @param d Dimension of the inner predictor envelope subspace. It should hold
#' that \eqn{0\le d \le r}.
#' @param initial Starting values for \eqn{\Gamma_1}. If set, will directly go
#' to the Non-Grassmann stage. When \code{d=0}, this will not be used, so no
#' specifying is fine.
#' @param NonGrass Whether to implement the optimization?
#' @param ftol.nonG Used in \code{ftol} in \code{\link{obj.IPE}}.
#'
#' @return
#' \item{Gamma1}{Will be \code{NULL} if \eqn{d=0}.}
#' \item{Gamma0}{Will be an orthogonal matrix if \eqn{d=0}.}
#' \item{betaIPE}{The IPE estimator of \eqn{\beta}.}
#' \item{MLE_otherpars}{Estimator for the other parameters. See details in
#' \code{\link{mle_after_Gamma1.IPE}}.}
#' \item{alphaIPE}{The intercept.}
#' \item{loglik}
IPE<-function(X, Y, d, initial=NULL, tol=.Machine$double.eps^0.75, NonGrass=F, ftol.nonG=0.1, PRINT=F){

  #---TEST---
#  d=0
  #--- ---

  # Preliminaries
  Y<-as.matrix(Y)
  X<-as.matrix(X)
  r<-ncol(Y)
  p<-dim(X)[2]
  n<-dim(X)[1]

  if (nrow(Y)!=n){
    stop("X and Y should have the same number of observations.")
  }

  SY<-stats::cov(Y)*(n-1)/n
  SYX<-stats::cov(Y, X)*(n-1)/n
  SXY<-t(SYX)
  SX<-stats::cov(X)*(n-1)/n
  invSX<-chol2inv(chol(SX))
  invSY<-chol2inv(chol(SY))
  betaOLS<-invSX%*%t(SYX)
  SYvX<-SY-t(SXY)%*%chol2inv(chol(SX))%*%SXY

  eig.SX<-eigen(SX, symmetric=T)
  invSX<-sweep(eig.SX$vectors, MARGIN=2, 1/eig.SX$values, "*")%*%t(eig.SX$vectors)
  sqrSX<-sweep(eig.SX$vectors, MARGIN=2, sqrt(eig.SX$values), "*")%*%t(eig.SX$vectors)
  msqrSX<-sweep(eig.SX$vectors, MARGIN=2, 1/sqrt(eig.SX$values), "*")%*%t(eig.SX$vectors)

  eig.SYvX<-eigen(SYvX, symmetric=T)
  invSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, 1/eig.SYvX$values, "*")%*%t(eig.SYvX$vectors)
  sqrSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, sqrt(eig.SYvX$values), "*")%*%t(eig.SYvX$vectors)
  msqrSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, 1/sqrt(eig.SYvX$values), "*")%*%t(eig.SYvX$vectors)

  if(d==0){

    Gam1<-NULL
    Gam0<-rsomat(p, p)
    MLE_otherpars<-list(betaIPE=betaOLS, betaOLS=betaOLS)
    alphahat<-colMeans(Y)-crossprod(MLE_otherpars$betaIPE, colMeans(X))

    ## In this case, we have the following D and E matrices.
#    D<-SYvX
#    E<-t(SXY)%*%invSX%*%SXY
    e4<-eigen((msqrSYvX%*%t(SXY)%*%invSX%*%SXY%*%msqrSYvX), symmetric=T, only.values=T)
    objfun<-sum(log(eig.SX$values))+sum(log(eig.SYvX$values))#+sum(log(1+abs(e4$values[(r-d+1):r])))

  } else{

    # Starting Values
    if(is.null(initial)){
      sv<-init.ienv(M=SX, U=SXY%*%t(SXY), d=d, p=r, PRINT=PRINT)
      Gam0.sv<-sv$init
      Gam1.sv<-qr.Q(qr(Gam0.sv), complete=T)[,(p-d+1):p,drop=F]
    } else{
      if(nrow(initial)!=p||ncol(initial)!=d){
          stop("The dimension of initial is wrong.")
      }
      Gam1.sv<-initial
      Gam0.sv<-qr.Q(qr(Gam1.sv), complete=T)[,(d+1):p,drop=F]
    }

    # Optimization of obj
    if(NonGrass){
      nonG<-obj.IPE(SX=SX, SXY=SXY, SYvX=SYvX, d=d, initial=Gam1.sv, ftol=ftol.nonG, PRINT=PRINT, tol=tol)
      Gam1<-nonG$Gamma1
      Gam0<-nonG$Gamma0

    } else{
      Gam1<-Gam1.sv
      Gam0<-Gam0.sv
    }

    # Estimation beta and so on
    MLE_otherpars<-mle_after_Gamma1.IPE(SX=SX, SXY=SXY, SYvX=SYvX, Gam1=Gam1.sv, r=r)

    # Intercept
    alphahat<-colMeans(Y)-crossprod(MLE_otherpars$betaIPE, colMeans(X))

    # Value of the objective function

    e1<-eigen(t(Gam1)%*%SX%*%Gam1, symmetric=T, only.values=T)
    e2<-eigen(t(Gam1)%*%invSX%*%Gam1, symmetric=T, only.values=T)

    QSxgo<-diag(p)-Pr(sqrSX%*%Gam1, tol=tol)
    P_QSxgoSxgz<-Pr(QSxgo%*%sqrSX%*%Gam0, tol=tol)
    Q_QSxgoSxgz<-diag(p)-P_QSxgoSxgz
    D<-SYvX+t(SXY)%*%msqrSX%*%QSxgo%*%Q_QSxgoSxgz%*%msqrSX%*%SXY
    E<-t(SXY)%*%msqrSX%*%P_QSxgoSxgz%*%msqrSX%*%SXY
    eig.D<-eigen(D, symmetric=T, only.values=F)
    msqrD<-sweep(eig.D$vectors, MARGIN=2, 1/sqrt(eig.D$values), "*")%*%t(eig.D$vectors)
    e4<-eigen((msqrD%*%E%*%msqrD), symmetric=T, only.values=T)

    objfun<-sum(log(eig.SX$values))+sum(log(e1$values))+sum(log(e2$values))+sum(log(eig.D$values))+sum(log(1+abs(e4$values[(r-d+1):r])))

  }

  loglik<--n/2*objfun

  return(list(Gamma1=Gam1, Gamma0=Gam0, betaIPE=MLE_otherpars$betaIPE,
              MLE_otherpars=MLE_otherpars, alphaIPE=alphahat, loglik=loglik))
}



#' Optimization of the Objective Function of IPE by Non-Grassmann Technique.
#'
#' Estimating inner predictor envelope subspace using non-Grassmann technique, which, on
#' on the same time, relies heavily on starting values.
#'
#' May need improvement.
#'
#' @param initial Starting values for \strong{\eqn{\Gamma_1}} and must be specified.
#' @param ftol Used in \code{abs(obj1-obj5)<ftol*abs(obj1)}.
#' @param PRINT Whether to print the value of the objective function in each iteration?
#' @param tol Used in \code{\link{Pr}}.
#' @return \code{Gamma1, Gamma0, objfun}.
obj.IPE<-function(SX, SXY, SYvX, d, initial, ftol=0.1, PRINT=F, tol=.Machine$double.eps^0.75){

  # FOR TEST
  #  SX=SX.out; SXY=SXY.out; SYvX=SYvX.out; d=d; initial=Gam1.sv; ftol=0.1; PRINT=T; tol=.Machine$double.eps^0.75

  # Preliminaries

  p<-dim(SX)[1]

  eig.SX<-eigen(SX, symmetric=T)
  invSX<-sweep(eig.SX$vectors, MARGIN=2, 1/eig.SX$values, "*")%*%t(eig.SX$vectors)
  sqrSX<-sweep(eig.SX$vectors, MARGIN=2, sqrt(eig.SX$values), "*")%*%t(eig.SX$vectors)
  msqrSX<-sweep(eig.SX$vectors, MARGIN=2, 1/sqrt(eig.SX$values), "*")%*%t(eig.SX$vectors)

  eig.SYvX<-eigen(SYvX, symmetric=T)
  invSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, 1/eig.SYvX$values, "*")%*%t(eig.SYvX$vectors)
  sqrSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, sqrt(eig.SYvX$values), "*")%*%t(eig.SYvX$vectors)
  msqrSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, 1/sqrt(eig.SYvX$values), "*")%*%t(eig.SYvX$vectors)

  Gam1<-initial
  Gam0<-qr.Q(qr(Gam1), complete=T)[,(d+1):p]

  # Get value of the objective function

  e1<-eigen(t(Gam1)%*%SX%*%Gam1, symmetric=T, only.values=T)
  e2<-eigen(t(Gam1)%*%invSX%*%Gam1, symmetric=T, only.values=T)

  QSxgo<-diag(p)-Pr(sqrSX%*%Gam1, tol=tol)
  P_QSxgoSxgz<-Pr(QSxgo%*%sqrSX%*%Gam0, tol=tol)
  Q_QSxgoSxgz<-diag(p)-P_QSxgoSxgz
  D<-SYvX+t(SXY)%*%msqrSX%*%QSxgo%*%Q_QSxgoSxgz%*%msqrSX%*%SXY
  E<-t(SXY)%*%msqrSX%*%P_QSxgoSxgz%*%msqrSX%*%SXY
  eig.D<-eigen(D, symmetric=T, only.values=F)
  #  eig.N<-eigen(N, symmetric=T, only.values=T)
  msqrD<-sweep(eig.D$vectors, MARGIN=2, 1/sqrt(eig.D$values), "*")%*%t(eig.D$vectors)
  e4<-eigen((msqrD%*%E%*%msqrD), symmetric=T, only.values=T)

  obj_added<-log(prod(eig.SX$values))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(eig.D$values))+sum(log(1+abs(e4$values[(r-d+1):r])))

  if(PRINT){
    cat("obj1: ", obj1, "\n")
  }


  # Iteration

  GEidx.go<-GE(Gam1)
  ca<-Gam1%*%solve(Gam1[GEidx.go[1:d],])  # This is the C_A without the rearrange of rows.
  caSXca<-crossprod(ca, (SX%*%ca))  # C_A^T SX C_A
  cainvSXca<-crossprod(ca, (invSX%*%ca))  # C_A^T SX^{-1} C_A
  #  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.go[(d+1):p],,drop=F], ca[GEidx.go[(d+1):p],,drop=F])+diag(d) # C_A^T C_A

  maxiter<-100
  i<-1
  while(i<maxiter){

    for(j in GEidx.go[(d+1):p]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.SX<-crossprod(ca[-j,], as.matrix(SX[-j,j]))/SX[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invSX<-crossprod(ca[-j,], as.matrix(invSX[-j,j]))/invSX[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      #      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.SX<-a+mcm.SX
      W.SX<-caSXca-tcrossprod(amcm.SX, amcm.SX)*SX[j,j]  # W_1
      amcm.invSX<-a+mcm.invSX
      W.invSX<-cainvSXca-tcrossprod(amcm.invSX, amcm.invSX)*invSX[j,j]
      #      amcm.U<-a+mcm.U
      #      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.SX<-chol2inv(chol(W.SX))
      invW.invSX<-chol2inv(chol(W.invSX))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.SX<-x+mcm.SX
        tmp.invSX<-x+mcm.invSX
        invWmcm.SX<-invW.SX%*%tmp.SX

        #        tmp.U<-x+mcm.U
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]

        cax<-ca  # Imagine, replace \ba in \ca with x. Then it is cax
        cax[j,]<-x

        invcaxSXcax<-invW.SX-tcrossprod(invWmcm.SX,invWmcm.SX)*SX[j,j]/as.numeric(1+crossprod(tmp.SX,invWmcm.SX)*SX[j,j])
        PSxG<-sqrSX%*%cax%*%invcaxSXcax%*%t(cax)%*%sqrSX
        QSxG<-diag(p)-PSxG

        invcaxcax<-invcaocao-tcrossprod(invcaocao%*%x)/as.numeric(1+x%*%(invcaocao%*%x))
        GG<-cax%*%invcaxcax%*%t(cax)

        P_QSxGSxGG<-Pr(QSxG%*%sqrSX%*%(diag(p)-GG), tol=tol)
        Q_QSxGSxGG<-diag(p)-P_QSxGSxGG

        D<-SYvX+t(SXY)%*%msqrSX%*%QSxG%*%(diag(p)-P_QSxGSxGG)%*%msqrSX%*%SXY
        E<-t(SXY)%*%msqrSX%*%P_QSxGSxGG%*%msqrSX%*%SXY
        eig.D<-eigen(D, symmetric=T, only.values=F)
        invD<-sweep(eig.D$vectors, MARGIN=2, 1/eig.D$values, "*")%*%t(eig.D$vectors)

        -2*log(1+x%*%(invcaocao%*%x))+log(1+SX[j,j]*crossprod(tmp.SX, invW.SX%*%tmp.SX))+log(1+invSX[j,j]*crossprod(tmp.invSX, (invW.invSX%*%tmp.invSX)))+sum(log(eig.D$values))+sum(log(1+abs(eigen((E%*%invD), symmetric=F,only.values=T)$values[(r-d+1):r])))

        #+sum(log(1+eigen((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*%(invW.M-tcrossprod(invWmcm.M,invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j])),symmetric=F,only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj, method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.SX<-a+mcm.SX
      caSXca<-W.SX+tcrossprod(amcm.SX, amcm.SX)*SX[j,j]  # Get back to C_A^T M C_A
      amcm.invSX<-a+mcm.invSX
      cainvSXca<-W.invSX+tcrossprod(amcm.invSX, amcm.invSX)*invSX[j,j]
      #      amcm.U<-a+mcm.U
      #      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    # Determine whether to continue the iteration
    Gam1<-qr.Q(qr(ca))
    Gam0<-qr.Q(qr(Gam1), complete=T)[,(d+1):p]
    e1<-eigen(t(Gam1)%*%SX%*%Gam1, symmetric=T, only.values=T)
    e2<-eigen(t(Gam1)%*%invSX%*%Gam1, symmetric=T, only.values=T)

    QSxgo<-diag(p)-Pr(sqrSX%*%Gam1, tol=tol)
    P_QSxgoSxgz<-Pr(QSxgo%*%sqrSX%*%Gam0, tol=tol)
    Q_QSxgoSxgz<-diag(p)-P_QSxgoSxgz
    D<-SYvX+t(SXY)%*%msqrSX%*%QSxgo%*%Q_QSxgoSxgz%*%msqrSX%*%SXY
    E<-t(SXY)%*%msqrSX%*%P_QSxgoSxgz%*%msqrSX%*%SXY
    eig.D<-eigen(D, symmetric=T, only.values=F)
    #    eig.N<-eigen(N, symmetric=T, only.values=T)
    msqrD<-sweep(eig.D$vectors, MARGIN=2, 1/sqrt(eig.D$values), "*")%*%t(eig.D$vectors)
    e4<-eigen((msqrD%*%E%*%msqrD), symmetric=T, only.values=T)

    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(eig.D$values))+sum(log(1+abs(e4$values[(r-d+1):r])))

    if(PRINT){
      cat("obj5: ", obj5, "\n")
    }

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  # FOR TEST
  #  frbnorm(Gam1.sv, dat$Gamma1); frbnorm(Gam1, dat$Gamma1)

  return(
    list(Gamma1=Gam1, Gamma0=Gam0, objfun=obj1)
  )
}


#' MLE of other parameters other than Gamma_1 of IPE.
#'
#' The \eqn{\beta}, OLS and IPE ones, are of the transposed version. That is,
#' the one used in Cook's notation for predictor envelopes.
#'
#' @param r Must be specified.
#'
#' @return
#' \item{betaIPE}
#' \item{betaOLS}
#' \item{SiYvX.IPE}
#' \item{U}{The \eqn{U} matrix. Semi-orthogonal.}
#' \item{C.IPE}{Will be \code{NA} if \eqn{d=r}.}
mle_after_Gamma1.IPE<-function(SX, SXY, SYvX, Gam1, p=dim(Gam1)[1], d=dim(Gam1)[2], r){

  eig.SX<-eigen(SX, symmetric=T)
  invSX<-sweep(eig.SX$vectors, MARGIN=2, 1/eig.SX$values, "*")%*%t(eig.SX$vectors)
  sqrSX<-sweep(eig.SX$vectors, MARGIN=2, sqrt(eig.SX$values), "*")%*%t(eig.SX$vectors)
  msqrSX<-sweep(eig.SX$vectors, MARGIN=2, 1/sqrt(eig.SX$values), "*")%*%t(eig.SX$vectors)

  eig.SYvX<-eigen(SYvX, symmetric=T)
  invSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, 1/eig.SYvX$values, "*")%*%t(eig.SYvX$vectors)
  sqrSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, sqrt(eig.SYvX$values), "*")%*%t(eig.SYvX$vectors)
  msqrSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, 1/sqrt(eig.SYvX$values), "*")%*%t(eig.SYvX$vectors)

  betaOLS<-invSX%*%SXY  # This is the transposed version of beta. That is, the one used in Cook's notation for predictor envelopes.

  if(d==r){
    betaIPE<-Gam1%*%chol2inv(chol(t(Gam1)%*%SX%*%Gam1))%*%t(Gam1)%*%SXY
    C.st<-NA
    SiYvX.st<-NA  # NEEDS REVISION because it exists.
#    cat("1\n")
    Gam0<-qr.Q(qr(Gam1), complete=T)[,(d+1):p,drop=F]

  } else{

    Gam0<-qr.Q(qr(Gam1), complete=T)[,(d+1):p,drop=F]

    # MLE of other parameters

    #  eig.goSXgo<-eigen(crossprod(Gam1, SX%*%Gam1), symmetric=T)
    #  invgoSXgo<-sweep(eig.goSXgo$vectors, MARGIN=2, 1/eig.goSXgo$values, "*")%*%t(eig.goSXgo$vectors)
    PgoSx<-Gam1%*%chol2inv(chol(crossprod(Gam1, SX%*%Gam1)))%*%t(Gam1)%*%SX  # P_{Gam1(SX)}, non-standard inner product!
    QgoSx<-diag(p)-PgoSx
    QSxgo<-diag(p)-Pr(sqrSX%*%Gam1)  # Q_{sqrSX%*%Gam1}
    PQSxgz<-Pr(QSxgo%*%sqrSX%*%Gam0)
    QQSxgz<-diag(p)-PQSxgz

    # SiYvX
    M<-SYvX+t(SXY)%*%msqrSX%*%QSxgo%*%QQSxgz%*%msqrSX%*%SXY
    N<-t(SXY)%*%msqrSX%*%PQSxgz%*%msqrSX%*%SXY
    eig.M<-eigen(M, symmetric=T)
    invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)
    sqrM<-sweep(eig.M$vectors, MARGIN=2, sqrt(eig.M$values), "*")%*%t(eig.M$vectors)
    msqrM<-sweep(eig.M$vectors, MARGIN=2, 1/sqrt(eig.M$values), "*")%*%t(eig.M$vectors)

    term4VK<-msqrM%*%N%*%msqrM
    eig.term4VK<-eigen(term4VK, symmetric=T)
    V<-eig.term4VK$vectors
    K<-diag(c(rep(0, r-d), eig.term4VK$values[(r-d+1):r]))
    # !".st" means star, namely, the MLE of that parameter.
    SiYvX.st<-M+sqrM%*%V%*%K%*%t(V)%*%sqrM

    term4U<-PQSxgz%*%msqrSX%*%SXY%*%chol2inv(chol(SiYvX.st))%*%t(SXY)%*%msqrSX%*%PQSxgz  # ! One more QSxgo added
    eig.term4U<-eigen(term4U, symmetric=T)
    U<-eig.term4U$vectors[,1:(r-d)]

    betaIPE<-(PgoSx+msqrSX%*%tcrossprod(U)%*%sqrSX)%*%betaOLS  # Notice beta should be the transposed version.

    C.st<-qr.Q(qr(crossprod(Gam0, msqrSX%*%U)))

  }

  return(list(betaIPE=betaIPE, betaOLS=betaOLS, SiYvX.IPE=SiYvX.st, C.IPE=C.st,
              Gam1.IPE=Gam1, Gam0.IPE=Gam0))
}


#' The Inner PLS Algorithm
#'
#' Invokes \code{\link{SIMPLS_gamma}} and \code{\link{IPE}}.
#'
#' @param d Dimension of the inner predictor envelope subspace. It should hold
#' that \eqn{0 \leq d \leq r}.
#'
#' @return
#' \item{Gamma1}
#' \item{Gamma0}
#' \item{betaiPLS}
#' \item{alphaiPLS}
#' \item{Wu}{The original \eqn{W_u} matrix from \code{SIMPLS_gamma}.}
#' \item{paras}{The results from \code{\link{IPE}}. Includes some other parameters not so often used.}
innerPLS<-function(X, Y, d){
  #---TEST---
#  d=0; r=6; p=10; n<-200; upper.eta=10; PRINT=F; ftol=1e-1
#  eigSiYvX=rep(50, r); eigSiX=c(rep(5e1,d), seq(from=5e-2, to=1e0, length.out=p-d))
#  dat<-data_gen.IPE(SAM=T, d=d, p=p, r=r, n=n, eig.SiX=eigSiX, eig.SiYvX=eigSiYvX, upper.eta=upper.eta, special_case=T)
#  X<-dat$X; Y<-dat$Y
  #--- ---

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  dimY<-dim(Y)
  n<-dimY[1]
  r<-dimY[2]
  p<-ncol(X)

  SX<-stats::cov(X)*(n-1)/n
  SXY<-stats::cov(X, Y)*(n-1)/n

  if(d!=0){
    resu.SIMPLS<-SIMPLS_gamma(M=SX, U=diag(p)-Pr(SXY), u=p-d)
    Gam0.iPLS<-resu.SIMPLS$Gammahat
    Gam1.iPLS<-qr.Q(qr(Gam0.iPLS), complete=T)[,(p-d+1):p,drop=F]
    paras<-IPE(X=X, Y=Y, d=d, initial=Gam1.iPLS, NonGrass=F)
    return(list(Gamma1=Gam1.iPLS, Gamma0=Gam0.iPLS, betaiPLS=paras$betaIPE,
                alphaiPLS=paras$alphaIPE, Wu=resu.SIMPLS$Wu, paras=paras))

  } else{
    Gam1.iPLS<-NULL
    Gam0.iPLS<-rsomat(p, p)
    paras<-IPE(X=X, Y=Y, d=d, initial=Gam1.iPLS, NonGrass=F)
    return(list(Gamma1=Gam1.iPLS, Gamma0=Gam0.iPLS, betaiPLS=paras$betaIPE,
                alphaiPLS=paras$alphaIPE, paras=paras))

  }

}


#' Dimension Selection of IPE by AIC&BIC.
#'
#' @param NonGrass Whether to implement further optimization in each dimension?
#'
#' @return A list.
#' \item{d.aic}
#' \item{d.bic}
#' \item{loglik.seq}
#' \item{aic.seq}
#' \item{bic.seq}
aicbic.IPE<-function(X, Y, NonGrass=F, trace=0){

  #  X <- as.matrix(X) #?
  dimY<-dim(Y)
  n<-dimY[1]
  r<-dimY[2]
  p<-ncol(X)
  loglik.seq<-unlist(lapply(0:r, function(x) IPE(X=X, Y=Y, d=x, NonGrass=NonGrass)$loglik))
  npara.seq<-r^2+p*(p+1)/2+(p-r)*(r:0) # r^2+(r-d)(p-r)+p(p+1)/2
  aic.seq<--2*loglik.seq+2*npara.seq
  bic.seq <--2*loglik.seq+log(n)*npara.seq
  d.aic<-which.min(aic.seq)-1
  d.bic<-which.min(bic.seq)-1


  if(trace==-1){
    cat("-")
  }
  return(list(d.aic=d.aic, d.bic=d.bic, loglik.seq=loglik.seq,
              aic.seq=aic.seq, bic.seq=bic.seq))
}


#' Dimension Selection of inner PLS, IPE, PLS and Predictor Envelope by Bootstrap.
#'
#' Dimension Selection of inner PLS, IPE, PLS and Predictor Envelope by bootstrapping
#' the data to acquire difference bases. Finally choose a dimension with relatively
#' the smallest variability.
#'
#' Implement bootstrap to get the average distance, which is measured by vector
#' correlation, between each bootstrapped basis and the basis obtained from the
#' original data. Then select the largest(predictor envelope and PLS) or
#' smallest(inner PLS and IPE) dimension with corresponding distance over
#' \code{thres}.
#'
#' @param thres The threshold to use.
#' @param estmrs A character string containing the estimators of which to choose
#' dimension.
#' @param trace If 1, information is printed. If -1, a "-" is printed after the
#' whole function ends. The -1 part is useful when replicating this function.
#'
#' @return A list. If the method is not performed, then no return values. If no
#' estimators is specified, then \code{NULL} is returned.
#' \item{d.boot.IPE}{The selected dimension. If not any of the vector correlations
#' exceeded the \code{thres}, \code{NA} will be returned}
#' \item{dis.IPE}{The bootstrap mean of vector correlation value of each dimension.}
#' \item{d.boot.iPLS}
#' \item{dis.iPLS}
#' \item{u.boot.PLS}
#' \item{dis.PLS}
#' \item{u.boot.PredEnv}
#' \item{dis.PredEnv}
d.boot.sim<-function(X, Y, B=50, thres=0.9, estmrs=c("iPLS", "IPE", "PredEnv", "SIMPLS"), trace=0){

  #---TEST---
#    d=4; r=6; p=10; n<-200; upper.eta=10; PRINT=F; ftol=1e-1
#    eigSiYvX=rep(50, r); eigSiX=c(rep(5e1,d), seq(from=5e-2, to=1e0, length.out=p-d))
#    dat<-data_gen.IPE(SAM=T, d=d, p=p, r=r, n=n, eig.SiX=eigSiX, eig.SiYvX=eigSiYvX, upper.eta=upper.eta, special_case=T)
#    X<-dat$X; Y<-dat$Y
#    B=50; thres=0.9; trace=0
#    estmrs=c("iPLS", "IPE", "PredEnv", "SIMPLS")
  #--- ---

  esti_nm<-c("iPLS", "IPE", "PredEnv", "SIMPLS")

  # If no estimators is specified, then end directly.
  if(!length(estmrs)){
#    cat("No estimators by bootstrap!\n")
    return(NULL)  # Some better function than return() here?
  }

  if(!all(estmrs%in%esti_nm)){
    stop("Specify the right estimator in estmrs!")
  }


  dimY<-dim(Y)
  n<-dimY[1]
  r<-dimY[2]
  p<-ncol(X)

  SX<-stats::cov(X)*(n-1)/n
  SXY<-stats::cov(X, Y)*(n-1)/n
  SYX<-t(SXY)
  SY<-stats::cov(Y)*(n-1)/n
  invSY<-chol2inv(chol(SY))
  U<-crossprod(SYX, invSY)%*%SYX
  M<-SX-U

  X.list<-list(NULL)
  Y.list<-list(NULL)

  for(i in 1:B){
    idx<-sample(1:n, n, replace = T)
    X.list[[i]]<-X[idx,]
    Y.list[[i]]<-Y[idx,]
  }

  dim_opt<-list()

  # IPE

  if("IPE"%in%estmrs){
    if(trace==1){
      cat("\nIPE\n")
    }

    dis.SV<-numeric(r-1)

#    bootSV<-function(j){
#      return(vec_cor(Gam1.ori.SV, IPE(X.list[[j]], Y.list[[j]], d.boot.SV, NonGrass=F)$Gamma1))
#    }

    for(d.boot.SV in 1:r){
      Gam1.ori.SV<-IPE(X, Y, d.boot.SV, NonGrass=F)$Gamma1
      dis.boot.SV<-unlist(lapply(1:B, function(j) vec_cor(Gam1.ori.SV, IPE(X.list[[j]], Y.list[[j]], d.boot.SV, NonGrass=F)$Gamma1)))
      dis.SV[d.boot.SV]<-mean(dis.boot.SV)

      if(trace==1){
        cat(d.boot.SV, ":", dis.boot.SV, "\n\n")
      }
    }

    d.boot.SV<-sort(which(dis.SV>thres), decreasing=F)[1]
    if(is.na(d.boot.SV)){
      d.boot.SV<-0
    }

    dim_opt$d.boot.IPE<-d.boot.SV
    dim_opt$dis.IPE<-dis.SV

  }

  # Inner PLS
  if("iPLS"%in%estmrs){

    if(trace==1){
      cat("\nInnerPLS\n")
    }

    dis.inp<-numeric(r-1)

#    bootinp<-function(j){
#      return(vec_cor(Gam0.ori.inp, SIMPLS_gamma(M=stats::cov(X.list[[j]])*(n-1)/n, U=diag(p)-Pr(stats::cov(X.list[[j]], Y.list[[j]])*(n-1)/n), u=p-d.boot.inp)$Gammahat))
#    }

    for(d.boot.inp in 1:r){
      Gam0.ori.inp<-SIMPLS_gamma(M=SX, U=diag(p)-Pr(SXY), u=p-d.boot.inp)$Gammahat
      dis.boot.inp<-unlist(lapply(1:B, function(j) vec_cor(Gam0.ori.inp, SIMPLS_gamma(M=stats::cov(X.list[[j]])*(n-1)/n, U=diag(p)-Pr(stats::cov(X.list[[j]], Y.list[[j]])*(n-1)/n), u=p-d.boot.inp)$Gammahat)))
      dis.inp[d.boot.inp]<-mean(dis.boot.inp)

      if(trace==1){
        cat(d.boot.inp, ":", dis.boot.inp, "\n\n")
      }
    }

    d.boot.inp<-sort(which(dis.inp>thres), decreasing=F)[1]
    if(is.na(d.boot.inp)){
      d.boot.inp<-0
    }

    dim_opt$d.boot.iPLS<-d.boot.inp
    dim_opt$dis.iPLS<-dis.inp

  }

  # PLS

  if("SIMPLS"%in%estmrs){
    if(trace==1){
      cat("\nPLS\n")
    }

    dis.PLS<-numeric(p-1)

#    bootPLS<-function(j){
#      return(vec_cor(Gam.ori.PLS, SIMPLS_gamma(M=stats::cov(X.list[[j]])*(n-1)/n, U=tcrossprod(stats::cov(X.list[[j]], Y.list[[j]])*(n-1)/n), u=u.boot.PLS)$Gammahat))
#    }
    for(u.boot.PLS in 1:(p-1)){
      Gam.ori.PLS<-SIMPLS_gamma(M=SX, U=tcrossprod(SXY), u=u.boot.PLS)$Gammahat
      dis.boot.PLS<-unlist(lapply(1:B, function(j) vec_cor(Gam.ori.PLS, SIMPLS_gamma(M=stats::cov(X.list[[j]])*(n-1)/n, U=tcrossprod(stats::cov(X.list[[j]], Y.list[[j]])*(n-1)/n), u=u.boot.PLS)$Gammahat)))
      dis.PLS[u.boot.PLS]<-mean(dis.boot.PLS)

      if(trace==1){
        cat(u.boot.PLS, ":", dis.boot.PLS, "\n\n")
      }
    }

    u.boot.PLS<-sort(which(dis.PLS>thres), decreasing=T)[1]

    dim_opt$u.boot.PLS<-u.boot.PLS
    dim_opt$dis.PLS<-dis.PLS

  }

  # PredEnv

  if("PredEnv"%in%estmrs){
    if(trace==1){
      cat("\nPredictor Envelope\n")
    }

    dis.predenv<-numeric(p-1)

    bootpredenv<-function(j){
      SX.temp<-stats::cov(X.list[[j]])*(n-1)/n
      SYX.temp<-stats::cov(Y.list[[j]], X.list[[j]])*(n-1)/n
      SY.temp<-stats::cov(Y.list[[j]])*(n-1)/n
      invSY.temp<-chol2inv(chol(SY.temp))
#      U.temp<-crossprod(SYX.temp,invSY.temp)%*%SYX.temp
#      M.temp<-SX.temp-U.temp
      #    Gam.boot<-Renvlp::xenv(X=X.list[[j]], Y=Y.list[[j]], u=u.boot.predenv, asy=F)$Gamma
#      Gam.boot<-stv(M=SX.temp-crossprod(SYX.temp,invSY.temp)%*%SYX.temp, U=crossprod(SYX.temp,invSY.temp)%*%SYX.temp, u=u.boot.predenv)$init
      return(vec_cor(Gam.ori.predenv, stv(M=SX.temp-crossprod(SYX.temp,invSY.temp)%*%SYX.temp, U=crossprod(SYX.temp,invSY.temp)%*%SYX.temp, u=u.boot.predenv)$init))
    }

    for(u.boot.predenv in 1:(p-1)){
      Gam.ori.predenv<-stv(M=M, U=U, u=u.boot.predenv)$init
      dis.boot.predenv<-unlist(lapply(X=1:B, FUN=bootpredenv))
      dis.predenv[u.boot.predenv]<-mean(dis.boot.predenv)

      if(trace==1){
        cat(u.boot.predenv, ":", dis.boot.predenv, "\n\n")
      }
    }

    u.boot.predenv<-sort(which(dis.predenv>thres), decreasing=T)[1]

    dim_opt$u.boot.PredEnv<-u.boot.predenv
    dim_opt$dis.PredEnv<-dis.predenv

  }

  if(trace==-1){
    cat("-")
  }

  return(dim_opt)
}


#' Dimension selection of SIMPLS by cross-validation.
#'
#' @param trace If 1, a "-" is printed after each fold is completed.
#' If -1, a "-" is printed after the whole function ends.
#'
#' @return
#' \item{u.cv.SIMPLS}
#' \item{pred_err}
u.cv.SIMPLS<-function(X, Y, kfold=3, trace=0){

  #test---
  #  kfold<-3
  #--- ---

  X<-as.matrix(X)
  a<-dim(Y)
  n<-a[1]
  r<-a[2]
  p<-ncol(X)

  prederr<-rep(0, p)
  for(u in 1:p){
    id<-sample(n, n)
    Xn<-as.matrix(X[id, ])
    Yn<-as.matrix(Y[id, ])

    for(j in 1:kfold){
      id.test<-(floor((j-1)*n/kfold)+1):ceiling(j*n/kfold)
      id.train<-setdiff(1:n, id.test)
      X.train<-Xn[id.train, ]
      Y.train<-Yn[id.train, ]
      X.test<-Xn[id.test, ]
      Y.test<-Yn[id.test, ]
      n.test<-length(id.test)

      #      fit <- env(X.train, Y.train, u, asy = F)
      #      betahat <- fit$beta
      #      muhat <- fit$mu
      #      resi <- as.matrix(Y.test - matrix(1, n.test, 1) %*% t(muhat) - as.matrix(X.test) %*% t(betahat))

      fit.PLS<-SIMPLS_env(X.train, Y.train, u)
      resi<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.PLS$mu.PLS)-as.matrix(X.test)%*%fit.PLS$beta.PLS)
      sprederr<-apply(resi, 1, function(x) sum(x^2))
      prederr[u]<-prederr[u]+sum(sprederr)

    }
    if(trace==1){
      cat("-")
    }

  }

  if(trace==1){
    cat("\n")
  }

  if(trace==-1){
    cat("-")
  }

  prederr<-sqrt(prederr/n)
  u.cv<-which.min(prederr)

  return(list(u.cv.SIMPLS=u.cv, pred_err=prederr))
}


#' Dimension Selection of Predictor Envelope by Cross Validation
#'
#' Based on \code{\link{u.cv.sim}}.
#'
#' @return
#' \item{u.cv.PredEnv}
#' \item{prederr.PredEnv}
u.cv.PredEnv<-function(X, Y, kfold=3, trace=0){

  #---TEST---
#  d=5; r=6; p=10; n<-200; upper.eta=10; PRINT=F; ftol=1e-1
#  eigSiYvX=rep(50, r); eigSiX=c(rep(5e1,d), seq(from=5e-2, to=1e0, length.out=p-d))
#  dat<-data_gen.IPE(SAM=T, d=d, p=p, r=r, n=n, eig.SiX=eigSiX, eig.SiYvX=eigSiYvX, upper.eta=upper.eta, special_case=T)
#  X<-dat$X; Y<-dat$Y
#  kfold=3; NonGrass=F; trace=-1
  #--- ---


  X<-as.matrix(X)
  a<-dim(Y)
  n<-a[1]
  r<-a[2]
  p<-ncol(X)

  prederr.PredEnv<-rep(0, p)
  for(dimsn in 1:p){
    id<-sample(n, n)
    Xn<-as.matrix(X[id, ])
    Yn<-as.matrix(Y[id, ])

    for(j in 1:kfold){
      id.test<-(floor((j-1)*n/kfold)+1):ceiling(j*n/kfold)
      id.train<-setdiff(1:n, id.test)
      X.train<-Xn[id.train, ]
      Y.train<-Yn[id.train, ]
      X.test<-Xn[id.test, ]
      Y.test<-Yn[id.test, ]
      n.test<-length(id.test)

      # Predictor Envelope
      fit.predenv<-Renvlp::xenv(X.train, Y.train, u=dimsn, asy=F)
      resi.predenv<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.predenv$mu)-as.matrix(X.test)%*%fit.predenv$beta)
      sprederr.predenv<-apply(resi.predenv, 1, function(x) sum(x^2))
      prederr.PredEnv[dimsn]<-prederr.PredEnv[dimsn]+sum(sprederr.predenv)

    }
    if(trace==1){
      cat("-")
    }

  }

  if(trace==1){
    cat("\n")
  }

  if(trace==-1){
    cat("-")
  }

  prederr.PredEnv<-sqrt(prederr.PredEnv/n)
  u.cv.PredEnv<-which.min(prederr.PredEnv)

  return(
    list(u.cv.PredEnv=u.cv.PredEnv, prederr.PredEnv=prederr.PredEnv)
  )
}


#' Dimension selection of inner PLS, IPE by cross-validation.
#'
#' The range is from \eqn{0} to \eqn{r}.
#'
#' @param kfold The folds to use in CV.
#' @param trace If 1, a "-" is printed after each fold is completed.
#' If -1, a "-" is printed after the whole function ends.
#'
#' @return
#' \item{d.cv.iPLS}{The selected dimension.}
#' \item{prederr.iPLS}{Prediction errors of each dimension.}
#' \item{d.cv.IPE}
#' \item{prederr.IPE}
u.cv.sim<-function(X, Y, kfold=3, trace=0, NonGrass=F){

  #---TEST---
#  d=5; r=6; p=10; n<-200; upper.eta=10; PRINT=F; ftol=1e-1
#  eigSiYvX=rep(50, r); eigSiX=c(rep(5e1,d), seq(from=5e-2, to=1e0, length.out=p-d))
#  dat<-data_gen.IPE(SAM=T, d=d, p=p, r=r, n=n, eig.SiX=eigSiX, eig.SiYvX=eigSiYvX, upper.eta=upper.eta, special_case=T)
#  X<-dat$X; Y<-dat$Y
#  kfold=3; NonGrass=F; trace=-1
  #--- ---


  X<-as.matrix(X)
  a<-dim(Y)
  n<-a[1]
  r<-a[2]
  p<-ncol(X)

  prederr.inp<-prederr.IPE<-rep(0, r+1)
  for(dimsn in 0:r){
    id<-sample(n, n)
    Xn<-as.matrix(X[id, ])
    Yn<-as.matrix(Y[id, ])

    for(j in 1:kfold){
      id.test<-(floor((j-1)*n/kfold)+1):ceiling(j*n/kfold)
      id.train<-setdiff(1:n, id.test)
      X.train<-Xn[id.train, ]
      Y.train<-Yn[id.train, ]
      X.test<-Xn[id.test, ]
      Y.test<-Yn[id.test, ]
      n.test<-length(id.test)

      # SIMPLS
#      fit.PLS<-SIMPLS_env(X.train, Y.train, u)
#      resi.SIMPLS<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.PLS$mu.PLS)-as.matrix(X.test)%*%fit.PLS$beta.PLS)
#      sprederr.SIMPLS<-apply(resi.SIMPLS, 1, function(x) sum(x^2))
#      prederr.SIMPLS[u]<-prederr.SIMPLS[u]+sum(sprederr.SIMPLS)

      # inner PLS
      fit.inp<-innerPLS(X=X.train, Y=Y.train, d=dimsn)
      resi.inp<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.inp$alphaiPLS)-as.matrix(X.test)%*%fit.inp$betaiPLS)
      sprederr.inp<-apply(resi.inp, 1, function(x) sum(x^2))
      prederr.inp[dimsn+1]<-prederr.inp[dimsn+1]+sum(sprederr.inp)

      # IPE

      fit.IPE<-IPE(X=X.train, Y=Y.train, d=dimsn, NonGrass=NonGrass)
      resi.IPE<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.IPE$alphaIPE)-as.matrix(X.test)%*%fit.IPE$betaIPE)
      sprederr.IPE<-apply(resi.IPE, 1, function(x) sum(x^2))
      prederr.IPE[dimsn+1]<-prederr.IPE[dimsn+1]+sum(sprederr.IPE)

    }

    if(trace==1){
      cat("-")
    }

  }

  if(trace==1){
    cat("\n")
  }

  if(trace==-1){
    cat("-")
  }

  prederr.inp<-sqrt(prederr.inp/n)
  d.cv.inp<-which.min(prederr.inp)-1

  prederr.IPE<-sqrt(prederr.IPE/n)
  d.cv.IPE<-which.min(prederr.IPE)-1

  return(
    list(d.cv.iPLS=d.cv.inp, prederr.iPLS=prederr.inp, d.cv.IPE=d.cv.IPE, prederr.IPE=prederr.IPE)
    )
}


#' Cross-validation prediction error for inner PLS, IPE, Predictor Envelope and SIMPLS
#'
#' Compute the prediction error for the envelope estimator using \code{kfold}
#' cross-validation.
#'
#' @param X,Y The data matrix.
#'
#' @return
#' \item{pred_error.inp}
#' \item{pred_error.IPE}
#' \item{pred_error.predenv}
#' \item{pred_error.PLS}
#' \item{pred_error.predenv2}
#' \item{pred_error.PLS2}
prederr.cv.sim<-function(X, Y, d.inp, d.IPE, u.predenv, u.PLS, kfold=3, nperm=50, NonGrass=F, PRINT=F){

# Test ---
#d.inp=d.IPE=1; u.predenv=u.PLS=6;
#NonGrass=F; PRINT=T
#kfold=5
#nperm=50
#--- ---
  Y<-as.matrix(Y)
  a<-dim(Y)
  n<-a[1]
  r<-a[2]
  p<-ncol(X)

  prederr.inp<-rep(0, nperm)
  prederr.IPE<-rep(0, nperm)
  prederr.predenv<-rep(0, nperm)
  prederr.PLS<-rep(0, nperm)
  prederr.predenv2<-rep(0, nperm)
  prederr.PLS2<-rep(0, nperm)

  for(i in 1:nperm){
    id<-sample(n, n)
    Xn<-X[id, ]
    Yn<-Y[id, ]
    Yn<-as.matrix(Yn)

    for(j in 1:kfold){
      id.test<-(floor((j-1)*n/kfold)+1):ceiling(j*n/kfold)
      id.train<-setdiff(1:n, id.test)
      X.train<-Xn[id.train, ]
      Y.train<-Yn[id.train, ]
      X.test<-Xn[id.test, ]
      Y.test<-Yn[id.test, ]
      n.test<-length(id.test)

      # inner PLS
      fit.inp<-innerPLS(X=X.train, Y=Y.train, d=d.inp)
      resi.inp<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.inp$alphainp)-as.matrix(X.test)%*%fit.inp$betainp)
      sprederr.inp<-apply(resi.inp, 1, function(x) sum(x^2))
      prederr.inp[i]<-prederr.inp[i]+sum(sprederr.inp)

      # IPE

      fit.IPE<-IPE(X=X.train, Y=Y.train, d=d.IPE, NonGrass=NonGrass)
      resi.IPE<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.IPE$alphaIPE)-as.matrix(X.test)%*%fit.IPE$betaIPE)
      sprederr.IPE<-apply(resi.IPE, 1, function(x) sum(x^2))
      prederr.IPE[i]<-prederr.IPE[i]+sum(sprederr.IPE)

      # Predictor Envelope

      fit.predenv<-Renvlp::xenv(X.train, Y.train, u=u.predenv, asy=F)
      resi.predenv<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.predenv$mu)-as.matrix(X.test)%*%fit.predenv$beta)
      sprederr.predenv<-apply(resi.predenv, 1, function(x) sum(x^2))
      prederr.predenv[i]<-prederr.predenv[i]+sum(sprederr.predenv)

      ## u=u.predenv+1
      fit.predenv2<-Renvlp::xenv(X.train, Y.train, u=u.predenv+1, asy=F)
      resi.predenv2<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.predenv2$mu)-as.matrix(X.test)%*%fit.predenv2$beta)
      sprederr.predenv2<-apply(resi.predenv2, 1, function(x) sum(x^2))
      prederr.predenv2[i]<-prederr.predenv2[i]+sum(sprederr.predenv2)

      # SIMPLS

      fit.PLS<-SIMPLS(X=X.train, Y=Y.train, u=u.PLS)
      resi.PLS<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.PLS$mu.PLS)-as.matrix(X.test)%*%fit.PLS$beta.PLS)
      sprederr.PLS<-apply(resi.PLS, 1, function(x) sum(x^2))
      prederr.PLS[i]<-prederr.PLS[i]+sum(sprederr.PLS)

      ## u=u.PLS+1
      fit.PLS2<-SIMPLS(X=X.train, Y=Y.train, u=u.PLS+1)
      resi.PLS2<-as.matrix(Y.test-matrix(1, n.test, 1)%*%t(fit.PLS2$mu.PLS)-as.matrix(X.test)%*%fit.PLS2$beta.PLS)
      sprederr.PLS2<-apply(resi.PLS2, 1, function(x) sum(x^2))
      prederr.PLS2[i]<-prederr.PLS2[i]+sum(sprederr.PLS2)

    }

    if(PRINT){
      cat("-")
    }
  }

  pred_error.inp<-mean(sqrt(prederr.inp/n))
  pred_error.IPE<-mean(sqrt(prederr.IPE/n))
  pred_error.predenv<-mean(sqrt(prederr.predenv/n))
  pred_error.PLS<-mean(sqrt(prederr.PLS/n))
  pred_error.predenv2<-mean(sqrt(prederr.predenv2/n))
  pred_error.PLS2<-mean(sqrt(prederr.PLS2/n))

  return(list(pred_error.inp=pred_error.inp, pred_error.IPE=pred_error.IPE,
              pred_error.predenv=pred_error.predenv, pred_error.PLS=pred_error.PLS,
              pred_error.predenv2=pred_error.predenv2, pred_error.PLS2=pred_error.PLS2))

}


#' Bootstrap Standard Errors for inner PLS
#'
#' Basically same to \code{\link{boot_res.inp}}, except that this is by
#' bootstrapping the original X and Y matrices.
#'
#' @return The bootstrap standard errors.
boot_dat.inp<-function(X, Y, d, B){

  #--- Test ---
  #  d<-1; B<-50
  #--- ---

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  dimY<-dim(Y)
  n<-dimY[1]
  r<-dimY[2]
  p<-ncol(X)

  bootinp<-function(i) {

    idx<-sample(1:n, n, replace=T)
    X.boot<-X[idx,]
    Y.boot<-Y[idx,]
    resu.inp<-innerPLS(X=X.boot, Y=Y.boot, d=d)

    return(resu.inp$betainp)
  }

  bootbeta<-lapply(1:B, function(i) bootinp(i))  # A list of B hatbetas.
  bootbeta<-matrix(unlist(bootbeta), nrow=B, byrow=T)
  bootse<-matrix(apply(bootbeta, 2, stats::sd), nrow=p)

  t.test_pval<-function(vec) stats::t.test(vec)$p.value
  bootpval<-matrix(apply(bootbeta, 2, t.test_pval), nrow=p)

  #test
  #  mean(bootbeta[,1])-sd(bootbeta[,1])/sqrt(B)*qt(p=0.025, df=B-1)
  #

  return(list(bootse=bootse, bootpval=bootpval))

}

#' Bootstrap Standard Errors for inner PLS
#'
#' Computes the bootstrap standard errors for each element of the regression
#' coefficients by bootstrapping the residuals.
#'
#' Learned from \code{\link[Renvlp]{boot.env}}.
#'
#' @param d Dimension of the inner envelope subspace. Can be 0? But may need
#' further debugging.
#'
#' @return The bootstrap standard errors.
boot_res.inp<-function(X, Y, d, B){

  #--- Test ---
  #  d<-1; B<-50
  #--- ---

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  dimY<-dim(Y)
  n<-dimY[1]
  r<-dimY[2]
  p<-ncol(X)

  resu.inp<-innerPLS(X=X, Y=Y, d=d)
  Yfit<-X%*%resu.inp$betainp
  res<-Y-Yfit

  bootinp<-function(i) {

    res.boot<-res[sample(1:n, n, replace=T),]
    Y.boot<-Yfit+res.boot
    resu.inp<-innerPLS(X=X, Y=Y.boot, d=d)

    return(resu.inp$betainp)
  }

  bootbeta<-lapply(1:B, function(i) bootinp(i))  # A list of B hatbetas.
  bootbeta<-matrix(unlist(bootbeta), nrow=B, byrow=T)
  bootse<-matrix(apply(bootbeta, 2, stats::sd), nrow=p)

  t.test_pval<-function(vec) stats::t.test(vec)$p.value
  bootpval<-matrix(apply(bootbeta, 2, t.test_pval), nrow=p)

  return(list(bootse=bootse, bootpval=bootpval))


}


#' Bootstrap Standard Errors for the inner PLS, IPE, Predictor envelopes and PLS.
#'
#' Computes the bootstrap standard errors for each element of the regression
#' coefficients by bootstrapping the original data matrices X and Y.
#'
#' @param d Dimension of the inner predictor envelope subspace.
#' @param u Dimension of the predictor envelope subspace.
#' @param u.PLS Dimension of the SIMPLS algorithm.
#' @param NonGrass Whether to use Non-Grassmann technique in \code{IPE} when
#' estimating the IPE.
#' @param ftol.nonG Used in \code{IPE}.
#'
#' @return
#' \item{bootse.inp}
#' \item{bootpval.inp}
#' \item{bootse.IPE}
#' \item{bootpval.IPE}
#' \item{bootse.predenv}
#' \item{bootpval.predenv}
#' \item{bootse.PLS}
#' \item{bootpval.PLS}
boot_dat<-function(X, Y, d, u, u.PLS, B, NonGrass=F, ftol.nonG=0.1, PRINT=F){

  #--- Test ---
#  d<-1; B<-50; u<-2
  #--- ---

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  dimY<-dim(Y)
  n<-dimY[1]
  r<-dimY[2]
  p<-ncol(X)

  idx_list<-lapply(1:B, function(i) sample(1:n, n, replace=T))
  t.test_pval<-function(vec) stats::t.test(vec)$p.value

  bootinp<-function(i) {
    return(innerPLS(X=X[idx_list[[i]],], Y=Y[idx_list[[i]],], d=d)$betaiPLS)
  }

  bootIPE<-function(i) {
    return(IPE(X=X[idx_list[[i]],], Y=Y[idx_list[[i]],], d=d, NonGrass=NonGrass, ftol.nonG=ftol.nonG)$betaIPE)
  }

  bootpredenv<-function(i) {
    return(Renvlp::xenv(X=X[idx_list[[i]],], Y=Y[idx_list[[i]],], u=u, asy=F)$beta)
  }

  bootPLS<-function(i) {
    SX.tmp<-stats::cov(X[idx_list[[i]],])*(n-1)/n
    SXY.tmp<-stats::cov(X[idx_list[[i]],], Y[idx_list[[i]],])*(n-1)/n
    Gam.pls.tmp<-SIMPLS_gamma(M=SX.tmp, U=tcrossprod(SXY.tmp), u=u.PLS)$Gammahat
    return(Gam.pls.tmp%*%chol2inv(chol(t(Gam.pls.tmp)%*%SX.tmp%*%Gam.pls.tmp))%*%t(Gam.pls.tmp)%*%SXY.tmp)
  }

  beta.inp<-innerPLS(X=X, Y=Y, d=d)$betaiPLS
  bootbeta.inp<-lapply(1:B, function(i) bootinp(i))  # A list of B hatbetas.
  bootbeta.inp<-matrix(unlist(bootbeta.inp), nrow=B, byrow=T)
  bootse.inp<-matrix(apply(bootbeta.inp, 2, stats::sd), nrow=p)
  bootpval.inp<-2*pnorm(abs(beta.inp)/bootse.inp, lower.tail=F)


  if(PRINT){
    cat("-")
  }

  beta.IPE<-IPE(X=X, Y=Y, d=d, NonGrass=NonGrass, ftol.nonG=ftol.nonG)$betaIPE
  bootbeta.IPE<-lapply(1:B, function(i) bootIPE(i))  # A list of B hatbetas.
  bootbeta.IPE<-matrix(unlist(bootbeta.IPE), nrow=B, byrow=T)
  bootse.IPE<-matrix(apply(bootbeta.IPE, 2, stats::sd), nrow=p)
  bootpval.IPE<-2*pnorm(abs(beta.IPE)/bootse.IPE, lower.tail=F)

  if(PRINT){
    cat("-")
  }

  beta.predenv<-Renvlp::xenv(X=X, Y=Y, u=u, asy=F)$beta
  bootbeta.predenv<-lapply(1:B, function(i) bootpredenv(i))  # A list of B hatbetas.
  bootbeta.predenv<-matrix(unlist(bootbeta.predenv), nrow=B, byrow=T)
  bootse.predenv<-matrix(apply(bootbeta.predenv, 2, stats::sd), nrow=p)
  bootpval.predenv<-2*pnorm(abs(beta.predenv)/bootse.predenv, lower.tail=F)

  if(PRINT){
    cat("-")
  }

  SX<-stats::cov(X)*(n-1)/n
  SXY<-stats::cov(X, Y)*(n-1)/n
  Gam.PLS<-SIMPLS_gamma(M=SX, U=tcrossprod(SXY), u=u.PLS)$Gammahat
  beta.PLS<-Gam.PLS%*%chol2inv(chol(t(Gam.PLS)%*%SX%*%Gam.PLS))%*%t(Gam.PLS)%*%SXY
  bootbeta.PLS<-lapply(1:B, function(i) bootPLS(i))  # A list of B hatbetas.
  bootbeta.PLS<-matrix(unlist(bootbeta.PLS), nrow=B, byrow=T)
  bootse.PLS<-matrix(apply(bootbeta.PLS, 2, stats::sd), nrow=p)
  bootpval.PLS<-2*pnorm(abs(beta.PLS)/bootse.PLS, lower.tail=F)


  return(list(bootse.inp=bootse.inp, bootpval.inp=bootpval.inp,
              bootse.IPE=bootse.IPE, bootpval.IPE=bootpval.IPE,
              bootse.predenv=bootse.predenv, bootpval.predenv=bootpval.predenv,
              bootse.PLS=bootse.PLS, bootpval.PLS=bootpval.PLS))
}


#' The iteration algorithm. Obsolete now.
#'
#' This algorithm is no good and is deprecated.
#'
#' Note that \eqn{\beta}(input and output) is the non-transposed version: same to response envelopes.
#'
#' @param ftol Should be less than 10.
#'
#' @return
#' \item{beta.iter}{}
#' \item{eta1}
#' \item{eta2}
#' \item{C}
iter_alg<-function(d, r, Gam1, Gam0, X, Y, ftol=1e-5, beta.t, PRINT){

  X1<-X%*%Gam1
  X2<-X%*%Gam0
  SYX<-stats::cov(Y, X)*(n-1)
  SXY<-t(SYX)
  SX<-stats::cov(X)*(n-1)/n
  SY<-stats::cov(Y)*(n-1)/n
  invSY<-chol2inv(chol(SY))

  # SV
  Omeo<-crossprod(Gam1, SX%*%Gam1)
  Omez<-crossprod(Gam0, SX%*%Gam0)
  eig.Omez<-eigen(Omez, symmetric=T)
  invOmez<-sweep(eig.Omez$vectors, MARGIN=2, 1/eig.Omez$values, "*")%*%t(eig.Omez$vectors)
  msqrOmez<-sweep(eig.Omez$vectors, MARGIN=2, 1/sqrt(eig.Omez$values), "*")%*%t(eig.Omez$vectors)

  term4C.span<-t(Gam0)%*%SXY%*%invSY%*%SYX%*%Gam0
  C.span<-msqrOmez%*%eigen(msqrOmez%*%term4C.span%*%msqrOmez, symmetric=T)$vectors[,1:(r-d)]
  C<-qr.Q(qr(C.span))

  eto<-chol2inv(chol(crossprod(Gam1, SX%*%Gam1)))%*%t(Gam1)%*%SXY
  ett<-chol2inv(chol(crossprod(Gam0%*%C, SX%*%Gam0%*%C)))%*%t(Gam0%*%C)%*%SXY

  # Iteration
  dis<-10

  if(PRINT){
    cat("## Iteration Algorithm ##\n\n")
  }

  while(dis>ftol){

    beta.old<-t(Gam1%*%eto+Gam0%*%C%*%ett)

    if(PRINT){
      cat(norm(beta.old-beta.t), "\n")
    }
    # R1
    R1<-Y-X2%*%C%*%ett
    SXR1<-stats::cov(X, R1)*(n-1)/n

    eto<-chol2inv(chol(crossprod(Gam1, SX%*%Gam1)))%*%t(Gam1)%*%SXR1

    # R2
    R2<-Y-X1%*%eto
    SXR2<-stats::cov(X, R2)*(n-1)/n
    SR2X<-t(SXR2)
    SR2<-stats::cov(R2)*(n-1)/n
    invSR2<-chol2inv(chol(SR2))

    term4C.span<-t(Gam0)%*%SXR2%*%invSR2%*%SR2X%*%Gam0
    C.span<-msqrOmez%*%eigen(msqrOmez%*%term4C.span%*%msqrOmez, symmetric=T)$vectors[,1:(r-d)]
    C<-qr.Q(qr(C.span))
    ett<-chol2inv(chol(crossprod(Gam0%*%C, SX%*%Gam0%*%C)))%*%t(Gam0%*%C)%*%SXR2

    beta.new<-t(Gam1%*%eto+Gam0%*%C%*%ett)

    dis<-norm(beta.old-beta.new)  # Criterion for the loop.
  }

  if(PRINT){
    cat("\n")
  }

  return(list(beta.iter=beta.new, eta1=eto, eta2=ett, C=C))
}
