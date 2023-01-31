#' Inner envelope estimation. USE ME!
#'
#' Estimating inner envelope subspace using non-Grassmann technique, which, on
#' on the same time, relies heavily on starting values. Provide MLE for other
#' parameters, especially \eqn{\beta_{IE}}, at the same time.
#'
#' The inverse in \code{fobj} is improved. \eqn{log|M|} is additionally added to
#'  \code{obj1, obj5}. \code{symmetric=F} is set in \code{eigen}.
#'
#' @param initial Starting values for \strong{Gamma0} and must be specified.
#' @param ftol Used in \code{abs(obj1-obj5)<ftol*abs(obj1)}.
#' @param PRINT Print the value of the objective function in each iteration?
#' @return \code{Gamma1, Gamma0, B, Omega0, objfun, betaIE}.
ienvMU<-function(M, U, d, initial, betaOLS=NULL, ftol=0.01, PRINT=F){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
  obj_added<-log(det(M))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))

  if(PRINT){
    cat("obj1: ", obj1, "\n")
  }


  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  #  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        invWmcm.M<-invW.M%*%tmp.M
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]
        #  H<-1+crossprod(tmp.M, invWmcm.M)*M[j,j]
        #  invT3.M<-invW.M-tcrossprod(invWmcm.M, invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M, invWmcm.M)*M[j,j])
        -2*log(1+x%*%(invcaocao%*%x))+log(1+M[j,j]*crossprod(tmp.M, invWmcm.M))+log(1+invM[j,j]*crossprod(tmp.invM, (invW.invM%*%tmp.invM)))+sum(log(1+eigen((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*%(invW.M-tcrossprod(invWmcm.M,invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j])),symmetric=F,only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj, method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))

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

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  invOmega0<-chol2inv(chol(Omega0))
  beta.ie<-(tcrossprod(Gam1, Gam1)+Gam0%*%B%*%chol2inv(chol(crossprod(B, (invOmega0%*%B))))%*%t(B)%*%tcrossprod(invOmega0, Gam0))%*%betaOLS

  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1, betaIE=beta.ie)
  )
}




#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv0}}.
#'
#' \code{fobj} is in its original form. \eqn{log|M|} is additionally added to
#' \code{obj1, obj5}.
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv0_logM<-function(M, U, d, initial){
  #while(h<=100){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
  obj_added<-log(det(M))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
  #  cat("obj1: ", obj1, "\n")

  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        T3.M<-W.M+tcrossprod(tmp.M, tmp.M)*M[j,j]
        T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]

        -2*log(1+x%*%(invcaocao%*%x))+log(1+M[j,j]*crossprod(tmp.M,(invW.M%*%tmp.M)))+log(1+invM[j,j]*crossprod(tmp.invM,(invW.invM%*%tmp.invM)))+sum(log(1+eigen((W.U+tcrossprod(tmp.U, tmp.U)*U[j,j])%*%chol2inv(chol(W.M+tcrossprod(tmp.M, tmp.M)*M[j,j])), only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj, method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
    #    cat("obj5: ", obj5, "\n")

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  #  cat("-")
  #  h<-h+1
  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}


#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv0}}.
#'
#' The inverse in \code{fobj} is improved.
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv1<-function(M, U, d, initial){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
  obj1<-sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
  #  cat("obj1: ", obj1, "\n")

  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        invWmcm.M<-invW.M%*%tmp.M
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]
        #  H<-1+crossprod(tmp.M, invWmcm.M)*M[j,j]
        #  invT3.M<-invW.M-tcrossprod(invWmcm.M, invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M, invWmcm.M)*M[j,j])
        -2*log(1+x%*%(invcaocao%*%x)) +log(1+M[j,j]*crossprod(tmp.M, invWmcm.M)) +log(1+invM[j,j]*crossprod(tmp.invM, (invW.invM%*%tmp.invM))) +sum(log(1+eigen((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*% (invW.M-tcrossprod(invWmcm.M,invWmcm.M)*M[j,j] /as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j])),only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj,  method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
    obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
#    cat("obj5: ", obj5, "\n")

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  #  cat("-")
  #  h<-h+1
  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}

#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv1}}.
#'
#' The inverse in \code{fobj} is improved. \eqn{log|M|} is additionally added to
#'  \code{obj1, obj5}.
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv1_logM<-function(M, U, d, initial){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
  obj_added<-log(det(M))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
  #  cat("obj1: ", obj1, "\n")

  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        invWmcm.M<-invW.M%*%tmp.M
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]
        #  H<-1+crossprod(tmp.M, invWmcm.M)*M[j,j]
        #  invT3.M<-invW.M-tcrossprod(invWmcm.M, invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M, invWmcm.M)*M[j,j])
        -2*log(1+x%*%(invcaocao%*%x)) +log(1+M[j,j]*crossprod(tmp.M, invWmcm.M)) +log(1+invM[j,j]*crossprod(tmp.invM, (invW.invM%*%tmp.invM))) +sum(log(1+eigen((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*% (invW.M-tcrossprod(invWmcm.M,invWmcm.M)*M[j,j] /as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j])),only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj, method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
    #    cat("obj5: ", obj5, "\n")

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}

#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv1}}.
#'
#' The inverse in \code{fobj} is improved. \eqn{log|M|} is additionally added to
#'  \code{obj1, obj5}. \code{symmetric=F} is set in \code{eigen}. This has been the
#' draft for my star function.
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv1_logM_sym<-function(M, U, d, initial, ftol=0.001, PRINT=F){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
  obj_added<-log(det(M))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))

  if(PRINT){
    cat("obj1: ", obj1, "\n")
  }


  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
#  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        invWmcm.M<-invW.M%*%tmp.M
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]
        #  H<-1+crossprod(tmp.M, invWmcm.M)*M[j,j]
        #  invT3.M<-invW.M-tcrossprod(invWmcm.M, invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M, invWmcm.M)*M[j,j])
        -2*log(1+x%*%(invcaocao%*%x))+log(1+M[j,j]*crossprod(tmp.M, invWmcm.M))+log(1+invM[j,j]*crossprod(tmp.invM, (invW.invM%*%tmp.invM)))+sum(log(1+eigen((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*%(invW.M-tcrossprod(invWmcm.M,invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j])),symmetric=F,only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj, method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))

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

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}


#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv1_logM_sym}}.
#'
#' Further improve \code{fobj}. Seems to add lots of brackets in fobj? I have
#' forgotten...
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv1_logM_sym_more<-function(M, U, d, initial){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
  obj_added<-log(det(M))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
  #  cat("obj1: ", obj1, "\n")

  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        invWmcm.M<-invW.M%*%tmp.M
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]
        #  H<-1+crossprod(tmp.M, invWmcm.M)*M[j,j]
        #  invT3.M<-invW.M-tcrossprod(invWmcm.M, invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M, invWmcm.M)*M[j,j])
        -2*log((1+x%*%(invcaocao%*%x)),base=exp(1))+ log((1+M[j,j]*crossprod(tmp.M, invWmcm.M)),base=exp(1))+ log((1+invM[j,j]*crossprod(tmp.invM, (invW.invM%*%tmp.invM))),base=exp(1))+ sum(log((1+eigen(((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*%(invW.M-(tcrossprod(invWmcm.M,invWmcm.M)*M[j,j])/as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j]))),symmetric=F,only.values=T)$values[(p-u1+1):(r-u1)]),base=exp(1)),na.rm=F)
      }

      res<-stats::optim(ca[j,], fobj, method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
    #    cat("obj5: ", obj5, "\n")

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}

#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv1}}.
#'
#' The inverse in \code{fobj} is improved. \eqn{log|M|} is additionally added to
#'  \code{obj1, obj5}. Moreover, \code{nlm} is used instead of \code{optim}.
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv1_logM_nlm<-function(M, U, d, initial){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
  obj_added<-log(det(M))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
  #  cat("obj1: ", obj1, "\n")

  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        invWmcm.M<-invW.M%*%tmp.M
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]
        #  H<-1+crossprod(tmp.M, invWmcm.M)*M[j,j]
        #  invT3.M<-invW.M-tcrossprod(invWmcm.M, invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M, invWmcm.M)*M[j,j])
        -2*log(1+x%*%(invcaocao%*%x)) +log(1+M[j,j]*crossprod(tmp.M, invWmcm.M)) +log(1+invM[j,j]*crossprod(tmp.invM, (invW.invM%*%tmp.invM))) +sum(log(1+eigen((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*% (invW.M-tcrossprod(invWmcm.M,invWmcm.M)*M[j,j] /as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j])),only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::nlm(fobj, ca[j,])
      ca[j,]<-res$estimate
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(t(Gam0)%*%M%*%Gam0, only.values=T)
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    e3<-eigen(crossprod(Gam0, (U%*%Gam0))%*%chol2inv(chol(crossprod(Gam0,  (M%*%Gam0)))), only.values=T)
    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
    #    cat("obj5: ", obj5, "\n")

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}


#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv0}}.
#'
#' The inverse in \code{fobj} is improved. And the symmetry in \code{obj1, obj5}
#' is improved.
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv2<-function(M, U, d, initial){
  #while(h<=100){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(crossprod(Gam0, (M%*%Gam0)))
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  Gm0resGm0.msqr<-sweep(e1$vectors, MARGIN=2, 1/sqrt(e1$values), "*")%*%t(e1$vectors)
  e3<-eigen(Gm0resGm0.msqr%*%crossprod(Gam0, (U%*%Gam0))%*%Gm0resGm0.msqr, symmetric=T, only.values=T)
  obj1<-sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
  #  cat("obj1: ", obj1, "\n")

  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        invWmcm.M<-invW.M%*%tmp.M
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]
        #  H<-1+crossprod(tmp.M, invWmcm.M)*M[j,j]
        #  invT3.M<-invW.M-tcrossprod(invWmcm.M, invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M, invWmcm.M)*M[j,j])
        -2*log(1+x%*%(invcaocao%*%x)) +log(1+M[j,j]*crossprod(tmp.M, invWmcm.M)) +log(1+invM[j,j]*crossprod(tmp.invM, (invW.invM%*%tmp.invM))) +sum(log(1+eigen((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*% (invW.M-tcrossprod(invWmcm.M,invWmcm.M)*M[j,j] /as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j])),only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj, method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(crossprod(Gam0, (M%*%Gam0)))
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    Gm0resGm0.msqr<-sweep(e1$vectors, MARGIN=2, 1/sqrt(e1$values), "*")%*%t(e1$vectors)
    e3<-eigen(Gm0resGm0.msqr%*%crossprod(Gam0, (U%*%Gam0))%*%Gm0resGm0.msqr, symmetric=T, only.values=T)
    obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
    #    cat("obj5: ", obj5, "\n")

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  #  cat("-")
  #  h<-h+1
  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}

#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv2}}.
#'
#' The inverse in \code{fobj} is improved. The symmetry in \code{obj1, obj5}
#' is improved. And \eqn{log|M|} is additionally added to \code{obj1, obj5}.
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv2_logM<-function(M, U, d, initial){
  #while(h<=100){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(crossprod(Gam0, (M%*%Gam0)))
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  Gm0resGm0.msqr<-sweep(e1$vectors, MARGIN=2, 1/sqrt(e1$values), "*")%*%t(e1$vectors)
  e3<-eigen(Gm0resGm0.msqr%*%crossprod(Gam0, (U%*%Gam0))%*%Gm0resGm0.msqr, symmetric=T, only.values=T)
  obj_added<-log(det(M))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
  #  cat("obj1: ", obj1, "\n")

  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        invWmcm.M<-invW.M%*%tmp.M
        #  T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]
        #  H<-1+crossprod(tmp.M, invWmcm.M)*M[j,j]
        #  invT3.M<-invW.M-tcrossprod(invWmcm.M, invWmcm.M)*M[j,j]/as.numeric(1+crossprod(tmp.M, invWmcm.M)*M[j,j])
        -2*log(1+x%*%(invcaocao%*%x)) +log(1+M[j,j]*crossprod(tmp.M, invWmcm.M)) +log(1+invM[j,j]*crossprod(tmp.invM, (invW.invM%*%tmp.invM))) +sum(log(1+eigen((W.U+tcrossprod(tmp.U,tmp.U)*U[j,j])%*% (invW.M-tcrossprod(invWmcm.M,invWmcm.M)*M[j,j] /as.numeric(1+crossprod(tmp.M,invWmcm.M)*M[j,j])),only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj,  method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(crossprod(Gam0, (M%*%Gam0)))
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    Gm0resGm0.msqr<-sweep(e1$vectors, MARGIN=2, 1/sqrt(e1$values), "*")%*%t(e1$vectors)
    e3<-eigen(Gm0resGm0.msqr%*%crossprod(Gam0, (U%*%Gam0))%*%Gm0resGm0.msqr, symmetric=T, only.values=T)
    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
    #    cat("obj5: ", obj5, "\n")

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  #  cat("-")
  #  h<-h+1
  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}


#' Inner envelope estimation.
#'
#' Revised version of \code{\link{ienv3}}.
#'
#' \eqn{log|M|} is additionally added to \code{obj1, obj5}.
#'
#' @param initial Starting values for Gamma0!
#' @return \code{Gamma1, Gamma0, B, Omega0}.
ienv3_logM<-function(M, U, d, initial){
  #while(h<=100){

  r<-dim(M)[1]
  u1<-d
  u0<-r-u1

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)

  Gam0<-initial

  e1<-eigen(crossprod(Gam0, (M%*%Gam0)))
  e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
  Gm0resGm0.msqr<-sweep(e1$vectors, MARGIN=2, 1/sqrt(e1$values), "*")%*%t(e1$vectors)
  e3<-eigen(Gm0resGm0.msqr%*%crossprod(Gam0, (U%*%Gam0))%*%Gm0resGm0.msqr, symmetric=T, only.values=T)
  obj_added<-log(det(M))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
  #  cat("obj1: ", obj1, "\n")

  GEidx.g0<-Renvlp::GE(Gam0)
  ca<-Gam0%*%solve(Gam0[GEidx.g0[1:u0],])  # This is the C_A without the rearrange of rows.
  caMca<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  cainvMca<-crossprod(ca, (invM%*%ca))  # C_A^T M^{-1} C_A
  caUca<-crossprod(ca, (U%*%ca))
  caca<-crossprod(ca[GEidx.g0[(u0+1):r],,drop=F], ca[GEidx.g0[(u0+1):r],,drop=F])+diag(u0) # C_A^T C_A

  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g0[(u0+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T
      mcm.M<-crossprod(ca[-j,], as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      mcm.invM<-crossprod(ca[-j,], as.matrix(invM[-j,j]))/invM[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      mcm.U<-crossprod(ca[-j,], as.matrix(U[-j,j]))/U[j,j]

      caocao<-caca-tcrossprod(a, a) # C_{A_1}^T C_{A_1}
      amcm.M<-a+mcm.M
      W.M<-caMca-tcrossprod(amcm.M, amcm.M)*M[j,j]  # W_1
      amcm.invM<-a+mcm.invM
      W.invM<-cainvMca-tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      W.U<-caUca-tcrossprod(amcm.U, amcm.U)*U[j,j]

      invW.M<-chol2inv(chol(W.M))
      invW.invM<-chol2inv(chol(W.invM))
      invcaocao<-chol2inv(chol(caocao))

      fobj<-function(x){
        tmp.M<-x+mcm.M
        tmp.invM<-x+mcm.invM
        tmp.U<-x+mcm.U
        #        T3.M<-W.M+tcrossprod(tmp.M, tmp.M)*M[j,j]
        eig.T3.M<-eigen(W.M+tcrossprod(tmp.M, tmp.M)*M[j,j])
        T3.M.msqr<-sweep(eig.T3.M$vectors, MARGIN=2, 1/sqrt(eig.T3.M$values), "*")%*%t(eig.T3.M$vectors)
        #        T3.U<-W.U+tcrossprod(tmp.U, tmp.U)*U[j,j]

        -2*log(1+x%*%(invcaocao%*%x)) +log(1+M[j,j]*crossprod(tmp.M,(invW.M%*%tmp.M))) +log(1+invM[j,j]*crossprod(tmp.invM,(invW.invM%*%tmp.invM)))+sum(log(1+eigen(T3.M.msqr%*%(W.U+tcrossprod(tmp.U, tmp.U)*U[j,j])%*%T3.M.msqr, symmetric=T, only.values=T)$values[(p-u1+1):(r-u1)]))
      }

      res<-stats::optim(ca[j,], fobj,  method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])
      caca<-caocao+tcrossprod(a, a)  # Return to C_A^T C_A
      amcm.M<-a+mcm.M
      caMca<-W.M+tcrossprod(amcm.M, amcm.M)*M[j,j]  # Get back to C_A^T M C_A
      amcm.invM<-a+mcm.invM
      cainvMca<-W.invM+tcrossprod(amcm.invM, amcm.invM)*invM[j,j]
      amcm.U<-a+mcm.U
      caUca<-W.U+tcrossprod(amcm.U, amcm.U)*U[j,j]

    }

    Gam0<-qr.Q(qr(ca))
    e1<-eigen(crossprod(Gam0, (M%*%Gam0)))
    e2<-eigen(t(Gam0)%*%invM%*%Gam0, only.values=T)
    Gm0resGm0.msqr<-sweep(e1$vectors, MARGIN=2, 1/sqrt(e1$values), "*")%*%t(e1$vectors)
    e3<-eigen(Gm0resGm0.msqr%*%crossprod(Gam0, (U%*%Gam0))%*%Gm0resGm0.msqr, symmetric=T, only.values=T)
    obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+e3$values[(p-u1+1):(r-u1)]))
    #    cat("obj5: ", obj5, "\n")

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }

  Gam1<-qr.Q(qr(Gam0), complete=T)[,(u0+1):r,drop=F]

  #---Omega_0---
  Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
  Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
  eig.GMG<-eigen(Gm0resGm0)
  Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
  Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
  eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

  V<-eig.boss$vectors
  K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
  Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

  #---B---
  eig.Omega0<-eigen(Omega0)
  Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
  #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
  B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

  #  cat("-")
  #  h<-h+1
  return(
    list(Gamma1=Gam1, Gamma0=Gam0, B=B, Omega0=Omega0, objfun=obj1)
  )
}


#' MLE of other parameters other than Gamma_1.
#'
#' @param M,U The \code{M} and \code{U} matrix.
#' @param u1 The dimension of \eqn{S_1} subspace. For historical reason, this is named
#' as u1, but not \code{d}. Should I change it?
#' @param p Must be inputted.
#'
#' @return
#' \item{betaIRE}
#' \item{Omega0.IRE}
#' \item{B.IRE}
mle_after_Gamma1.IRE<-function(M, U, Gam1, betaOLS, p, u1=dim(Gam1)[2], r=dim(Gam1)[1]){

  # Test---
#  M=dat$M; U=dat$U; Gam1=Gam1.pls; betaOLS=betaOLS; p=p
#  u1=dim(Gam1)[2]; r=dim(Gam1)[1]
  #--- ---

  if(u1==p){
    betaIRE<-tcrossprod(Gam1)%*%betaOLS
  } else{

    Gam0<-qr.Q(qr(Gam1), complete=T)[,(u1+1):r,drop=F]

    #---Omega0---
    Gm0fitGm0<-crossprod(Gam0, (U%*%Gam0))
    Gm0resGm0<-crossprod(Gam0, (M%*%Gam0))
    eig.GMG<-eigen(Gm0resGm0)
    Gm0resGm0.sqr<-sweep(eig.GMG$vectors, MARGIN=2, sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{\frac12}
    Gm0resGm0.msqr<-sweep(eig.GMG$vectors, MARGIN=2, 1/sqrt(eig.GMG$values), "*")%*%t(eig.GMG$vectors) # (Gam_0^T Sigma_{res} Gam_0)^{-\frac12}
    eig.boss<-eigen((Gm0resGm0.msqr%*%Gm0fitGm0%*%Gm0resGm0.msqr), symmetric=T)  # The eigen-structure of the most complex term. NOTICE: Set the symmetric arg to be TRUE! Or complex eigenvalues and eigenvectors may be produced.

    V<-eig.boss$vectors
    K<-c(rep(0,p-u1), eig.boss$values[(p-u1+1):(r-u1)])
    Omega0<-Gm0resGm0+Gm0resGm0.sqr%*%sweep(V, MARGIN=2, K, "*")%*%t(V)%*%Gm0resGm0.sqr

    #---B---
    eig.Omega0<-eigen(Omega0)
    Omega0.msqr<-sweep(eig.Omega0$vectors, MARGIN=2, 1/sqrt(eig.Omega0$values),"*")%*%t(eig.Omega0$vectors)
    #  B_unm<-Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)]
    B<-qr.Q(qr(Omega0%*%Omega0.msqr%*%eigen(Omega0.msqr%*%Gm0fitGm0%*%Omega0.msqr, symmetric=T)$vectors[,1:(p-u1)])) # Same as above to set symmetric=T.

    invOmega0<-chol2inv(chol(Omega0))
    betaIRE<-(tcrossprod(Gam1, Gam1)+Gam0%*%B%*%chol2inv(chol(crossprod(B, (invOmega0%*%B))))%*%t(B)%*%tcrossprod(invOmega0, Gam0))%*%betaOLS
  }

  return(list(betaIRE=betaIRE, Omega0.IRE=Omega0, B.IRE=B, Gam1.IRE=Gam1, Gam0.IRE=Gam0))
}



#' Calculate asymptotic variance of IPE for special case
#'
#' @param data.pop.IPE Object of class inheriting from \code{\link{data_gen.IPE}}.
#' Population version is fine.
#' @param n The supposed sample size.
#'
#' @return
#' \item{avarbetaIPE}
#' \item{avarbetaOLS}
#' \item{ratio_diag_OLS_IPE}
avar_sc.IPE<-function(data.pop.IPE, n=1000, d, p){

  Gam1<-data.pop.IPE$Gamma1
  Gam0<-data.pop.IPE$Gamma0
  C<-data.pop.IPE$C
  eta1<-eto<-data.pop.IPE$eta1
  eta2<-ett<-data.pop.IPE$eta2
  sigma2<-data.pop.IPE$eig.SiX[1]
  Lam<-diag(data.pop.IPE$eig.SiX[(d+1):p])  # Suppose that \sigma_0^2=1.
  SiYvX<-data.pop.IPE$SiYvX
  SiX<-data.pop.IPE$SiX
  #Omeo<-sigma2*diag(d)
  #Omez<-Lam

  # Preliminaries

  eig.Lam<-eigen(Lam)
  msqrLam<-sweep(eig.Lam$vectors, MARGIN=2, 1/sqrt(eig.Lam$values), "*")%*%t(eig.Lam$vectors)
  invLam<-sweep(eig.Lam$vectors, MARGIN=2, 1/eig.Lam$values, "*")%*%t(eig.Lam$vectors)
  sqrLam<-sweep(eig.Lam$vectors, MARGIN=2, sqrt(eig.Lam$values), "*")%*%t(eig.Lam$vectors)

  eig.SiYvX<-eigen(SiYvX, symmetric=T)
  invSiYvX<-sweep(eig.SiYvX$vectors, MARGIN=2, 1/eig.SiYvX$values, "*")%*%t(eig.SiYvX$vectors)
  sqrSiYvX<-sweep(eig.SiYvX$vectors, MARGIN=2, sqrt(eig.SiYvX$values), "*")%*%t(eig.SiYvX$vectors)
  msqrSiYvX<-sweep(eig.SiYvX$vectors, MARGIN=2, 1/sqrt(eig.SiYvX$values), "*")%*%t(eig.SiYvX$vectors)

  eig.SiX<-eigen(SiX, symmetric=T)
  invSiX<-sweep(eig.SiX$vectors, MARGIN=2, 1/eig.SiX$values, "*")%*%t(eig.SiX$vectors)

  PCLam<-Pr(C, prodmat=Lam)
  QCLam<-diag(p-d)-PCLam


  # Begin calculation

  M12<-1/sigma2*(tcrossprod(Gam1))%x%SiYvX+(Gam0%*%PCLam%*%invLam%*%t(Gam0))%x%SiYvX

  M3<-(Gam0%*%QCLam%*%invLam%*%t(Gam0))%x%(t(ett)%*%chol2inv(chol(ett%*%invSiYvX%*%t(ett)))%*%ett)

  term4M4<-Lam%*%QCLam+sigma2*invLam-2*diag(p-d)+1/sigma2*Lam
  M4<-(Gam0%*%QCLam%*%chol2inv(chol(term4M4))%*%t(QCLam)%*%t(Gam0))%x%crossprod(eto)

  avarbetaIPE<-(M12+M3+M4)/n

  # OLS
  avarbetaOLS<-invSiX%x%SiYvX/n

  # Compare the diagonal elements
  diag_sum_IPE<-sum(diag(avarbetaIPE))
  diag_sum_OLS<-sum(diag(avarbetaOLS))

  ratio_diag_OLS_IPE<-sum(diag(avarbetaOLS))/sum(diag(avarbetaIPE))

  return(list(avarbetaIPE=avarbetaIPE, avarbetaOLS=avarbetaOLS,
              diag_sum_IPE=diag_sum_IPE, diag_sum_OLS=diag_sum_OLS,
              ratio_diag_OLS_IPE=ratio_diag_OLS_IPE))
}

#' Calculate asymptotic variance of IPE
#'
#' See \code{\link{avar_sc.IPE}} for details.
avar_man.IPE<-function(data.pop.IPE, n=1000, d, r, p){

  Gam1<-data.pop.IPE$Gamma1
  Gam0<-data.pop.IPE$Gamma0
  C<-data.pop.IPE$C
  eto<-data.pop.IPE$eta1
  ett<-data.pop.IPE$eta2
  SiX<-data.pop.IPE$SiX
  SiYvX<-data.pop.IPE$SiYvX

  Omeo<-t(Gam1)%*%SiX%*%Gam1
  Omez<-t(Gam0)%*%SiX%*%Gam0
  invOmeo<-chol2inv(chol(Omeo))
  invOmez<-chol2inv(chol(Omez))
  invSiX<-chol2inv(chol(SiX))
  invSiYvX<-chol2inv(chol(SiYvX))
  CCOCC<-C%*%chol2inv(chol((t(C)%*%Omez%*%C)))%*%t(C)
  eeSee<-t(ett)%*%chol2inv(chol(ett%*%invSiYvX%*%t(ett)))%*%ett

  # The 1st&2nd terms, relatively simple
  M1<-(Gam1%*%invOmeo%*%t(Gam1))%x%SiYvX+(Gam0%*%CCOCC%*%t(Gam0))%x%SiYvX

  # Nightmare Term
  S1o<-eto%*%invSiYvX%*%t(eto)-eto%*%invSiYvX%*%eeSee%*%invSiYvX%*%t(eto)
  S1t<-Gam0%*%(Omez-Omez%*%CCOCC%*%Omez)%*%t(Gam0)
  S1<-S1o%x%S1t
  S2<-Omeo%x%(Gam0%*%invOmez%*%t(Gam0))-2*diag(d)%x%tcrossprod(Gam0)+invOmeo%x%(Gam0%*%Omez%*%t(Gam0))
  S<-S1+S2

  Krp<-matrixcalc::K.matrix(r,p)
  A2o<-eto%*%(diag(r)-invSiYvX%*%eeSee)
  A2t<-Gam0%*%(diag(p-d)-Omez%*%CCOCC)%*%t(Gam0)
  A2<-(A2o%x%A2t)%*%Krp

  M2<-t(A2)%*%MASS::ginv(S)%*%A2

  # Anyway, the last one
  M3<-(Gam0%*%(invOmez-CCOCC)%*%t(Gam0))%x%eeSee

  avar.man<-M1+M2+M3

  avarbetaIPE<-avar.man/n
  avarbetaOLS<-invSiX%x%SiYvX/n

  # Compare the diagonal elements
  diag_sum_IPE<-sum(diag(avarbetaIPE))
  diag_sum_OLS<-sum(diag(avarbetaOLS))

  ratio_diag_OLS_IPE<-sum(diag(avarbetaOLS))/sum(diag(avarbetaIPE))

  return(list(avarbetaIPE=avarbetaIPE, avarbetaOLS=avarbetaOLS,
              diag_sum_IPE=diag_sum_IPE, diag_sum_OLS=diag_sum_OLS,
              ratio_diag_OLS_IPE=ratio_diag_OLS_IPE))
}


#' Calculate asymptotic variance of IPE
#'
#' See \code{\link{avar_man.IPE}} for details. The only difference is that the data, \code{data.pop.IPE}
#' must have a \eqn{d \geq 2}.
#'
#' @param data.pop.IPE In the data, \eqn{d \geq 2} is required, for \code{matrixcalc::K.matrix}
#' to calculate the communication matrix.
avar_ori.IPE<-function(data.pop.IPE, n=1000, d, r, p){

  Gam1<-data.pop.IPE$Gamma1
  Gam0<-data.pop.IPE$Gamma0
  C<-data.pop.IPE$C
  eta1<-eto<-data.pop.IPE$eta1
  eta2<-ett<-data.pop.IPE$eta2
  SiX<-data.pop.IPE$SiX
  SiYvX<-data.pop.IPE$SiYvX
  Omeo<-t(Gam1)%*%SiX%*%Gam1
  Omez<-t(Gam0)%*%SiX%*%Gam0

  Kpd<-matrixcalc::K.matrix(p,d)
  Kpr<-matrixcalc::K.matrix(p,r)
  invSiX<-chol2inv(chol(SiX))
  invSiYvX<-chol2inv(chol(SiYvX))

  H11<-cbind(diag(r)%x%Gam1, diag(r)%x%(Gam0%*%C), t(eta2)%x%Gam0)
  H1<-rbind(Kpr%*%H11, matrix(0, nrow=p*(p+1)/2, ncol=r*d+r*(r-d)+(p-d)*(r-d)))  # Kpr in the left.

  H2<-rbind(Kpr%*%(t(eta1)%x%diag(p)-((t(eta2)%*%t(Gam0%*%C))%x%Gam1)%*%Kpd), 2*contr_rd(p)%*%((Gam1%*%Omeo)%x%diag(p)-Gam1%x%(Gam0%*%Omez%*%t(Gam0))))  # Kpr in the left

  H32<-cbind(contr_rd(p)%*%(Gam1%x%Gam1)%*%expan_rd(d), contr_rd(p)%*%(Gam0%x%Gam0)%*%expan_rd(p-d))
  H3<-rbind(matrix(0, nrow=r*p, ncol=d*(d+1)/2+(p-d)*(p-d+1)/2), H32)
  H<-cbind(H1, H2, H3)

  J_row1<-cbind(SiX%x%invSiYvX, matrix(0, ncol=p*(p+1)/2, nrow=p*r))
  J_row2<-cbind(matrix(0, nrow=p*(p+1)/2, ncol=p*r),  1/2*t(expan_rd(p))%*%(invSiX%x%invSiX)%*%expan_rd(p))
  J<-rbind(J_row1, J_row2)

  avar.H<-(H%*%MASS::ginv(crossprod(H, J%*%H))%*%t(H))[1:(r*p),1:(r*p)]

  avarbetaIPE<-avar.H/n
  avarbetaOLS<-invSiX%x%SiYvX/n

  # Compare the diagonal elements
  diag_sum_IPE<-sum(diag(avarbetaIPE))
  diag_sum_OLS<-sum(diag(avarbetaOLS))

  ratio_diag_OLS_IPE<-sum(diag(avarbetaOLS))/sum(diag(avarbetaIPE))

  return(list(avarbetaIPE=avarbetaIPE, avarbetaOLS=avarbetaOLS,
              diag_sum_IPE=diag_sum_IPE, diag_sum_OLS=diag_sum_OLS,
              ratio_diag_OLS_IPE=ratio_diag_OLS_IPE))
}
