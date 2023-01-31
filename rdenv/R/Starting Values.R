#' Generic starting values computation function.
#'
#' This function is completely copied from \code{\link[Renvlp]{envMU}} only with
#' the names of the variable revised.
#'
#' @param U U is not of full rank, or it follows that u2=0.
#' @param u u<r.
#' @param Gam.t A redundant input, to make inputs conform to that of \code{\link{stv2}}.
#' @return The starting values and value of the minimized objective function:
#' \code{list(init=init1, obj=obj1)}.
stv<-function(M, U, u, Gam.t=NULL){

  #---M+U---

  #---Normal Version of U---

  MU<-M+U
  eig.MU<-eigen(MU)
  invMU<-sweep(eig.MU$vectors, MARGIN=2, 1/eig.MU$values, "*")%*%t(eig.MU$vectors)  # Get the inverse of MU.
  invMU2<-sweep(eig.MU$vectors, MARGIN=2, 1/sqrt(eig.MU$values), "*")%*%t(eig.MU$vectors)  # Inverse of the square root matrix

  startv<-function(g, midmatrix){
    t(g)%*%midmatrix%*%g
  }

  J1<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, U)
  J1.st<-sort(J1, decreasing=T, index.return=T)
  init1<-as.matrix(eig.MU$vectors[,J1.st$ix[1:u]])  # Choose u eigenvectors that maximize J_M^*(G). We get the starting value here.
  eig1<-eigen(t(init1)%*%M%*%init1)
  eig2<-eigen(t(init1)%*%invMU%*%init1)
  obj1<-sum(log(eig1$values))+sum(log(eig2$values))  # Here is L_u(G).

  #---Standardized version of U---
  U.stdMU<-invMU2%*%tcrossprod(U, invMU2)
  J2<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, U.stdMU)
  J2.st<-sort(J2, decreasing=T, index.return=T)
  init2<-as.matrix(eig.MU$vectors[,J2.st$ix[1:u]])
  e1<-eigen(t(init2)%*%M%*%init2)
  e2<-eigen(t(init2)%*%invMU%*%init2)
  obj2<-sum(log(e1$values))+sum(log(e2$values))  # Use the product of eigenvalues facilitate the calculation of determinant.

  if(obj2 < obj1){
    init1<-init2
    obj1<-obj2
  }

  #---M---

  eig.M<-eigen(M)
  J3<-apply(X=eig.M$vectors, MARGIN=2, FUN=startv, U)
  J3.st<-sort(J3, decreasing=T, index.return=T)
  init3<-as.matrix(eig.M$vectors[,J3.st$ix[1:u]])
  e1<-eigen(t(init3)%*%M%*%init3)
  e2<-eigen(t(init3)%*%invMU%*%init3)
  obj3<-sum(log(e1$values))+sum(log(e2$values))

  if(obj3 < obj1){
    init1<-init3
    obj1<-obj3
  }

  invM2<-sweep(eig.M$vectors, MARGIN=2, 1/sqrt(eig.M$values), "*")%*%t(eig.M$vectors)

  U.stdM<-invM2%*%tcrossprod(U, invM2)
  J4<-apply(X=eig.M$vectors, MARGIN=2, FUN=startv, U.stdM)
  J4.st<-sort(J4, decreasing=T, index.return=T)
  init4<-as.matrix(eig.M$vectors[,J4.st$ix[1:u]])
  e1<-eigen(t(init4)%*%M%*%init4)
  e2<-eigen(t(init4)%*%invMU%*%init4)
  obj4<-sum(log(e1$values))+sum(log(e2$values))

  if(obj4 < obj1){
    init1<-init4
    obj1<-obj4
  }

  return(list(init=init1, obj=obj1))

}



#' Revised version of \code{\link{stv}}.
#'
#' This function is for studying and testing.
#'
#' The 4 potential starting values are:
#' \itemize{
#'   \item Normal U by eigenvectors of \code{M+U},
#'   \item Standardized U by eigenvectors of \code{M+U},
#'   \item Normal U by eigenvectors of \code{M},
#'   \item Standardized U by eigenvectors of \code{M}.
#' }
#' We can choose to print the distance between the estimator and the true
#' parameter and the values of the objective function of each of the 4 sets
#' of choices.
#'
#' Moreover, the 4 candidates, the values of the objective function on them and
#'  the indices will all be returned.
#'
#' Finally, two other more candidates are added if \code{c} is not \code{NULL},
#'  for studying.
#'
#' @param c,p Scale added to \eqn{P_U} and the dimension of predictors, which is
#' only in need when specifying c.
#' @return \code{list(init=init_min, obj=obj_min, init1=init1, init2=init2, init3=init3,
#' init4=init4, obj1=obj1, obj2=obj2, obj3=obj3, obj4=obj4,
#' idx1=idx1, idx2=idx2, idx3=idx3, idx4=idx4)}
stv2<-function(M, U, u, Gam.t, PRINT=F, c=NULL, p=NULL){

  r<-dim(M)[1]

  if(PRINT){
    cat("## 4sv ##\n\n")
    #---M+U---
    cat("Eigenvectors from M+U.\n")

    #---Normal Version of U---
    cat("Normal U.\n")
  }

  MU<-M+U
  eig.MU<-eigen(MU)
  invMU<-sweep(eig.MU$vectors, MARGIN=2, 1/eig.MU$values, "*")%*%t(eig.MU$vectors)  # Get the inverse of MU.
  invMU2<-sweep(eig.MU$vectors, MARGIN=2, 1/sqrt(eig.MU$values), "*")%*%t(eig.MU$vectors)  # Inverse of the square root matrix

  startv<-function(g, midmatrix){
    t(g)%*%midmatrix%*%g
  }

  J1<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, U)
  J1.st<-sort(J1, decreasing=T, index.return=T)
  init1<-as.matrix(eig.MU$vectors[,J1.st$ix[1:u]])  # Choose u eigenvectors that maximize J_M^*(G). We get the starting value here.
  eig1<-eigen(t(init1)%*%M%*%init1)
  eig2<-eigen(t(init1)%*%invMU%*%init1)
  obj1<-sum(log(eig1$values))+sum(log(eig2$values))  # Here is L_u(G).
  idx1<-J1.st$ix[1:u]

  if(PRINT){
    cat(frbnorm(init1, Gam.t), "|", angle(init1, Gam.t), "|", obj1, "\n")
  }

  init_min<-init1
  obj_min<-obj1

  #---Standardized version of U---
  if(PRINT){
    cat("Standardized U.\n")
  }

  U.stdMU<-invMU2%*%tcrossprod(U, invMU2)
  J2<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, U.stdMU)
  J2.st<-sort(J2, decreasing=T, index.return=T)
  init2<-as.matrix(eig.MU$vectors[,J2.st$ix[1:u]])
  e1<-eigen(t(init2)%*%M%*%init2)
  e2<-eigen(t(init2)%*%invMU%*%init2)
  obj2<-sum(log(e1$values))+sum(log(e2$values))
  idx2<-J2.st$ix[1:u]

  if(PRINT){
    cat(frbnorm(init2, Gam.t), "|", angle(init2, Gam.t), "|", obj2, "\n")
  }

  if(obj2<obj_min){
    init_min<-init2
    obj_min<-obj2
  }

  if(PRINT){
    #---M---
    cat("Eigenvectors from M.\n")

    #---Normal U---
    cat("Normal U.\n")
  }


  eig.M<-eigen(M)
  J3<-apply(X=eig.M$vectors, MARGIN=2, FUN=startv, U)
  J3.st<-sort(J3, decreasing=T, index.return=T)
  init3<-as.matrix(eig.M$vectors[,J3.st$ix[1:u]])
  e1<-eigen(t(init3)%*%M%*%init3)
  e2<-eigen(t(init3)%*%invMU%*%init3)
  obj3<-sum(log(e1$values))+sum(log(e2$values))
  idx3<-J3.st$ix[1:u]

  if(PRINT){
    cat(frbnorm(init3, Gam.t), "|", angle(init3, Gam.t), "|", obj3, "\n")
  }

  if(obj3<obj_min){
    init_min<-init3
    obj_min<-obj3
  }

  #---Standardized U---
  if(PRINT){
    cat("Standardized U.\n")
  }

  invM2<-sweep(eig.M$vectors, MARGIN=2, 1/sqrt(eig.M$values), "*")%*%t(eig.M$vectors)

  U.stdM<-invM2%*%tcrossprod(U, invM2)
  J4<-apply(X=eig.M$vectors, MARGIN=2, FUN=startv, U.stdM)
  J4.st<-sort(J4, decreasing=T, index.return=T)
  init4<-as.matrix(eig.M$vectors[,J4.st$ix[1:u]])
  e1<-eigen(t(init4)%*%M%*%init4)
  e2<-eigen(t(init4)%*%invMU%*%init4)
  obj4<-sum(log(e1$values))+sum(log(e2$values))
  idx4<-J4.st$ix[1:u]

  if(PRINT){
    cat(frbnorm(init4, Gam.t), "|", angle(init4, Gam.t), "|", obj4, "\n")
  }

  if(obj4<obj_min){
    init_min<-init4
    obj_min<-obj4
  }

  #---Use eigenvectors from M+U and Pu for inner product.

#  if(PRINT){
#    cat("Eigenvectors from M+U and P_U for inner product.\n")
#  }

#  J7<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, Pr(U))
#  J7.st<-sort(J7, decreasing=T, index.return=T)
#  init7<-as.matrix(eig.MU$vectors[,J7.st$ix[1:u]])  # Choose u eigenvectors that maximize J_M^*(G). We get the starting value here.
#  eig1<-eigen(t(init7)%*%M%*%init7)
#  eig2<-eigen(t(init7)%*%invMU%*%init7)
#  obj7<-sum(log(eig1$values))+sum(log(eig2$values))  # Here is L_u(G).
#  idx7<-J7.st$ix[1:u]

#  if(PRINT){
#    cat(frbnorm(init7, Gam.t), "|", angle(init7, Gam.t), "|", obj7, "\n")
#  }

  #---A New Idea by scaling Pu---
  if(!is.null(c)){

    if(is.null(p)){
      stop("Input p!\n")
    }
    if(PRINT){

      cat("Eigenvectors from M+c*(min(u,p) eigvec.U)\n")

      cat("Normal U.\n")
    }

    eig.U<-eigen(U)
    MU<-M+c*tcrossprod(eig.U$vectors[,1:min(u,p)], NULL)
    eig.MU<-eigen(MU)
    invMU2<-sweep(eig.MU$vectors, MARGIN=2, 1/sqrt(eig.MU$values), "*")%*%t(eig.MU$vectors)  # Inverse of the square root matrix

    J5<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, U)
    J5.st<-sort(J1, decreasing=T, index.return=T)
    init5<-as.matrix(eig.MU$vectors[,J5.st$ix[1:u]])
    eig1<-eigen(t(init5)%*%M%*%init5)
    eig2<-eigen(t(init5)%*%invMU%*%init5)
    obj5<-sum(log(Re(eig1$values)))+sum(log(Re(eig2$values)))
    idx5<-J5.st$ix[1:u]

    if(PRINT){
      cat(frbnorm(init5, Gam.t), "|", angle(init5, Gam.t), "|", obj5, "\n")
    }

    if(obj5<obj_min){
      init_min<-init5
      obj_min<-obj5
    }

    #---Standardized version of U---
    if(PRINT){
      cat("Standardized U.\n")
    }

    U.stdMU<-invMU2%*%tcrossprod(U, invMU2)
    J6<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, U.stdMU)
    J6.st<-sort(J6, decreasing=T, index.return=T)
    init6<-as.matrix(eig.MU$vectors[,J6.st$ix[1:u]])
    e1<-eigen(t(init6)%*%M%*%init6)
    e2<-eigen(t(init6)%*%invMU%*%init6)
    obj6<-sum(log(Re(e1$values)))+sum(log(Re(e2$values)))
    idx6<-J6.st$ix[1:u]

    if(PRINT){
      cat(frbnorm(init6, Gam.t), "|", angle(init6, Gam.t), "|", obj6, "\n")
      cat("\n")
    }

    if(obj6<obj_min){
      init_min<-init6
      obj_min<-obj6
    }

  }

  return(list(init=init_min, obj=obj_min, init1=init1, init2=init2, init3=init3,
              init4=init4, obj1=obj1, obj2=obj2, obj3=obj3, obj4=obj4,
              idx1=idx1, idx2=idx2, idx3=idx3, idx4=idx4))
}


#' Improved version of \code{\link{stv}}
#'
#' Which of four candidates to use can be specified.
#'
#' The 4 potential starting values are:
#' \itemize{
#'   \item Normal U by eigenvectors of \code{M+U},
#'   \item Standardized U by eigenvectors of \code{M+U},
#'   \item Normal U by eigenvectors of \code{M},
#'   \item Standardized U by eigenvectors of \code{M}.
#' }
#'
#' @param USE A vector, each of whose element taking value from 1 to 4. Specify
#' which candidates to use.
stv3<-function(M, U, u, USE=1:4, Gam.t, PRINT=F){

  if(is.null(USE)){
    stop("Why using this function???")
  }
  if(any(USE>4) || any(USE<1)){
    stop("Wrong candidates specification!")
  }

  if(PRINT){
    cat("## 4sv ##\n\n")
  }

  r<-dim(M)[1]
  MU<-M+U
  eig.MU<-eigen(MU)
  eig.M<-eigen(M)
  invMU<-sweep(eig.MU$vectors, MARGIN=2, 1/eig.MU$values, "*")%*%t(eig.MU$vectors)  # Get the inverse of MU.
  invMU2<-sweep(eig.MU$vectors, MARGIN=2, 1/sqrt(eig.MU$values), "*")%*%t(eig.MU$vectors)  # Inverse of the square root matrix

  startv<-function(g, midmatrix){
    t(g)%*%midmatrix%*%g
  }

  init_cand<-list(NULL)
  obj_cand<-numeric(0)
  i<-1

  if(1 %in% USE){

    if(PRINT){
      cat("Eigenvectors from M+U.\n")
      cat("Normal U.\n")
    }

    J1<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, U)
    J1.st<-sort(J1, decreasing=T, index.return=T)
    init1<-as.matrix(eig.MU$vectors[,J1.st$ix[1:u]])  # Choose u eigenvectors that maximize J_M^*(G). We get the starting value here.
    eig1<-eigen(t(init1)%*%M%*%init1)
    eig2<-eigen(t(init1)%*%invMU%*%init1)
    obj1<-sum(log(eig1$values))+sum(log(eig2$values))  # Here is L_u(G).
    idx1<-J1.st$ix[1:u]

    init_cand[[i]]<-init1
    obj_cand[i]<-obj1
    i<-i+1

    if(PRINT){
      cat(frbnorm(init1, Gam.t), "|", angle(init1, Gam.t), "|", obj1, "\n")
    }
  }

  if(2 %in% USE){
    if(PRINT){
      cat("Eigenvectors from M+U.\n")
      cat("Standardized U.\n")
    }

    U.stdMU<-invMU2%*%tcrossprod(U, invMU2)
    J2<-apply(X=eig.MU$vectors, MARGIN=2, FUN=startv, U.stdMU)
    J2.st<-sort(J2, decreasing=T, index.return=T)
    init2<-as.matrix(eig.MU$vectors[,J2.st$ix[1:u]])
    e1<-eigen(t(init2)%*%M%*%init2)
    e2<-eigen(t(init2)%*%invMU%*%init2)
    obj2<-sum(log(e1$values))+sum(log(e2$values))
    idx2<-J2.st$ix[1:u]

    init_cand[[i]]<-init2
    obj_cand[i]<-obj2
    i<-i+1

    if(PRINT){
      cat(frbnorm(init2, Gam.t), "|", angle(init2, Gam.t), "|", obj2, "\n")
    }
  }

  if(3 %in% USE){
    if(PRINT){
      cat("Eigenvectors from M.\n")
      cat("Normal U.\n")
    }

    J3<-apply(X=eig.M$vectors, MARGIN=2, FUN=startv, U)
    J3.st<-sort(J3, decreasing=T, index.return=T)
    init3<-as.matrix(eig.M$vectors[,J3.st$ix[1:u]])
    e1<-eigen(t(init3)%*%M%*%init3)
    e2<-eigen(t(init3)%*%invMU%*%init3)
    obj3<-sum(log(e1$values))+sum(log(e2$values))
    idx3<-J3.st$ix[1:u]

    init_cand[[i]]<-init3
    obj_cand[i]<-obj3
    i<-i+1

    if(PRINT){
      cat(frbnorm(init3, Gam.t), "|", angle(init3, Gam.t), "|", obj3, "\n")
    }
  }

  if(4 %in% USE){

    if(PRINT){
      cat("Eigenvectors from M.\n")
      cat("Standardized U.\n")
    }

    invM2<-sweep(eig.M$vectors, MARGIN=2, 1/sqrt(eig.M$values), "*")%*%t(eig.M$vectors)

    U.stdM<-invM2%*%tcrossprod(U, invM2)
    J4<-apply(X=eig.M$vectors, MARGIN=2, FUN=startv, U.stdM)
    J4.st<-sort(J4, decreasing=T, index.return=T)
    init4<-as.matrix(eig.M$vectors[,J4.st$ix[1:u]])
    e1<-eigen(t(init4)%*%M%*%init4)
    e2<-eigen(t(init4)%*%invMU%*%init4)
    obj4<-sum(log(e1$values))+sum(log(e2$values))
    idx4<-J4.st$ix[1:u]

    init_cand[[i]]<-init4
    obj_cand[i]<-obj4

    if(PRINT){
      cat(frbnorm(init4, Gam.t), "|", angle(init4, Gam.t), "|", obj4, "\n")
    }
  }

  min_idx<-which.min(obj_cand)
  init_min<-init_cand[[min_idx]]
  obj_min<-obj_cand[min_idx]

  return(list(init=init_min, obj=obj_min))
}



#' Starting values for inner envelope subspace estimation!
#'
#' A fairly good algorithm by choosing from 4 candidates.
#'
#' When comparing the values of objective functions at the candidate points,
#' module is used as \code{eigen} always returns complex eigenvalues.
#'
#' The four candidates are:
#' \itemize{
#'   \item Standardized U by eigenvectors of \code{M},
#'   \item Normal U by eigenvectors of \code{M},
#'   \item Standardized U by eigenvectors of \code{M+Q_U},
#'   \item Normal U by eigenvectors of \code{M+Q_U}.
#' }
#'
#' @param d The dimension of the inner envelope subspace.
#' @param Gam0.t True value of Gamma0 as to print the process when \code{PRINT=TRUE}.
#' @param PRINT Logical. Whether to print? Will print the Frobenius distance
#' between each of the starting values candidates and the true value, and the
#' value of the objective function with \code{det(M)} added.
#' @return The starting values for Gamma0, value of the MLE of inner envelope
#'  at the chosen starting point and all the 4 candidates:
#'  \code{list(init=init_min, obj=obj_min, init1=init1, init2=init2, init3=init3,
#' init4=init4, init5=init5))}.
init.ienv<-function(M, U, d, p, Gam0.t=NULL, PRINT=F){

  r<-dim(M)[1]

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)
  msqrM<-sweep(eig.M$vectors, MARGIN=2, 1/sqrt(eig.M$values), "*")%*%t(eig.M$vectors)

  eig.U<-eigen(U, symmetric=T, only.values=F)
  QU<-tcrossprod(eig.U$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of U.
  eig.M_QU<-eigen(M+QU)
  invM_QU<-sweep(eig.M_QU$vectors, MARGIN=2, 1/eig.M_QU$values, "*")%*%t(eig.M_QU$vectors)
  msqrM_QU<-sweep(eig.M_QU$vectors, MARGIN=2, 1/sqrt(eig.M_QU$values), "*")%*%t(eig.M_QU$vectors)

  #---Use eigenvectors from M---
  if(PRINT){
    cat("\n## My New Method ##\n\n")
    cat("--Eigenvectors from M--\n")
  }
  #---First one: M^{-\frac12} U M^{-\frac12}---

  eig.final1<-eigen((msqrM%*%U%*%msqrM), symmetric=T, only.values=F)
  Pro_small1<-tcrossprod(eig.final1$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of M^{-\frac12} U M^{-\frac12}.

  # We need to find the r-d eigenvectors of M that is not orthogonal to the r-p smallest eigenvectors of M^{-\frac12} U M^{-\frac12} or U
  # So we use "angle" between each eigenvector of M and the space spanned by eig.final$vectors[,(p+1):r].
  # The "angle" is computed as the module of the projected eigenvector.
  ang1<-sqrt(colSums((Pro_small1%*%eig.M$vectors)^2))
  idx1<-sort(ang1, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init1<-as.matrix(eig.M$vectors[,idx1])

  e1<-eigen(t(init1)%*%M%*%init1, symmetric=T, only.values=T)
  e2<-eigen(t(init1)%*%invM%*%init1, only.values=T)
  e3<-eigen(crossprod(init1, (U%*%init1))%*%chol2inv(chol(crossprod(init1, (M%*%init1)))), symmetric=F, only.values=T)
  obj_added<-log(prod(eig.M$values))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  e2.env<-eigen(t(init1)%*%invM_QU%*%init1, symmetric=T, only.values=T)
  obj1.env<-sum(log(e1$values))+sum(log(e2.env$values))


  if(PRINT){
    cat("Standardized U.\n")
    cat(frbnorm(init1, Gam0.t), "|", angle(init1, Gam0.t), "|", obj1, "|", obj1.env, "\n")
  }

  init_min<-init.env_min<-init1
  obj_min<-obj1
  obj.env_min<-obj1.env
  idx<-idx1

  #---Second one: U---

  ang2<-sqrt(colSums((QU%*%eig.M$vectors)^2))
  idx2<-sort(ang2, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init2<-as.matrix(eig.M$vectors[,idx2])

  e1<-eigen(t(init2)%*%M%*%init2, symmetric=T, only.values=T)
  e2<-eigen(t(init2)%*%invM%*%init2, only.values=T)
  e3<-eigen(crossprod(init2, (U%*%init2))%*%chol2inv(chol(crossprod(init2, (M%*%init2)))), symmetric=F, only.values=T)
  obj2<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  e2.env<-eigen(t(init2)%*%invM_QU%*%init2, symmetric=T, only.values=T)
  obj2.env<-sum(log(e1$values))+sum(log(e2.env$values))

  if(PRINT){
    cat("Normal U.\n")
    cat(frbnorm(init2, Gam0.t), "|", angle(init2, Gam0.t), "|", obj2, "|", obj2.env, "\n")
  }

  if(obj2<obj_min){
    init_min<-init2
    obj_min<-obj2
    idx<-idx2
  }
  if(obj2.env<obj.env_min){
    init.env_min<-init2
    obj.env_min<-obj2.env
  }


  #---The fifth(Yes, no problem) candidate: M^{-\frac12} Q_{U} M^{-\frac12} by eigvec.M---

#  msqrM_QU_msqrM<-msqrM%*%QU%*%msqrM
#  startv<-function(g, midmatrix){
#    t(g)%*%midmatrix%*%g
#  }
#  J5<-apply(X=eig.M$vectors, MARGIN=2, FUN=startv, msqrM_QU_msqrM)
#  idx5<-sort(J5, decreasing=T, index.return=T)$ix[1:(r-d)]
#  init5<-as.matrix(eig.M$vectors[,idx5])


#  e1<-eigen(t(init5)%*%M%*%init5, symmetric=T, only.values=T)
#  e2<-eigen(t(init5)%*%invM%*%init5, only.values=T)
#  e3<-eigen(crossprod(init5, (U%*%init5))%*%chol2inv(chol(crossprod(init5, (M%*%init5)))), symmetric=F, only.values=T)
#  obj5<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

#  if(PRINT){
#    cat("msqrM_QU_msqrM.\n")
#    cat(frbnorm(init5, Gam0.t), "|", angle(init5, Gam0.t), "|", obj5, "\n")
#  }

  #---Use eigenvectors from M+Q_U---

  if(PRINT){
    cat("--Eigenvectors from M+Q_U--\n")
  }

  #---Third one: {M+Q_U}^{-\frac12} U {M+Q_U}^{-\frac12}---

  eig.final3<-eigen((msqrM_QU%*%U%*%msqrM_QU), symmetric=T, only.values=F)
  Pro_small3<-tcrossprod(eig.final3$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of M^{-\frac12} U M^{-\frac12}.

  ang3<-sqrt(colSums((Pro_small3%*%eig.M_QU$vectors)^2))
  idx3<-sort(ang3, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init3<-as.matrix(eig.M_QU$vectors[,idx3])

  e1<-eigen(t(init3)%*%M%*%init3, symmetric=T, only.values=T)
  e2<-eigen(t(init3)%*%invM%*%init3, only.values=T)
  e3<-eigen(crossprod(init3, (U%*%init3))%*%chol2inv(chol(crossprod(init3, (M%*%init3)))), symmetric=F, only.values=T)
  obj3<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  e2.env<-eigen(t(init3)%*%invM_QU%*%init3, symmetric=T, only.values=T)
  obj3.env<-sum(log(e1$values))+sum(log(e2.env$values))

  if(PRINT){
    cat("Standardized U.\n")
    cat(frbnorm(init3, Gam0.t), "|", angle(init3, Gam0.t), "|", obj3, "|", obj3.env, "\n")
  }

  if(obj3<obj_min){
    init_min<-init3
    obj_min<-obj3
    idx<-idx3
  }
  if(obj3.env<obj.env_min){
    init.env_min<-init3
    obj.env_min<-obj3.env
  }

  #---The fourth candidate: U---

  ang4<-sqrt(colSums((QU%*%eig.M_QU$vectors)^2))
  idx4<-sort(ang4, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init4<-as.matrix(eig.M_QU$vectors[,idx4])

  e1<-eigen(t(init4)%*%M%*%init4, symmetric=T, only.values=T)
  e2<-eigen(t(init4)%*%invM%*%init4, only.values=T)
  e3<-eigen(crossprod(init4, (U%*%init4))%*%chol2inv(chol(crossprod(init4, (M%*%init4)))), symmetric=F, only.values=T)
  obj4<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  e2.env<-eigen(t(init4)%*%invM_QU%*%init4, symmetric=T, only.values=T)
  obj4.env<-sum(log(e1$values))+sum(log(e2.env$values))

  if(PRINT){
    cat("Normal U.\n")
    cat(frbnorm(init4, Gam0.t), "|", angle(init4, Gam0.t), "|", obj4, "|", obj4.env, "\n")
    cat("\n")
  }

  if(obj4<obj_min){
    init_min<-init4
    obj_min<-obj4
    idx<-idx4
  }
  if(obj4.env<obj.env_min){
    init.env_min<-init4
    obj.env_min<-obj4.env
  }

  return(list(init=as.matrix(init_min), obj=obj_min, init.env=init.env_min, obj.env=obj.env_min,
              init1=init1, init2=init2, init3=init3, init4=init4,
              idx=idx, idx1=idx1, idx2=idx2, idx3=idx3, idx4=idx4,
              obj1=obj1, obj2=obj2, obj3=obj3, obj4=obj4,
              obj1.env=obj1.env, obj2.env=obj2.env, obj3.env=obj3.env, obj4.env=obj4.env))
}

#' Revised version of \code{\link{init.ienv}}
#'
#' Use \code{M+c*Qu} instead of \code{M+Qu}.
#'
#' The four candidates are:
#' \itemize{
#'   \item Standardized U by eigenvectors of \code{M},
#'   \item Normal U by eigenvectors of \code{M},
#'   \item Standardized U by eigenvectors of \code{M+c*Q_U},
#'   \item Normal U by eigenvectors of \code{M+c*Q_U}.
#' }
#'
#' @param c Scale added to \eqn{Q_U}.
#' @return
#' \item{init}{\eqn{\Gamma_0}.}
#' \item{M_cQU}{\eqn{M+c*Q_U}.}
init.ienv2<-function(M, U, d, p, c, Gam0.t=NULL, PRINT=F){

  r<-dim(M)[1]

  #---Use eigenvectors from M---
  if(PRINT){
    cat("## My New Method ##\n\n")
    cat("--Eigenvectors from M--\n")
  }

  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)
  msqrM<-sweep(eig.M$vectors, MARGIN=2, 1/sqrt(eig.M$values), "*")%*%t(eig.M$vectors)

  #---First one: M^{-\frac12} U M^{-\frac12}---

  eig.final1<-eigen((msqrM%*%U%*%msqrM), symmetric=T, only.values=F)
  Pro_small1<-tcrossprod(eig.final1$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of M^{-\frac12} U M^{-\frac12}.

  # We need to find the r-d eigenvectors of M that is not orthogonal to the r-p smallest eigenvectors of M^{-\frac12} U M^{-\frac12} or U
  # So we use "angle" between each eigenvector of M and the space spanned by eig.final$vectors[,(p+1):r].
  # The "angle" is computed as the module of the projected eigenvector.
  ang1<-sqrt(colSums((Pro_small1%*%eig.M$vectors)^2))
  idx1<-sort(ang1, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init1<-as.matrix(eig.M$vectors[,idx1])

  e1<-eigen(t(init1)%*%M%*%init1, symmetric=T, only.values=T)
  e2<-eigen(t(init1)%*%invM%*%init1, only.values=T)
  e3<-eigen(crossprod(init1, (U%*%init1))%*%chol2inv(chol(crossprod(init1, (M%*%init1)))), symmetric=F, only.values=T)
  obj_added<-log(prod(eig.M$values))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("Standardized U.\n")
    cat(frbnorm(init1, Gam0.t), "|", angle(init1, Gam0.t), "|", obj1, "\n")
  }

  init_min<-init1
  obj_min<-obj1


  #---Second one: U---

  eig.U<-eigen(U, symmetric=T, only.values=F)
  QU<-tcrossprod(eig.U$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of U.
  ang2<-sqrt(colSums((QU%*%eig.M$vectors)^2))
  idx2<-sort(ang2, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init2<-as.matrix(eig.M$vectors[,idx2])

  e1<-eigen(t(init2)%*%M%*%init2, symmetric=T, only.values=T)
  e2<-eigen(t(init2)%*%invM%*%init2, only.values=T)
  e3<-eigen(crossprod(init2, (U%*%init2))%*%chol2inv(chol(crossprod(init2, (M%*%init2)))), symmetric=F, only.values=T)
  obj2<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("Normal U.\n")
    cat(frbnorm(init2, Gam0.t), "|", angle(init2, Gam0.t), "|", obj2, "\n")
  }

  if(obj2<obj_min){
    init_min<-init2
    obj_min<-obj2
  }

  ##---Use eigenvectors from M+c*Q_U---

  if(PRINT){
    cat("--Eigenvectors from M+c*Q_U--\n")
  }

  M_cQU<-M+c*QU
  eig.M_QU<-eigen(M_cQU)
  msqrM_QU<-sweep(eig.M_QU$vectors, MARGIN=2, 1/sqrt(eig.M_QU$values), "*")%*%t(eig.M_QU$vectors)

  #---Third one: {M+Q_U}^{-\frac12} U {M+Q_U}^{-\frac12}---

  eig.final3<-eigen((msqrM_QU%*%U%*%msqrM_QU), symmetric=T, only.values=F)
  Pro_small3<-tcrossprod(eig.final3$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of M^{-\frac12} U M^{-\frac12}.

  ang3<-sqrt(colSums((Pro_small3%*%eig.M_QU$vectors)^2))
  idx3<-sort(ang3, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init3<-as.matrix(eig.M_QU$vectors[,idx3])

  e1<-eigen(t(init3)%*%M%*%init3, symmetric=T, only.values=T)
  e2<-eigen(t(init3)%*%invM%*%init3, only.values=T)
  e3<-eigen(crossprod(init3, (U%*%init3))%*%chol2inv(chol(crossprod(init3, (M%*%init3)))), symmetric=F, only.values=T)
  obj3<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("Standardized U.\n")
    cat(frbnorm(init3, Gam0.t), "|", angle(init3, Gam0.t), "|", obj3, "\n")
  }

  if(obj3<obj_min){
    init_min<-init3
    obj_min<-obj3
  }

  #---The fourth candidate: U---

  ang4<-sqrt(colSums((QU%*%eig.M_QU$vectors)^2))
  idx4<-sort(ang4, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init4<-as.matrix(eig.M_QU$vectors[,idx4])

  e1<-eigen(t(init4)%*%M%*%init4, symmetric=T, only.values=T)
  e2<-eigen(t(init4)%*%invM%*%init4, only.values=T)
  e3<-eigen(crossprod(init4, (U%*%init4))%*%chol2inv(chol(crossprod(init4, (M%*%init4)))), symmetric=F, only.values=T)
  obj4<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("Normal U.\n")
    cat(frbnorm(init4, Gam0.t), "|", angle(init4, Gam0.t), "|", obj4, "\n")
  }

  if(obj4<obj_min){
    init_min<-init4
    obj_min<-obj4
  }


  return(list(init=init_min, obj=obj_min, init1=init1, init2=init2, init3=init3,
              init4=init4, obj1=obj1, obj2=obj2, obj3=obj3, obj4=obj4,
              M_cQU=M_cQU))
}

#' OBSOLETE
#'
#' Added a new candidate, selected by the new idea: compress \eqn{U} to make it
#' reduce \eqn{M}.
#'
#' Invocation: \code{\link{stv3}}. And to estimate \eqn{A_0}, only the normal
#' \eqn{U} in \code{stv3} is used.
#'
#' Does not work when \eqn{d=p}.
init.ienv3<-function(M, U, d, p, Gam0.t=NULL, PRINT=F, c){

  if(PRINT){
    cat("## SV 4 INNER ENVE ##\n\n")
  }

  r<-dim(M)[1]
  eig.M<-eigen(M)
  invM<-sweep(eig.M$vectors, MARGIN=2, 1/eig.M$values, "*")%*%t(eig.M$vectors)
  msqrM<-sweep(eig.M$vectors, MARGIN=2, 1/sqrt(eig.M$values), "*")%*%t(eig.M$vectors)

#---New Method---

  W1<-eigen(U)$vectors[,1:p]

  M1<-crossprod(W1, M%*%W1)
  eig.M1<-eigen(M1)
  invM1<-sweep(eig.M1$vectors, MARGIN=2, 1/eig.M1$values, "*")%*%t(eig.M1$vectors)

  M2<-crossprod(W1, invM%*%W1)
  eig.M2<-eigen(M2)
  invM2<-sweep(eig.M2$vectors, MARGIN=2, 1/eig.M2$values, "*")%*%t(eig.M2$vectors)

  if(!is.null(Gam0.t)){
    Gam1.t<-qr.Q(qr(Gam0.t), complete=T)[,(r-d+1):r]
    A_real<-t(W1)%*%Gam1.t
    A0_real<-qr.Q(qr(A_real), complete=T)[,(d+1):p]
  } else{
    A0_real<-NULL
  }

  # CORE: Invocation, using Cook's 4-sv-method.
  if(PRINT){
    cat("----\n")
  }

  A0<-stv3(invM1, M2-invM1, p-d, USE=c(1,3), Gam.t=A0_real, PRINT=PRINT)$init

  if(PRINT){
    cat("----\n")
  }

  A<-qr.Q(qr(A0), complete=T)[,(p-d+1):p]
  sv1_c<-W1%*%A  # Our new sv candidate
  sv1<-qr.Q(qr(sv1_c), complete=T)[,(d+1):r]

  e1<-eigen(t(sv1)%*%M%*%sv1, symmetric=T, only.values=T)
  e2<-eigen(t(sv1)%*%invM%*%sv1, only.values=T)
  e3<-eigen(crossprod(sv1, (U%*%sv1))%*%chol2inv(chol(crossprod(sv1, (M%*%sv1)))), symmetric=F, only.values=T)
  obj_added<-log(prod(eig.M$values))
  obj_sv1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("--New Idea--\n")
    cat(frbnorm(sv1, Gam0.t), "|", angle(sv1, Gam0.t), "|", obj_sv1, "\n")
  }

  init_min<-sv1
  obj_min<-obj_sv1
#---Typical Method---

  #---Use eigenvectors from M---
  if(PRINT){
    cat("--Eigenvectors from M--\n")
  }

  #---First one: M^{-\frac12} U M^{-\frac12}---

  eig.final1<-eigen((msqrM%*%U%*%msqrM), symmetric=T, only.values=F)
  Pro_small1<-tcrossprod(eig.final1$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of M^{-\frac12} U M^{-\frac12}.

  # We need to find the r-d eigenvectors of M that is not orthogonal to the r-p smallest eigenvectors of M^{-\frac12} U M^{-\frac12} or U
  # So we use "angle" between each eigenvector of M and the space spanned by eig.final$vectors[,(p+1):r].
  # The "angle" is computed as the module of the projected eigenvector.
  ang1<-sqrt(colSums((Pro_small1%*%eig.M$vectors)^2))
  idx1<-sort(ang1, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init1<-as.matrix(eig.M$vectors[,idx1])

  e1<-eigen(t(init1)%*%M%*%init1, symmetric=T, only.values=T)
  e2<-eigen(t(init1)%*%invM%*%init1, only.values=T)
  e3<-eigen(crossprod(init1, (U%*%init1))%*%chol2inv(chol(crossprod(init1, (M%*%init1)))), symmetric=F, only.values=T)
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("Standardized U.\n")
    cat(frbnorm(init1, Gam0.t), "|", angle(init1, Gam0.t), "|", obj1, "\n")
  }

  if(obj1<obj_min){
    init_min<-init1
    obj_min<-obj1
  }

  #---Second one: U---

  eig.U<-eigen(U, symmetric=T, only.values=F)
  QU<-tcrossprod(eig.U$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of U.
  ang2<-sqrt(colSums((QU%*%eig.M$vectors)^2))
  idx2<-sort(ang2, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init2<-as.matrix(eig.M$vectors[,idx2])

  e1<-eigen(t(init2)%*%M%*%init2, symmetric=T, only.values=T)
  e2<-eigen(t(init2)%*%invM%*%init2, only.values=T)
  e3<-eigen(crossprod(init2, (U%*%init2))%*%chol2inv(chol(crossprod(init2, (M%*%init2)))), symmetric=F, only.values=T)
  obj2<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("Normal U.\n")
    cat(frbnorm(init2, Gam0.t), "|", angle(init2, Gam0.t), "|", obj2, "\n")
  }

  if(obj2<obj_min){
    init_min<-init2
    obj_min<-obj2
  }

  ##---Use eigenvectors from M+c*Q_U---

  if(PRINT){
    cat("--Eigenvectors from M+c*Q_U--\n")
  }

  eig.M_QU<-eigen(M+c*QU)
  msqrM_QU<-sweep(eig.M_QU$vectors, MARGIN=2, 1/sqrt(eig.M_QU$values), "*")%*%t(eig.M_QU$vectors)

  #---Third one: {M+Q_U}^{-\frac12} U {M+Q_U}^{-\frac12}---

  eig.final3<-eigen((msqrM_QU%*%U%*%msqrM_QU), symmetric=T, only.values=F)
  Pro_small3<-tcrossprod(eig.final3$vectors[,(p+1):r], NULL)  # Projection matrix of the smallest r-p eigenvectors of M^{-\frac12} U M^{-\frac12}.

  ang3<-sqrt(colSums((Pro_small3%*%eig.M_QU$vectors)^2))
  idx3<-sort(ang3, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init3<-as.matrix(eig.M_QU$vectors[,idx3])

  e1<-eigen(t(init3)%*%M%*%init3, symmetric=T, only.values=T)
  e2<-eigen(t(init3)%*%invM%*%init3, only.values=T)
  e3<-eigen(crossprod(init3, (U%*%init3))%*%chol2inv(chol(crossprod(init3, (M%*%init3)))), symmetric=F, only.values=T)
  obj3<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("Standardized U.\n")
    cat(frbnorm(init3, Gam0.t), "|", angle(init3, Gam0.t), "|", obj3, "\n")
  }

  if(obj3<obj_min){
    init_min<-init3
    obj_min<-obj3
  }

  #---The fourth candidate: U---

  ang4<-sqrt(colSums((QU%*%eig.M_QU$vectors)^2))
  idx4<-sort(ang4, decreasing=T, index.return=T)$ix[1:(r-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init4<-as.matrix(eig.M_QU$vectors[,idx4])

  e1<-eigen(t(init4)%*%M%*%init4, symmetric=T, only.values=T)
  e2<-eigen(t(init4)%*%invM%*%init4, only.values=T)
  e3<-eigen(crossprod(init4, (U%*%init4))%*%chol2inv(chol(crossprod(init4, (M%*%init4)))), symmetric=F, only.values=T)
  obj4<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(p-d+1):(r-d)])))

  if(PRINT){
    cat("Normal U.\n")
    cat(frbnorm(init4, Gam0.t), "|", angle(init4, Gam0.t), "|", obj4, "\n")
  }

  if(obj4<obj_min){
    init_min<-init4
    obj_min<-obj4
  }


  return(list(init=init_min, obj=obj_min, init1=init1, init2=init2, init3=init3,
              init4=init4, obj1=obj1, obj2=obj2, obj3=obj3, obj4=obj4,
              eig.M_cQU=eig.M_QU))
}

#' Starting Values for Inner Predictor Envelopes.
#'
#' Only two types of starting values: both using eigenvectors from SX.
#' Mainly copied from \code{\link{init.ienv2}}.
#'
#' The two candidates are:
#' \itemize{
#'   \item Standardized U by eigenvectors of \code{M},
#'   \item Normal U by eigenvectors of \code{M}.
#' }
#'
#' @param SYvX \eqn{S_{Y|X}}.
#' @param c Temporarily obsolete.
#' @param tol Used in \code{\link{Pr}}.
#' @return
#' \item{init}{\eqn{\Gamma_1}.}
sv.IPE<-function(SX, SXY, SYvX, d, r, c, Gam1.t=NULL, PRINT=F, tol=.Machine$double.eps^0.75){

  #---TEST---
  #  d=5; r=7; p=15; eig.SiX=1.5^(1:p); eig.SiYvX=rep(100, r); n=400; SAM=T
  #PRINT=T
  #dat<-data_gen.ixenv(SAM=SAM, d=d, p=p, r=r, eig.SiX=eig.SiX, eig.SiYvX=eig.SiYvX, n=n)

  #if(!SAM){
  #  SX.out<-dat$SiX
  #  SYvX.out<-dat$SiYvX
  #  SXY.out<-dat$SiX%*%dat$beta

  #} else {
  #  X<-dat$X
  #  Y<-dat$Y
  #  SY.out<-stats::cov(Y)*(n-1)/n
  #  SXY.out<-stats::cov(X, Y)*(n-1)/n
  #  SYX.out<-t(SXY.out)
  #  SX.out<-stats::cov(X)*(n-1)/n
  #  SYvX.out<-SY.out-SYX.out%*%chol2inv(chol(SX.out))%*%SXY.out
  #  SYvX.out2<-stats::cov(dat$epsilon)*(n-1)/n
  #}

  #  SX<-SX.out
  #  SXY<-SXY.out
  #  SYvX<-SYvX.out
  #  Gam0.t<-dat$Gamma0
  #--- ---

  p<-dim(SX)[1]

  #---Use eigenvectors from M---
  if(PRINT){
    cat("## SV of IPE ##\n\n")
    cat("--Eigenvectors from SX--\n")
  }

  eig.SX<-eigen(SX, symmetric=T)
  invSX<-sweep(eig.SX$vectors, MARGIN=2, 1/eig.SX$values, "*")%*%t(eig.SX$vectors)
  sqrSX<-sweep(eig.SX$vectors, MARGIN=2, sqrt(eig.SX$values), "*")%*%t(eig.SX$vectors)
  msqrSX<-sweep(eig.SX$vectors, MARGIN=2, 1/sqrt(eig.SX$values), "*")%*%t(eig.SX$vectors)

  eig.SYvX<-eigen(SYvX, symmetric=T)
  invSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, 1/eig.SYvX$values, "*")%*%t(eig.SYvX$vectors)
  msqrSYvX<-sweep(eig.SYvX$vectors, MARGIN=2, 1/sqrt(eig.SYvX$values), "*")%*%t(eig.SYvX$vectors)

  #---First one: M^{-\frac12} U M^{-\frac12}---

  eig.final1<-eigen((msqrSX%*%SXY%*%invSYvX%*%t(SXY)%*%msqrSX), symmetric=T, only.values=F)
  Pro_small1<-tcrossprod(eig.final1$vectors[,(r+1):p], NULL)  # Projection matrix of the smallest r-p eigenvectors of M^{-\frac12} U M^{-\frac12}.

  ang1<-sqrt(colSums((Pro_small1%*%eig.SX$vectors)^2))
  idx1<-sort(ang1, decreasing=T, index.return=T)$ix[1:(p-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init1_com<-as.matrix(eig.SX$vectors[,idx1])
  init1<-qr.Q(qr(init1_com), complete=T)[,(p-d+1):p,drop=F]

  e1<-eigen(t(init1)%*%SX%*%init1, symmetric=T, only.values=T)
  e2<-eigen(t(init1)%*%invSX%*%init1, symmetric=T, only.values=T)

  QSxG<-diag(p)-Pr(sqrSX%*%init1, tol=tol)
  e3<-eigen((msqrSYvX%*%t(SXY)%*%msqrSX%*%QSxG%*%msqrSX%*%SXY%*%msqrSYvX), symmetric=T, only.values=T)

  obj_added<-log(prod(eig.SX$values))
  obj1<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(r-d+1):r])))

  if(PRINT){
    cat("Direct.\n")
    cat(frbnorm(init1, Gam1.t), "|", angle(init1, Gam1.t), "|", obj1, "\n")
  }

  init_min<-init1
  obj_min<-obj1


  #---Second one: SXY%*%invSYvX%*%t(SXY)---

  eig.mid<-eigen(SXY%*%invSYvX%*%t(SXY), symmetric=T, only.values=F)
  Qmid<-tcrossprod(eig.mid$vectors[,(r+1):p], NULL)  # Projection matrix of the smallest r-p eigenvectors of U.
  ang2<-sqrt(colSums((Qmid%*%eig.SX$vectors)^2))
  idx2<-sort(ang2, decreasing=T, index.return=T)$ix[1:(p-d)]  # Indexes of the eigenvectors of M that is closest to the subspace of Pro_small
  init2_com<-as.matrix(eig.SX$vectors[,idx2])
  init2<-qr.Q(qr(init2_com), complete=T)[,(p-d+1):p,drop=F]

  e1<-eigen(t(init1)%*%SX%*%init1, symmetric=T, only.values=T)
  e2<-eigen(t(init1)%*%invSX%*%init1, symmetric=T, only.values=T)

  QSxG<-diag(p)-Pr(sqrSX%*%init2, tol=tol)
  e3<-eigen((msqrSYvX%*%t(SXY)%*%msqrSX%*%QSxG%*%msqrSX%*%SXY%*%msqrSYvX), symmetric=T, only.values=T)

  obj_added<-log(prod(eig.SX$values))
  obj2<-obj_added+sum(log(e1$values))+sum(log(e2$values))+sum(log(1+abs(e3$values[(r-d+1):r])))

  if(PRINT){
    cat("msqrSX removed.\n")
    cat(frbnorm(init2, Gam1.t), "|", angle(init2, Gam1.t), "|", obj2, "\n")
  }

  if(obj2<obj_min){
    init_min<-init2
    obj_min<-obj2
  }

  return(list(init=init_min, obj=obj_min, init1=init1, init2=init2, obj1=obj1, obj2=obj2))
}


#' Get starting values for RD Envelope.
#'
#' Get the starting values for Gamma1 and Gamma12.
#'
#' Invocation: \code{\link{stv}}, \code{\link[Renvlp]{ginv}}.
#'
#' \eqn{Q_\beta} is used when estimating Gamma1. And \code{tol=.Machine$double.eps^0.8} is used in \code{ginv} during this computation.
#'
#' @param u1,u2,r u1+u2=r-1 is allowed!
#' @return Gamma1, Gamma2, Gamma0. Now Gamma2 is strictly inside Gamma12, while two versions of Gamma1 are provided.
starval<-function(M, U, beta, invMU=NULL, u1, u2, r){

  Pbeta<-beta%*%Renvlp::ginv(crossprod(beta), tol=.Machine$double.eps^0.8)%*%t(beta)

  resu1<-stv(M, diag(r)-Pbeta, r-u1, r)  # CORE INVOCATION
  Gam0<-resu1$init
  Gam1<-as.matrix(qr.Q(qr(Gam0), complete=T)[,(r-u1+1):r])

  resu2<-stv(M, U, u1+u2, r)
  Gam12<-resu2$init

  Phi<-qr.solve(Gam12, Gam1)  # I BET IT NEEDS IMPROVEMENT HERE!!!
  Ph<-qr.Q(qr(Phi), complete=T)
  Gam1.12<-Gam12%*%Ph[,1:u1]  # A Gamma1 which is strictly inside Gamma12
  Gam2<-Gam12%*%Ph[,(u1+1):(u1+u2)]


  # For the privilege of use outside udfenvMU().
  if(is.null(invMU)){
    eig.MU<-eigen(M+U)
    invMU<-sweep(eig.MU$vectors, MARGIN=2, 1/eig.MU$values, "*")%*%t(eig.MU$vectors)
  }

  #  e1<-eigen(t(Gam2)%*%M%*%Gam2)
  #  e2<-eigen(t(Gam2)%*%invMU%*%Gam2)
  #  obj2<-sum(log(e1$values))+sum(log(e2$values))  # Value of the minimized objective function of \Gamma_2

  return(list(Gamma1=Gam1, Gamma1.12=Gam1.12, Gamma2=Gam2, Gamma12=Gam12))
}


#' Revised version of \code{\link{starval}}.
#'
#' The only difference is that this function invocates \code{\link{stv2}} instead of \code{\link{stv}}, with corresponding changes in inputs.
#' Things printed can be seen in \code{\link{stv2}}.
#'
#' Gam1 from Pu is used.
starval2<-function(M, U, beta, invMU=NULL, u1, u2, r, Gam1.t, Gam12.t){

  Pbeta<-beta%*%Renvlp::ginv(crossprod(beta), tol=.Machine$double.eps^0.8)%*%t(beta)
  Pu<-U%*%Renvlp::ginv(crossprod(U), tol=.Machine$double.eps^0.8)%*%t(U)
  cat("Gamma1:\nPbeta\n")
  Gam0.t<-as.matrix(qr.Q(qr(Gam1.t), complete=T)[,(u1+1):r])
  resu1.beta<-stv2(M, diag(r)-Pbeta, r-u1, r, Gam0.t)  # CORE INVOCATION
  cat("--------\nPu\n")
  resu1.U<-stv2(M, diag(r)-Pu, r-u1, r, Gam0.t)
  Gam0.beta<-resu1.beta$init
  Gam0.U<-resu1.U$init
  Gam1<-as.matrix(qr.Q(qr(Gam0.U), complete=T)[,(r-u1+1):r])

  cat("-----------------------\nGamma12:\n")
  resu2<-stv2(M, U, u1+u2, r, Gam12.t)
  Gam12<-resu2$init

  Phi<-qr.solve(Gam12, Gam1)  # I BET IT NEEDS IMPROVEMENT HERE!!!
  Ph<-qr.Q(qr(Phi), complete=T)
  Gam1.12<-Gam12%*%Ph[,1:u1]  # A Gamma1 which is strictly inside Gamma12
  Gam2<-Gam12%*%Ph[,(u1+1):(u1+u2)]


  # For the privilege of use outside udfenvMU().
  if(is.null(invMU)){
    eig.MU<-eigen(M+U)
    invMU<-sweep(eig.MU$vectors, MARGIN=2, 1/eig.MU$values, "*")%*%t(eig.MU$vectors)
  }

  return(list(Gamma1=Gam1, Gamma1.12=Gam1.12, Gamma2=Gam2, Gamma12=Gam12))
}

#' Using GEE to estimate inner envelope subspace.
#'
#' Written by Linquan Ma.
#'
#' @param u is the dimension of the inner envelope subspace. u!=0 is required.
#' @examples
#' \dontrun{
#' #---Usage---
#' out<-get_GEE(u, dat$X, dat$Y, 50)
#' out1<-out[-(1:(u*(p-u)))]
#' out2<-out[1:(u*(p-u))]
#' Gamma0<-qr.Q(qr(rbind(diag(r-p), matrix(out1, ncol=r-p))))
#' Gamma1<-qr.Q(qr(Gamma0), complete=T)[, (r-p+1):r]
#' B<-qr.Q(qr(rbind(diag(u), matrix(out2, ncol=u))))
#' Gam1<-Gamma1%*%B  # This is the Estimator of Gamma1
#' }
get_GEE<-function(u, X, Y, num=200){

  p<-ncol(X)
  r<-ncol(Y)
  d<-r-p
  q<-(r-u)*u+(p-u)*(r-p)  # The number of the true parameters needed to specify Gamma_1 and B.

  value<-numeric(num)
  out<-matrix(nrow=num, ncol=q)

  f_GEE<-function(t){

    #    t1<-t[1:(u*(p-u))]  # BUG revised.
    #    t2<-t[-(1:(u*(p-u)))]  # BUG revised.
    Gamma0<-qr.Q(qr(rbind(diag(d), matrix(t[-(1:(u*(p-u)))], ncol=d))))  # Gamma0 is actually Gamma_0 B_0

    #    Gamma1<-qr.Q(qr(Gamma0), complete=T)[, (d+1):r]  # Gamma1 is of r*p
    g_Y<-Y%*%Gamma0
    #    B <- qr.Q(qr(rbind(diag(u), matrix(t1, ncol = u))))
    #    h_Y<-Y%*%Gamma1%*%qr.Q(qr(rbind(diag(u), matrix(t1, ncol=u))))

    return(sum(abs(c(cor(g_Y, X), cor(g_Y, (Y%*%qr.Q(qr(Gamma0), complete=T)[, (d+1):r]%*%qr.Q(qr(rbind(diag(u), matrix(t[1:(u*(p-u))], ncol=u))))))))))  # Sum of the absolute value of each element in the two covariance matrices
  }


  for(i in 1:num){
    start<-runif(q, -5, 5)  # Randomly choose starting values
    tmp<-optim(start, f_GEE, gr=NULL)  # BUG revised.
    value[i]<-tmp$value
    out[i,]<-tmp$par
    #    cat('-')  # For testing.
  }
  #  index<-which(value==min(value))
  index<-which.min(value)

  return(out[index,])
}
