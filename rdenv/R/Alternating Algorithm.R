#' Revised version of \code{\link{alt_alg}}.
#'
#' Additionally print the value of the objective function after each round of optimization.
alt_alg.p<-function(u1, u2, r, M, invMU, ca, ce, cb1, GEidx.g1, GEidx.g12, GEidx.cb1, obj1){

  maxiter=100
  ftol=1e-04  # Change???
  i<-1
  while(i<maxiter){

    cat("ca\n")
    #------Lu_ca------

    obj3<-obj1
    # These are the addends in the original objective function (before partition)
    C1<-crossprod(ca[GEidx.g12[(u1+u2+1):r], , drop=F])+diag(u1+u2) # C_A^T C_A
    C2<-crossprod(ca, (M%*%ca))  # C_A^T M^{-1} C_A
    C3<-crossprod(ca, (invMU%*%ca))  # C_A^T (M+U) C_A
    C4<-crossprod(cb1, (C1%*%cb1)) # C_E^T C_E=C_{B_1}^T C_A^T C_A C_{B_1}
    C5<-crossprod(cb1, (C2%*%cb1))  # C_{B_1}^T C_A^T M C_A C_{B_1}

    maxiter<-100
    ftol<-0.001
    i<-1
    while(i<maxiter){

      for(j in GEidx.g12[(u1+u2+1):r]){

        a<-as.matrix(ca[j,])  # The "last row" a^T

        t2<-crossprod(ca[-j,],as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
        t3<-crossprod(ca[-j,],as.matrix(invMU[-j,j]))/invMU[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
        W1<-C1-tcrossprod(a, a)  # C_{A_1}^T C_{A_1}
        W2<-C2-tcrossprod((a+t2), (a+t2))*M[j,j]  # W_2
        W3<-C3-tcrossprod((a+t3), (a+t3))*invMU[j,j]

        t5<-crossprod(ca[-j,],as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}

        W4<-C4-tcrossprod(crossprod(cb1, a))  # C_{B_1}^T C_{A_1}^T C_{A_1} C_{B_1}
        W5<-C5-tcrossprod(crossprod(cb1,(a+t5)))*M[j,j]  # W_1

        invW1<-chol2inv(chol(W1))
        invW2<-chol2inv(chol(W2))
        invW3<-chol2inv(chol(W3))
        invW4<-chol2inv(chol(W4))
        invW5<-chol2inv(chol(W5))

        cdinvW4<-cb1%*%tcrossprod(invW4, cb1)  # The middle term of the second addend in L
        cdinvW5<-cb1%*%tcrossprod(invW5, cb1)


        fobj<-function(x){  # Mind the vectors here!!!
          tmp2<-x+t2
          tmp3<-x+t3
          tmp5<-x+t5
          #        A1<-log(1+x%*%(invW1%*%x))
          #        A2.unlog<-1+M[j,j]*crossprod(tmp2, (invW2%*%tmp2))
          #        A2<-log(A2.unlog)
          #        A3<-log(1+invMU[j,j]*crossprod(tmp3, (invW3%*%tmp3)))
          #        A4<-log(1+x%*%(cdinvW4%*%x))
          #        A5<-log(1+M[j,j]*crossprod(tmp5, (cdinvW5%*%tmp5)))

          tmp61<-invW2%*%(W1+tcrossprod(as.matrix(x)))%*%cb1
          mid62<-(W2-M[j,j]*tcrossprod(tmp2)/as.numeric(1+M[j,j]*crossprod(tmp2, (invW2%*%tmp2))))
          #        A6<-sum(log(eigen(crossprod(tmp61, (mid62%*%tmp61)))$values))

          -2*log(1+x%*%(invW1%*%x))+log(1+M[j,j]*crossprod(tmp2, (invW2%*%tmp2)))+log(1+invMU[j,j]*crossprod(tmp3, (invW3%*%tmp3)))-2*log(1+x%*%(cdinvW4%*%x))+log(1+M[j,j]*crossprod(tmp5, (cdinvW5%*%tmp5)))+sum(log(eigen(crossprod(tmp61, (mid62%*%tmp61)))$values))
        }

        res<-stats::optim(ca[j,], fobj, method="BFGS")
        ca[j,]<-res$par
        a<-as.matrix(ca[j,])

        C1<-W1+tcrossprod(a, a)  # Refresh the matrices
        C2<-W2+tcrossprod((a+t2), (a+t2))*M[j,j]
        C3<-W3+tcrossprod((a+t3), (a+t3))*invMU[j,j]
        C4<-W4+tcrossprod(crossprod(cb1, a), crossprod(cb1, a))
        C5<-W5+tcrossprod(crossprod(cb1,a+t5), crossprod(cb1,a+t5))*M[j,j]
      }

      Gam12<-qr.Q(qr(ca))
      ce<-ca%*%cb1
      Gam1<-qr.Q(qr(ce))
      Phi<-qr.solve(Gam12, Gam1)  # SLOW HERE?
      Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
      Gam2<-Gam12%*%Phi0   # A rather coarse way of choosing Gamma2... NEEDS IMPROVEMENT!!!

      e1<-eigen(t(Gam1)%*%M%*%Gam1)
      e2<-eigen(t(Gam2)%*%M%*%Gam2)
      e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
      obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))
      cat("obj5: ", obj5, "\n")

      if(abs(obj3-obj5)<ftol*abs(obj3)){
        obj3<-obj5
        break
      }
      else{
        obj3<-obj5
        i<-i+1
        cat(i-1, "\n")
      }

    }

    cat("cb1\n")
    #------Lu_cb1------

    maxiter<-100
    ftol<-0.001

    M1<-crossprod(ca[GEidx.g12[(u1+u2+1):r], , drop=F])+diag(u1+u2)  # C_A^T C_A
    CM1C<-crossprod(cb1, (M1%*%cb1)) # The first addend
    M2<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
    CM2C<-crossprod(cb1, (M2%*%cb1))  # Second

    eig.caMca<-eigen(M2)
    invcaMca<-sweep(eig.caMca$vectors, MARGIN=2, 1/eig.caMca$values, "*")%*%t(eig.caMca$vectors)
    M3<-crossprod(M1, (invcaMca%*%M1))
    CM3C<-crossprod(cb1, (M3%*%cb1))  # Third


    i<-1
    while(i<maxiter){

      for(j in GEidx.cb1[(u1+1):(u1+u2)]){

        a<-as.matrix(cb1[j,])  # The "last row" a^T
        t1<-crossprod(cb1[-j,],as.matrix(M1[-j,j]))/M1[j,j]
        t2<-crossprod(cb1[-j,],as.matrix(M2[-j,j]))/M2[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
        t3<-crossprod(cb1[-j,],as.matrix(M3[-j,j]))/M3[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}

        W1<-CM1C-tcrossprod((a+t1), (a+t1))*M1[j,j]
        W2<-CM2C-tcrossprod((a+t2), (a+t2))*M2[j,j]
        W3<-CM3C-tcrossprod((a+t3), (a+t3))*M3[j,j]

        invW1<-chol2inv(chol(W1))
        invW2<-chol2inv(chol(W2))
        invW3<-chol2inv(chol(W3))

        fobj<-function(x){  # Objective function
          tmp1<-x+t1
          tmp2<-x+t2
          tmp3<-x+t3
          -2*log(1+M1[j,j]*crossprod(tmp1,invW1%*%tmp1))+log(1+M2[j,j]*crossprod(tmp2,invW2%*%tmp2))+log(1+M3[j,j]*crossprod(tmp3,invW3%*%tmp3))
        }

        gobj<-function(x){  # Derivative
          tmp1<-x+t1
          tmp2<-x+t2
          tmp3<-x+t3
          T1<-invW1%*%tmp1
          T2<-invW2%*%tmp2
          T3<-invW3%*%tmp3
          -4*T1/as.numeric(1/M1[j,j]+crossprod(tmp1,T1))+2*T2/as.numeric(1/M2[j,j]+crossprod(tmp2,T2))+2*T3/as.numeric(1/M3[j,j]+crossprod(tmp3,T3))
        }

        res<-stats::optim(cb1[j,], fobj, gobj, method="BFGS")
        cb1[j,]<-res$par
        a<-as.matrix(cb1[j,])
        CM1C<-W1+tcrossprod((a+t1), (a+t1))*M1[j,j]  # Refresh the matrices
        CM2C<-W2+tcrossprod((a+t2), (a+t2))*M2[j,j]
        CM3C<-W3+tcrossprod((a+t3), (a+t3))*M3[j,j]
      }

      Gam12<-qr.Q(qr(ca))
      ce<-ca%*%cb1
      Gam1<-qr.Q(qr(ce))
      Phi<-qr.solve(Gam12, Gam1)
      Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
      Gam2<-Gam12%*%Phi0   # A rather coarse way of choosing Gamma2... NEEDS IMPROVEMENT!!!

      e1<-eigen(t(Gam1)%*%M%*%Gam1)
      e2<-eigen(t(Gam2)%*%M%*%Gam2)
      e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
      obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))
      cat("obj5: ", obj5, "\n")  # DELETE AFTER TEST

      if(abs(obj3-obj5)<ftol*abs(obj3)){
        obj3<-obj5
        break
      }
      else{
        obj3<-obj5
        i<-i+1
        cat(i-1, "\n")  # DELETE AFTER TEST
      }
    }

    # The outer loop
    if(abs(obj1-obj3)<ftol){
      break
    }
    else{
      obj1<-obj3
      i<-i+1
    }
  }

  return(
    list(ca=ca, ce=ce, cb1=cb1, Gamma1=Gam1, Gamma2=Gam2, Gamma12=Gam12, obj5=obj3)
  )
}

#' Alternating Algorithm.
#'
#' @param u1,u2 u1+u2<r.
alt_alg<-function(u1, u2, r, M, invMU, ca, ce, cb1, GEidx.g1, GEidx.g12, GEidx.cb1, obj1){

  maxiter=100
  ftol=1e-04  # Change???
  i<-1
  while(i<maxiter){

    #------Lu_ca------

    obj3<-obj1
    # These are the addends in the original objective function (before partition)
    C1<-crossprod(ca[GEidx.g12[(u1+u2+1):r], , drop=F])+diag(u1+u2) # C_A^T C_A
    C2<-crossprod(ca, (M%*%ca))  # C_A^T M^{-1} C_A
    C3<-crossprod(ca, (invMU%*%ca))  # C_A^T (M+U) C_A
    C4<-crossprod(cb1, (C1%*%cb1)) # C_E^T C_E=C_{B_1}^T C_A^T C_A C_{B_1}
    C5<-crossprod(cb1, (C2%*%cb1))  # C_{B_1}^T C_A^T M C_A C_{B_1}

    maxiter<-100
    ftol<-0.001
    i<-1
    while(i<maxiter){

      for(j in GEidx.g12[(u1+u2+1):r]){

        a<-as.matrix(ca[j,])  # The "last row" a^T

        t2<-crossprod(ca[-j,],as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
        t3<-crossprod(ca[-j,],as.matrix(invMU[-j,j]))/invMU[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
        W1<-C1-tcrossprod(a, a)  # C_{A_1}^T C_{A_1}
        W2<-C2-tcrossprod((a+t2), (a+t2))*M[j,j]  # W_2
        W3<-C3-tcrossprod((a+t3), (a+t3))*invMU[j,j]

        t5<-crossprod(ca[-j,],as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}

        W4<-C4-tcrossprod(crossprod(cb1, a))  # C_{B_1}^T C_{A_1}^T C_{A_1} C_{B_1}
        W5<-C5-tcrossprod(crossprod(cb1,(a+t5)))*M[j,j]  # W_1

        invW1<-chol2inv(chol(W1))
        invW2<-chol2inv(chol(W2))
        invW3<-chol2inv(chol(W3))
        invW4<-chol2inv(chol(W4))
        invW5<-chol2inv(chol(W5))

        cdinvW4<-cb1%*%tcrossprod(invW4, cb1)  # The middle term of the second addend in L
        cdinvW5<-cb1%*%tcrossprod(invW5, cb1)


        fobj<-function(x){  # Mind the vectors here!!!
          tmp2<-x+t2
          tmp3<-x+t3
          tmp5<-x+t5
          #        A1<-log(1+x%*%(invW1%*%x))
          #        A2.unlog<-1+M[j,j]*crossprod(tmp2, (invW2%*%tmp2))
          #        A2<-log(A2.unlog)
          #        A3<-log(1+invMU[j,j]*crossprod(tmp3, (invW3%*%tmp3)))
          #        A4<-log(1+x%*%(cdinvW4%*%x))
          #        A5<-log(1+M[j,j]*crossprod(tmp5, (cdinvW5%*%tmp5)))

          tmp61<-invW2%*%(W1+tcrossprod(as.matrix(x)))%*%cb1
          mid62<-(W2-M[j,j]*tcrossprod(tmp2)/as.numeric(1+M[j,j]*crossprod(tmp2, (invW2%*%tmp2))))
          #        A6<-sum(log(eigen(crossprod(tmp61, (mid62%*%tmp61)))$values))

          -2*log(1+x%*%(invW1%*%x))+log(1+M[j,j]*crossprod(tmp2, (invW2%*%tmp2)))+log(1+invMU[j,j]*crossprod(tmp3, (invW3%*%tmp3)))-2*log(1+x%*%(cdinvW4%*%x))+log(1+M[j,j]*crossprod(tmp5, (cdinvW5%*%tmp5)))+sum(log(eigen(crossprod(tmp61, (mid62%*%tmp61)))$values))
        }

        res<-stats::optim(ca[j,], fobj, method="BFGS")
        ca[j,]<-res$par
        a<-as.matrix(ca[j,])

        C1<-W1+tcrossprod(a, a)  # Refresh the matrices
        C2<-W2+tcrossprod((a+t2), (a+t2))*M[j,j]
        C3<-W3+tcrossprod((a+t3), (a+t3))*invMU[j,j]
        C4<-W4+tcrossprod(crossprod(cb1, a), crossprod(cb1, a))
        C5<-W5+tcrossprod(crossprod(cb1,a+t5), crossprod(cb1,a+t5))*M[j,j]
      }

      Gam12<-qr.Q(qr(ca))
      ce<-ca%*%cb1
      Gam1<-qr.Q(qr(ce))
      Phi<-qr.solve(Gam12, Gam1)  # SLOW HERE?
      Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
      Gam2<-Gam12%*%Phi0   # A rather coarse way of choosing Gamma2... NEEDS IMPROVEMENT!!!

      e1<-eigen(t(Gam1)%*%M%*%Gam1)
      e2<-eigen(t(Gam2)%*%M%*%Gam2)
      e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
      obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))

      if(abs(obj3-obj5)<ftol*abs(obj3)){
        obj3<-obj5
        break
      }
      else{
        obj3<-obj5
        i<-i+1
      }
    }

    #------Lu_cb1------

    maxiter<-100
    ftol<-0.001

    M1<-crossprod(ca[GEidx.g12[(u1+u2+1):r], , drop=F])+diag(u1+u2)  # C_A^T C_A
    CM1C<-crossprod(cb1, (M1%*%cb1)) # The first addend
    M2<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
    CM2C<-crossprod(cb1, (M2%*%cb1))  # Second

    eig.caMca<-eigen(M2)
    invcaMca<-sweep(eig.caMca$vectors, MARGIN=2, 1/eig.caMca$values, "*")%*%t(eig.caMca$vectors)
    M3<-crossprod(M1, (invcaMca%*%M1))
    CM3C<-crossprod(cb1, (M3%*%cb1))  # Third


    i<-1
    while(i<maxiter){

      for(j in GEidx.cb1[(u1+1):(u1+u2)]){

        a<-as.matrix(cb1[j,])  # The "last row" a^T
        t1<-crossprod(cb1[-j,],as.matrix(M1[-j,j]))/M1[j,j]
        t2<-crossprod(cb1[-j,],as.matrix(M2[-j,j]))/M2[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
        t3<-crossprod(cb1[-j,],as.matrix(M3[-j,j]))/M3[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}

        W1<-CM1C-tcrossprod((a+t1), (a+t1))*M1[j,j]
        W2<-CM2C-tcrossprod((a+t2), (a+t2))*M2[j,j]
        W3<-CM3C-tcrossprod((a+t3), (a+t3))*M3[j,j]

        invW1<-chol2inv(chol(W1))
        invW2<-chol2inv(chol(W2))
        invW3<-chol2inv(chol(W3))

        fobj<-function(x){  # Objective function
          tmp1<-x+t1
          tmp2<-x+t2
          tmp3<-x+t3
          -2*log(1+M1[j,j]*crossprod(tmp1,invW1%*%tmp1))+log(1+M2[j,j]*crossprod(tmp2,invW2%*%tmp2))+log(1+M3[j,j]*crossprod(tmp3,invW3%*%tmp3))
        }

        gobj<-function(x){  # Derivative
          tmp1<-x+t1
          tmp2<-x+t2
          tmp3<-x+t3
          T1<-invW1%*%tmp1
          T2<-invW2%*%tmp2
          T3<-invW3%*%tmp3
          -4*T1/as.numeric(1/M1[j,j]+crossprod(tmp1,T1))+2*T2/as.numeric(1/M2[j,j]+crossprod(tmp2,T2))+2*T3/as.numeric(1/M3[j,j]+crossprod(tmp3,T3))
        }

        res<-stats::optim(cb1[j,], fobj, gobj, method="BFGS")
        cb1[j,]<-res$par
        a<-as.matrix(cb1[j,])
        CM1C<-W1+tcrossprod((a+t1), (a+t1))*M1[j,j]  # Refresh the matrices
        CM2C<-W2+tcrossprod((a+t2), (a+t2))*M2[j,j]
        CM3C<-W3+tcrossprod((a+t3), (a+t3))*M3[j,j]
      }

      Gam12<-qr.Q(qr(ca))
      ce<-ca%*%cb1
      Gam1<-qr.Q(qr(ce))
      Phi<-qr.solve(Gam12, Gam1)
      Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
      Gam2<-Gam12%*%Phi0   # A rather coarse way of choosing Gamma2... NEEDS IMPROVEMENT!!!

      e1<-eigen(t(Gam1)%*%M%*%Gam1)
      e2<-eigen(t(Gam2)%*%M%*%Gam2)
      e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
      obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))

      if(abs(obj3-obj5)<ftol*abs(obj3)){
        obj3<-obj5
        break
      }
      else{
        obj3<-obj5
        i<-i+1
      }
    }

    # The outer loop
    if(abs(obj1-obj3)<ftol){
      break
    }
    else{
      obj1<-obj3
      i<-i+1
    }
  }

  return(
    list(ca=ca, ce=ce, cb1=cb1, Gamma1=Gam1, Gamma2=Gam2, Gamma12=Gam12, obj5=obj3)
  )
}


#' Compute value of the final objective function over C_A with C_{B_1} fixed.
#'
#' For other information, see \code{\link{Lu_cb1}}. These two functions are paired,
#' having the same inputs and outputs.
Lu_ca<-function(u1, u2, r, M, invMU, ca, ce, cb1, GEidx.g1, GEidx.g12, GEidx.cb1, obj1){

  # These are the addends in the original objective function (before partition)
  C1<-crossprod(ca[GEidx.g12[(u1+u2+1):r], , drop=F])+diag(u1+u2) # C_A^T C_A
  C2<-crossprod(ca, (M%*%ca))  # C_A^T M^{-1} C_A
  C3<-crossprod(ca, (invMU%*%ca))  # C_A^T (M+U) C_A
  C4<-crossprod(cb1, (C1%*%cb1)) # C_E^T C_E=C_{B_1}^T C_A^T C_A C_{B_1}
  C5<-crossprod(cb1, (C2%*%cb1))  # C_{B_1}^T C_A^T M C_A C_{B_1}


  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    for(j in GEidx.g12[(u1+u2+1):r]){

      a<-as.matrix(ca[j,])  # The "last row" a^T

      t2<-crossprod(ca[-j,],as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      t3<-crossprod(ca[-j,],as.matrix(invMU[-j,j]))/invMU[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
      W1<-C1-tcrossprod(a)  # C_{A_1}^T C_{A_1}
      W2<-C2-tcrossprod(a+t2)*M[j,j]  # W_2
      W3<-C3-tcrossprod(a+t3)*invMU[j,j]

      t5<-crossprod(ca[-j,],as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}

      W4<-C4-tcrossprod(crossprod(cb1, a))  # C_{B_1}^T C_{A_1}^T C_{A_1} C_{B_1}
      W5<-C5-tcrossprod(crossprod(cb1,a+t5))*M[j,j]  # W_1

      invW1<-chol2inv(chol(W1))
      invW2<-chol2inv(chol(W2))
      invW3<-chol2inv(chol(W3))
      invW4<-chol2inv(chol(W4))
      invW5<-chol2inv(chol(W5))

      cdinvW4<-cb1%*%tcrossprod(invW4, cb1)  # The middle term of the second addend in L
      cdinvW5<-cb1%*%tcrossprod(invW5, cb1)


      fobj<-function(x){  # Mind the vectors here!!!
        tmp2<-x+t2
        tmp3<-x+t3
        tmp5<-x+t5
        A1<-log(1+x%*%(invW1%*%x))
        A2.unlog<-1+M[j,j]*crossprod(tmp2, (invW2%*%tmp2))
        A2<-log(A2.unlog)
        A3<-log(1+invMU[j,j]*crossprod(tmp3, (invW3%*%tmp3)))
        A4<-log(1+x%*%(cdinvW4%*%x))
        A5<-log(1+M[j,j]*crossprod(tmp5, (cdinvW5%*%tmp5)))

        tmp61<-invW2%*%(W1+tcrossprod(as.matrix(x)))%*%cb1
        mid62<-(W2-M[j,j]*tcrossprod(tmp2)/as.numeric(A2.unlog))
        A6.unlog<-crossprod(tmp61, (mid62%*%tmp61))
        eig.A6<-eigen(A6.unlog)
        A6<-sum(log(eig.A6$values))

        return(-2*A1+A2+A3-2*A4+A5+A6)
      }

      res<-stats::optim(ca[j,], fobj, method="BFGS")
      ca[j,]<-res$par
      a<-as.matrix(ca[j,])

      C1<-W1+tcrossprod(a)  # Refresh the matrices
      C2<-W2+tcrossprod(a+t2)*M[j,j]
      C3<-W3+tcrossprod(a+t3)*invMU[j,j]
      C4<-W4+tcrossprod(crossprod(cb1, a))
      C5<-W5+tcrossprod(crossprod(cb1,a+t5))*M[j,j]
    }

    Gam12<-qr.Q(qr(ca))
    ce<-ca%*%cb1
    Gam1<-qr.Q(qr(ce))
    Phi<-qr.solve(Gam12, Gam1)
    Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
    Gam2<-Gam12%*%Phi0   # A rather coarse way of choosing Gamma2... NEEDS IMPROVEMENT!!!

    e1<-eigen(t(Gam1)%*%M%*%Gam1)
    e2<-eigen(t(Gam2)%*%M%*%Gam2)
    e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
    obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }
  return(list(ca=ca, ce=ce, cb1=cb1, Gamma1=Gam1, Gamma2=Gam2, Gamma12=Gam12, obj5=obj5))
}


#' Compute value of the final objective function over C_{B_1} with C_A fixed.
#'
#' This function is paired with \code{\link{Lu_ca}}.
#'
#' @param u1,u2 u1+u2!=r. While u1+u2<r-1 is also allowed.
#' @return C_A, C_{B_1}, C_E, Gamma1/2/12 and the value of the minimized objective function.
Lu_cb1<-function(u1, u2, r, M, invMU, ca, ce, cb1, GEidx.g1, GEidx.g12, GEidx.cb1, obj1){

  maxiter<-100
  ftol<-0.001

  M1<-crossprod(ca[GEidx.g12[(u1+u2+1):r], , drop=F])+diag(u1+u2)  # C_A^T C_A
  CM1C<-crossprod(cb1, (M1%*%cb1)) # The first addend
  M2<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  CM2C<-crossprod(cb1, (M2%*%cb1))  # Second

  eig.caMca<-eigen(M2)
  invcaMca<-sweep(eig.caMca$vectors, MARGIN=2, 1/eig.caMca$values, "*")%*%t(eig.caMca$vectors)
  M3<-crossprod(M1, (invcaMca%*%M1))
  CM3C<-crossprod(cb1, (M3%*%cb1))  # Third


  i<-1
  while(i<maxiter){

    for(j in GEidx.cb1[(u1+1):(u1+u2)]){

      a<-as.matrix(cb1[j,])  # The "last row" a^T
      t1<-crossprod(cb1[-j,],as.matrix(M1[-j,j]))/M1[j,j]
      t2<-crossprod(cb1[-j,],as.matrix(M2[-j,j]))/M2[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
      t3<-crossprod(cb1[-j,],as.matrix(M3[-j,j]))/M3[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}

      W1<-CM1C-tcrossprod(a+t1)*M1[j,j]
      W2<-CM2C-tcrossprod(a+t2)*M2[j,j]
      W3<-CM3C-tcrossprod(a+t3)*M3[j,j]

      invW1<-chol2inv(chol(W1))
      invW2<-chol2inv(chol(W2))
      invW3<-chol2inv(chol(W3))

      fobj<-function(x){  # Objective function
        tmp1<-x+t1
        tmp2<-x+t2
        tmp3<-x+t3
        A1<-log(1+M1[j,j]*crossprod(tmp1,invW1%*%tmp1))
        A2<-log(1+M2[j,j]*crossprod(tmp2,invW2%*%tmp2))
        A3<-log(1+M3[j,j]*crossprod(tmp3,invW3%*%tmp3))
        return(-2*A1+A2+A3)
      }

      gobj<-function(x){  # Derivative
        tmp1<-x+t1
        tmp2<-x+t2
        tmp3<-x+t3
        T1<-invW1%*%tmp1
        T2<-invW2%*%tmp2
        T3<-invW3%*%tmp3
        -4*T1/as.numeric(1/M1[j,j]+crossprod(tmp1,T1))+2*T2/as.numeric(1/M2[j,j]+crossprod(tmp2,T2))+2*T3/as.numeric(1/M3[j,j]+crossprod(tmp3,T3))
      }

      res<-stats::optim(cb1[j,], fobj, gobj, method="BFGS")
      cb1[j,]<-res$par
      a<-as.matrix(cb1[j,])
      CM1C<-W1+tcrossprod(a+t1)*M1[j,j]  # Refresh the matrices
      CM2C<-W2+tcrossprod(a+t2)*M2[j,j]
      CM3C<-W3+tcrossprod(a+t3)*M3[j,j]

    }

    Gam12<-qr.Q(qr(ca))
    ce<-ca%*%cb1
    Gam1<-qr.Q(qr(ce))
    Phi<-qr.solve(Gam12, Gam1)
    Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
    Gam2<-Gam12%*%Phi0   # A rather coarse way of choosing Gamma2... NEEDS IMPROVEMENT!!!

    e1<-eigen(t(Gam1)%*%M%*%Gam1)
    e2<-eigen(t(Gam2)%*%M%*%Gam2)
    e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
    obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }

  }
  return(list(ca=ca, ce=ce, cb1=cb1, Gamma1=Gam1, Gamma2=Gam2, Gamma12=Gam12, obj5=obj5))
}


#' Compute value of the final objective function over C_A with C_{B_1} fixed.
#'
#' Paired with \code{\link{Lu_cb1.rm1}}.
Lu_ca.rm1<-function(u1, u2, r, M, invMU, ca, ce, cb1, GEidx.g1, GEidx.g12, GEidx.cb1, obj1){

  # These are the addends in the original objective function (before partition)
  C1<-crossprod(ca[GEidx.g12[(u1+u2+1):r],,drop=F])+diag(u1+u2) # C_A^T C_A
  C2<-crossprod(ca, (M%*%ca))  # C_A^T M^{-1} C_A
  C3<-crossprod(ca, (invMU%*%ca))  # C_A^T (M+U) C_A
  C4<-crossprod(cb1, (C1%*%cb1)) # C_E^T C_E=C_{B_1}^T C_A^T C_A C_{B_1}
  C5<-crossprod(cb1, (C2%*%cb1))  # C_{B_1}^T C_A^T M C_A C_{B_1}


  maxiter<-100
  ftol<-0.001
  i<-1
  while(i<maxiter){

    j<-GEidx.g12[r]

    a<-as.matrix(ca[j,])  # The "last row" a^T

    t2<-crossprod(ca[-j,],as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
    t3<-crossprod(ca[-j,],as.matrix(invMU[-j,j]))/invMU[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}
    W1<-C1-tcrossprod(a)  # C_{A_1}^T C_{A_1}
    W2<-C2-tcrossprod(a+t2)*M[j,j]  # W_2
    W3<-C3-tcrossprod(a+t3)*invMU[j,j]

    t5<-crossprod(ca[-j,],as.matrix(M[-j,j]))/M[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}

    W4<-C4-tcrossprod(crossprod(cb1, a))  # C_{B_1}^T C_{A_1}^T C_{A_1} C_{B_1}
    W5<-C5-tcrossprod(crossprod(cb1,a+t5))*M[j,j]  # W_1

#    invW1<-chol2inv(chol(W1))
    invW2<-chol2inv(chol(W2))
    invW3<-chol2inv(chol(W3))
    invW4<-chol2inv(chol(W4))
    invW5<-chol2inv(chol(W5))

    cdinvW4<-cb1%*%tcrossprod(invW4, cb1)  # The middle term of the second addend in L
    cdinvW5<-cb1%*%tcrossprod(invW5, cb1)


    fobj<-function(x){  # Mind the vectors here!!!
      tmp2<-x+t2
      tmp3<-x+t3
      tmp5<-x+t5
      A1<-log(1+sum(x^2))
      A2.unlog<-1+M[j,j]*crossprod(tmp2, (invW2%*%tmp2))
      A2<-log(A2.unlog)
      A3<-log(1+invMU[j,j]*crossprod(tmp3, (invW3%*%tmp3)))
      A4<-log(1+x%*%(cdinvW4%*%x))
      A5<-log(1+M[j,j]*crossprod(tmp5, (cdinvW5%*%tmp5)))

      tmp61<-invW2%*%(W1+tcrossprod(as.matrix(x)))%*%cb1
      mid62<-(W2-M[j,j]*tcrossprod(tmp2)/as.numeric(A2.unlog))
      A6.unlog<-crossprod(tmp61, (mid62%*%tmp61))
      eig.A6<-eigen(A6.unlog)
      A6<-sum(log(eig.A6$values))

      return(-2*A1+A2+A3-2*A4+A5+A6)
    }

    res<-stats::optim(ca[j,], fobj, method="BFGS")
    ca[j,]<-res$par
    a<-as.matrix(ca[j,])

    C1<-W1+tcrossprod(a)  # Refresh the matrices
    C2<-W2+tcrossprod(a+t2)*M[j,j]
    C3<-W3+tcrossprod(a+t3)*invMU[j,j]
    C4<-W4+tcrossprod(crossprod(cb1, a))
    C5<-W5+tcrossprod(crossprod(cb1,a+t5))*M[j,j]

    Gam12<-qr.Q(qr(ca))
    ce<-ca%*%cb1
    Gam1<-qr.Q(qr(ce))
    Phi<-qr.solve(Gam12, Gam1)
    Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
    Gam2<-Gam12%*%Phi0   # A rather coarse way of choosing Gamma2... NEEDS IMPROVEMENT!!!

    e1<-eigen(t(Gam1)%*%M%*%Gam1)
    e2<-eigen(t(Gam2)%*%M%*%Gam2)
    e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
    obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }
  }

  return(list(ca=ca, ce=ce, cb1=cb1, Gamma1=Gam1, Gamma2=Gam2, Gamma12=Gam12, obj5=obj5))
}


#' Compute value of the final objective function over C_{B_1} with C_A fixed.
#'
#' This function is paired with \code{\link{Lu_ca.rm1}}.
#'
#' @param u1,u2 u1+u2=r-1.
#' @return C_A, C_{B_1}, C_E, Gamma1/2/12 and the value of the minimized objective function.
Lu_cb1.rm1<-function(u1, u2, r, M, invMU, ca, ce, cb1, GEidx.g1, GEidx.g12, GEidx.cb1, obj1){

  maxiter<-100
  ftol<-0.001

  M1<-crossprod(ca[GEidx.g12[(u1+u2+1):r],,drop=F])+diag(u1+u2)  # C_A^T C_A
  CM1C<-crossprod(cb1, (M1%*%cb1)) # The first addend
  M2<-crossprod(ca, (M%*%ca))  # C_A^T M C_A
  CM2C<-crossprod(cb1, (M2%*%cb1))  # Second

  eig.caMca<-eigen(M2)
  invcaMca<-sweep(eig.caMca$vectors, MARGIN=2, 1/eig.caMca$values, "*")%*%t(eig.caMca$vectors)
  M3<-crossprod(M1, (invcaMca%*%M1))
  CM3C<-crossprod(cb1, (M3%*%cb1))  # Third


  i<-1
  while(i<maxiter){

    j<-GEidx.cb1[(u1+u2)]

    a<-as.matrix(cb1[j,])  # The "last row" a^T
    t1<-crossprod(cb1[-j,],as.matrix(M1[-j,j]))/M1[j,j]
    t2<-crossprod(cb1[-j,],as.matrix(M2[-j,j]))/M2[j,j]  # M_{22}^{-1} C{A_1}^T M_{12}
    t3<-crossprod(cb1[-j,],as.matrix(M3[-j,j]))/M3[j,j]  # V_{22}^{-1} C{A_1}^T V_{12}

    W1<-CM1C-tcrossprod(a+t1)*M1[j,j]
    W2<-CM2C-tcrossprod(a+t2)*M2[j,j]
    W3<-CM3C-tcrossprod(a+t3)*M3[j,j]

    invW1<-chol2inv(chol(W1))
    invW2<-chol2inv(chol(W2))
    invW3<-chol2inv(chol(W3))

    fobj<-function(x){  # Objective function
      tmp1<-x+t1
      tmp2<-x+t2
      tmp3<-x+t3
      A1<-log(1+M1[j,j]*crossprod(tmp1,invW1%*%tmp1))
      A2<-log(1+M2[j,j]*crossprod(tmp2,invW2%*%tmp2))
      A3<-log(1+M3[j,j]*crossprod(tmp3,invW3%*%tmp3))
      return(-2*A1+A2+A3)
    }

    gobj<-function(x){  # Derivative
      tmp1<-x+t1
      tmp2<-x+t2
      tmp3<-x+t3
      T1<-invW1%*%tmp1
      T2<-invW2%*%tmp2
      T3<-invW3%*%tmp3
      -4*T1/as.numeric(1/M1[j,j]+crossprod(tmp1,T1))+2*T2/as.numeric(1/M2[j,j]+crossprod(tmp2,T2))+2*T3/as.numeric(1/M3[j,j]+crossprod(tmp3,T3))
    }

    res<-stats::optim(cb1[j,], fobj, gobj, method="BFGS")
    cb1[j,]<-res$par
    a<-as.matrix(cb1[j,])
    CM1C<-W1+tcrossprod(a+t1)*M1[j,j]  # Refresh the matrices
    CM2C<-W2+tcrossprod(a+t2)*M2[j,j]
    CM3C<-W3+tcrossprod(a+t3)*M3[j,j]


    Gam12<-qr.Q(qr(ca))
    ce<-ca%*%cb1
    Gam1<-qr.Q(qr(ce))
    Phi<-qr.solve(Gam12, Gam1)
    Phi0<-as.matrix(qr.Q(qr(Phi), complete=T)[,(u1+1):(u1+u2)])
    Gam2<-Gam12%*%Phi0   # A rather coarse way of choosing Gamma2... NEEDS IMPROVEMENT!!!

    e1<-eigen(t(Gam1)%*%M%*%Gam1)
    e2<-eigen(t(Gam2)%*%M%*%Gam2)
    e3<-eigen(t(Gam12)%*%invMU%*%Gam12)
    obj5<-sum(log(e1$values))+sum(log(e2$values))+sum(log(e3$values))

    if(abs(obj1-obj5)<ftol*abs(obj1)){
      break
    }
    else{
      obj1<-obj5
      i<-i+1
    }
  }

  return(list(ca=ca, ce=ce, cb1=cb1, Gamma1=Gam1, Gamma2=Gam2, Gamma12=Gam12, obj5=obj5))
}
