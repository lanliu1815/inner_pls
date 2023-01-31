#' Save Parameters of a IPE Model as txt Files
#'
#' Saves \code{Gamma1, Gamma0, beta, eta1, eta2, C, eig.SiX, eig.SiYvX} as txt
#' files. Notice only the eigenvalues--\code{eig.SiX, eig.SiYvX}--are saved.
#'
#' Notice that \code{setwd} function is used! Be cautious of possible
#' side effects!
#'
#' @param dat The R object containing the parameters.
#' @param work_dir0 The working directory to write the data on.
#' @param SC Whether it is a special case?
#'
#' @return A folder at \code{work_dir0}, the name of which is
#' \code{"d", d, "r", r, "p", p, "n", n, " ", time_now}.
data_save<-function(dat, work_dir0="Y:/Research/Data Generated and Results/IPE", SC=F){

  WD<-getwd()

  time_now<-gsub(":", "-", Sys.time())

  p<-nrow(dat$Gamma1)
  r<-length(dat$eig.SiYvX)
  d<-ncol(dat$Gamma1)
  # Create a folder whose name composes the values used for each parameter. Notice that "dia" is followed by the first element used in dia.
  setwd(work_dir0)
  if(SC){
    fdnm_exp<-paste0("d", d, "r", r, "p", p, " ", time_now, " ", "SC")
  } else {
    fdnm_exp<-paste0("d", d, "r", r, "p", p, " ", time_now)
  }
  dir.create(fdnm_exp)
  work_dir1<-paste0(work_dir0, "/", fdnm_exp)
  setwd(work_dir1)

  # The parameters to export(to pass to Matlab).
  par=c("Gamma1", "Gamma0", "beta", "eta1", "eta2", "C", "eig.SiX", "eig.SiYvX")

  for (j in 1:length(par)) {
    write.table(dat[[par[j]]], file=paste(par[j], ".txt", sep=""), row.names=F, col.names=F)
  }

  setwd(WD)

}

#' Read Parameters of a IPE Model from a Folder
#'
#' Reads \code{"Gamma1", "Gamma0", "C", "beta", "eta1", "eta2", "eig.SiX", "eig.SiYvX"}.
#' Also restores \code{SiX} and \code{SiYvX}
#'
#' @param fdnm The name of the folder containing the txt files.
#' @param work_dir0 The working directory where the folder is.
#'
#' @return
#' \item{dat}{A list containing the parameters.}
data_READ<-function(fdnm, work_dir0="Y:/Research/Data Generated and Results/IPE"){

#  work_dir0<-"Y:/Research/Data Generated and Results/IPE"
#  fdnm_exp<-paste0("d", d, "r", r, "p", p, "n", n, further)
  WD<-getwd()

  work_dir1<-paste0(work_dir0, "/", fdnm)
  setwd(work_dir1)

  #par=c("X", "Y", "Gamma1", "Gamma0", "Gamma0_sv", "B", "beta", "betaOLS", "M")
  par=c("Gamma1", "Gamma0", "C", "beta", "eta1", "eta2", "eig.SiX", "eig.SiYvX")

  dat<-list()
  for (j in 1:length(par)){
    dat[[par[j]]]<-as.matrix(read.table(file=paste0(par[j], ".txt"), check.names=F))
    colnames(dat[[par[j]]])<-NULL
  }
#  names(dat)<-par

  # "Restore" SiX & SiYvX
  dat$Gamma<-cbind(dat$Gamma1, dat$Gamma0)
  p<-nrow(dat$Gamma1)
  r<-length(dat$eig.SiYvX)
  dat$SiX<-sweep(dat$Gamma, MARGIN=2, dat$eig.SiX[1:p], "*")%*%t(dat$Gamma)
  A<-rsomat(r, r)
  dat$SiYvX<-sweep(A, MARGIN=2, dat$eig.SiYvX[1:r], "*")%*%t(A)

  setwd(WD)

  return(dat)
}

#' Run the Simulation of inner PLS, IPE, Predictor Envelopes and SIMPLS
#'
#' OBSOLETE because too slow due to endless "if/else" loop.
#'
#' Four possible estimators: inner PLS, IPE, Predictor Envelopes and SIMPLS.
#'
#' @param dat A list containing the parameters of a inner predictor envelope.
#' @param rep_num Replication times of the simulation.
#' @param estmrs The estimators to use.
#' @param SELECT_DIM Whether to select dimensions of the estimators. The dimension
#' selection is done by \code{\link{d.boot.sim}}.
#' @param err_type The error type to use, as listed now. If more than one, the
#' first will be used.
#' @param coef_normal Coefficient to time the SiYvX matrix.
#' @param upper.uniform Will generate \code{U(0, upper.uniform)} errors for uniform
#' errors.
#' @param B The bootstrap replication used in dimension selection.
#' @param NonGrass Whether to further implement optimization after obtaining the
#' starting values of IPE.
#'
#' @return A list. Objects with \code{.iPLS} have counterparts with \code{.IPE, .predenv, .PLS}.
#' \item{dis}{Distance between estimators of beta and the real one in each replication.}
#' \item{dimen}{Dimension selected for estimators in each replication}
#' \item{mean_dis}
#' \item{mean_dimen}
#' \item{MSE}{The mean squared error of the estimators.}
#' \item{dimen.iPLS}{The bootstrap vector correlations of each dimension in each
#'  replication. Will be \code{NA} if \code{SELECT_DIM=F}.}
#' \item{mean_dimen.iPLS}
#' \item{time_used}
Run_Simulation0<-function(dat, rep_num=5, estmrs=c("iPLS", "IPE", "PredEnv", "SIMPLS"), SELECT_DIM=T, err_type=c("normal", "uniform"), coef_normal=1, upper.uniform, B=50, n=1000, NonGrass=F, ftol.nonG=0.5){

  #---TEST---
#  rep_num=5; NonGrass=F; B=50; n=1000
#  SELECT_DIM=F
#  err_type<-c("normal")
#  estmrs<-c("iPLS", "IPE", "PredEnv", "SIMPLS")
  #--- --- ---

  if(length(err_type)!=1){
    warning("Specify a single err_type!")
    err_type<-err_type[1]
  }
  if(!err_type%in%c("normal", "uniform")){
    stop("Specify the right err_type!")
  }

  esti_nm<-c("iPLS", "IPE", "PredEnv", "SIMPLS")
  # Make sure estmrs contain right names
  if(!all(estmrs%in%esti_nm)){
    stop("Specify the right estimator in estmrs!")
  }


  PRINT=F;

  resu_num<-length(estmrs)+1
  dim_num<-length(estmrs)

  d<-ncol(dat$Gamma1)
  p<-nrow(dat$Gamma1)
  r<-nrow(dat$SiYvX)

  dis<-matrix(NA, nrow=rep_num, ncol=resu_num, dimnames=list(NULL, c("OLS", esti_nm[esti_nm%in%estmrs])))
  dimen<-matrix(NA, nrow=rep_num, ncol=dim_num, dimnames=list(NULL, esti_nm[esti_nm%in%estmrs]))

  dimen.SV<-matrix(NA, nrow=rep_num, ncol=r)
  dimen.iPLS<-matrix(NA, nrow=rep_num, ncol=r)
  dimen.PLS<-matrix(NA, nrow=rep_num, ncol=p-1)
  dimen.predenv<-matrix(NA, nrow=rep_num, ncol=p-1)
  colnames(dimen.SV)<-1:r
  colnames(dimen.iPLS)<-1:r
  colnames(dimen.predenv)<-1:(p-1)
  colnames(dimen.PLS)<-1:(p-1)

  output<-list()

  ptm<-proc.time()
  for (i in 1:rep_num) {

    X<-MASS::mvrnorm(n, rep(0, p), dat$SiX)
    eps<-switch(err_type,
                "normal"=MASS::mvrnorm(n, rep(0, r), dat$SiYvX*coef_normal),
                "uniform"=matrix(runif(n=n*r, max=upper.uniform), nrow=n, ncol=r))
    # Normal
 #   eps<-MASS::mvrnorm(n, rep(0, r), dat$SiYvX); err_type<-"Normal"
    # Student's t
    #  eps<-matrix(rt(n=n*r, df=6), nrow=n, ncol=r); err_type<-"t6"
    # U(0,?)
    #  eps<-matrix(runif(n=n*r, max=upper.uniform), nrow=n, ncol=r); err_type<-"Uniform"
    # Chi Square
    #  eps<-matrix(rchisq(n=n*r, df=4), nrow=n, ncol=r); err_type<-"Chisq"

    Y<-X%*%dat$beta+eps

    # Sample Covariances

    SY<-stats::cov(Y)*(n-1)/n
    SYX<-stats::cov(Y, X)*(n-1)/n
    SXY<-t(SYX)
    SX<-stats::cov(X)*(n-1)/n
    invSX<-chol2inv(chol(SX))
    invSY<-chol2inv(chol(SY))
    betaOLS<-invSX%*%t(SYX)

    dis[i,"OLS"]<-norm(betaOLS-dat$beta, type="F")

    # Dimension Selection
    if(SELECT_DIM){

      dime<-d.boot.sim(X=X, Y=Y, VERBOSE=F, B=B, estmrs=estmrs)
    } else{
      dime<-list()
      dime$d.boot.SV<-dime$d.boot.iPLS<-d
      dime$u.boot.PLS<-dime$u.boot.predenv<-p
      dime$dis.SV<-dime$dis.iPLS<-dime$dis.PLS<-dime$dis.predenv<-NA
    }

    # innerPLS
    if("iPLS"%in%estmrs){

      dimen.iPLS[i,]<-dime$dis.iPLS
      dimen[i,"iPLS"]<-d.iPLS<-dime$d.boot.iPLS

      if(d.iPLS!=0){
        Gam0.pls<-SIMPLS_gamma(M=SX, U=diag(p)-Pr(SXY), u=p-d.iPLS)$Gammahat
        Gam1.pls<-qr.Q(qr(Gam0.pls), complete=T)[,(p-d.iPLS+1):p,drop=F]
      }

      innerpls<-IPE(X=X, Y=Y, d=d.iPLS, initial=Gam1.pls, NonGrass=F)
      dis[i,"iPLS"]<-norm(innerpls$MLE_otherpars$betaIPE-dat$beta, type="F")

    }

    # IPE
    if("IPE"%in%estmrs){

      dimen.SV[i,]<-dime$dis.SV
      dimen[i,"IPE"]<-d.IPE<-dime$d.boot.SV
      mle<-IPE(X=X, Y=Y, d=d.IPE, NonGrass=NonGrass)

      dis[i,"IPE"]<-norm(mle$MLE_otherpars$betaIPE-dat$beta, type="F")


    }

    # Predictor Envelope
    if("PredEnv"%in%estmrs){

      dimen.predenv[i,]<-dime$dis.predenv
      #  resu.predenv<-Renvlp::u.xenv(X=X, Y=Y)
      #  u.predenv<-resu.predenv$u.aic
      dimen[i,"PredEnv"]<-u.predenv<-dime$u.boot.predenv
      enve_pred<-Renvlp::xenv(X=X, Y=Y, u=u.predenv, asy=F)
      dis[i,"PredEnv"]<-norm(enve_pred$beta-dat$beta, type="F")

    }

    # PLS
    if("SIMPLS"%in%estmrs){

      dimen.PLS[i,]<-dime$dis.PLS
      dimen[i,"SIMPLS"]<-u.PLS<-dime$u.boot.PLS
      Gam.pls<-SIMPLS_gamma(M=SX, U=tcrossprod(SXY), u=u.PLS)$Gammahat
      beta.PLS<-Gam.pls%*%chol2inv(chol(t(Gam.pls)%*%SX%*%Gam.pls))%*%t(Gam.pls)%*%SXY
      dis[i,"SIMPLS"]<-norm(beta.PLS-dat$beta, type="F")

    }

    cat("-")
  }

  cat("\n")

  output$time_used<-proc.time()-ptm

  if(SELECT_DIM){

    if("iPLS"%in%estmrs){
      output$dimen.iPLS<-dimen.iPLS
      output$mean_dimen.iPLS<-colMeans(dimen.iPLS)
    }

    if("IPE"%in%estmrs){
      output$dimen.SV<-dimen.SV
      output$mean_dimen.SV<-colMeans(dimen.SV)
    }

    if("PredEnv"%in%estmrs){
      output$dimen.predenv<-dimen.predenv
      output$mean_dimen.predenv<-colMeans(dimen.predenv)
    }

    if("SIMPLS"%in%estmrs){
      output$dimen.PLS<-dimen.PLS
      output$mean_dimen.PLS<-colMeans(dimen.PLS)
    }

    output$dimen<-dimen
    output$mean_dimen<-colMeans(dimen)

  }

  output$dis<-dis
  output$mean_dis<-colMeans(dis)


#  MSE<-matrix(NA, ncol=resu_num, nrow=1, dimnames=list(NULL, c("OLS", esti_nm[esti_nm%in%estmrs])))
#  apply(X=dis^2, MARGIN=2, FUN=sum)/rep_num

#  for (k in 1:resu_num) {
#    MSE[1,k]<-sum(dis[,k]^2)/rep_num
#  }
  output$MSE<-apply(X=dis^2, MARGIN=2, FUN=sum)/rep_num

  # Print a boxplot?
#  dis_boxplot<-dis
#  colnames(dis_boxplot)<-c("OLS", esti_nm[esti_nm%in%estmrs])
#  boxplot(dis_boxplot)

  return(output)

}



#' Run the Simulation of inner PLS, IPE, Predictor Envelopes and SIMPLS
#'
#' Run the Simulation of inner PLS, IPE, Predictor Envelopes and SIMPLS.
#' Coded vectorizingly. Should be quick.
#'
#' Prints "-" each time \code{\link{d.boot.sim}} completes in a replication.
#'
#' @param dat A list containing the parameters of a inner predictor envelope.
#' @param rep_num Replication times of the simulation.
#' @param estmrs The estimators to use. \code{PredEnv} depends on package \code{Renvlp},
#' and \code{LASSO,Ridge} depend on package \code{glmnet}.
#' @param SELECT_DIM Whether to perform dimension selection. If \code{FALSE}, true
#' value is used.
#' @param u.PredEnv,u.SIMPLS Specified dimension to be used for predictor envelopes
#' and SIMPLS. Valid only when \code{SELECT_DIM=FALSE}.
#' @param dim.iPLS_IPE The type of dimension selection method to use for inner PLS and SIMPLS.
#' If set to \code{CV}, both estimators must be specified in \code{estmrs}.
#' If set to \code{BIC}, iPLS is proceeded by bootstrap.
#' @param dim.SIMPLS The type of dimension selection method to use for SIMPLS.
#' If "bootstrap", bootstrap is used for \code{SELECT_DIM=TRUE} and
#' true value otherwise. If \code{estmrs} does not contain \code{"SIMPLS"}, this
#' argument does not matter.
#' @param dim.PredEnv Similar to \code{dim.SIMPLS}.
#' @param err_type The error type to use, as listed now. If more than one, the
#' first will be used.
#' @param coef_normal Coefficient to time the SiYvX matrix.
#' @param upper.uniform Will generate \code{U(0, upper.uniform)} errors for uniform
#' errors.
#' @param B The bootstrap replication used in dimension selection.
#' @param kfold The cross-validation folds used in dimension selection, for SIMPLS.
#' @param NonGrass Whether to further implement optimization after obtaining the
#' starting values of IPE.
#' @param lambda_seq Provided to \code{\link[glmnet]{cv.glmnet}} as \code{lambda=lambda_seq}
#' for LASSO and Ridge.
#'
#' @return A list. Simulation results and arguments used. Important ones are:
#' \item{dis}{Distance between estimators(including OLS) of beta and the real
#' one in each replication.}
#' \item{dimen}{Dimension selected for estimators in each replication}
#' \item{mean_dis}
#' \item{mean_dimen}
#' \item{MSE}{The mean squared error of the estimators.}
#' \item{time_used}
#' \item{dimen_dis.iPLS}{The bootstrap vector correlations of each dimension in each
#'  replication. If CV is used, this does not appear. If dimension selection is
#'  performed, \code{NA} is returned.}
#' \item{mean_dimen_dis.iPLS}
Run_Simulation<-function(dat, rep_num=5, estmrs=c("iPLS", "IPE", "PredEnv", "SIMPLS", "LASSO", "Ridge"), SELECT_DIM=F, u.PredEnv=NULL, u.SIMPLS=NULL, dim.iPLS_IPE=c("bootstrap", "CV", "BIC"), dim.PredEnv=c("bootstrap", "CV", "BIC"), dim.SIMPLS=c("bootstrap", "CV"), err_type=c("Normal", "Uniform"), coef_normal=1, upper.uniform, B=50, kfold=3, n=1000, NonGrass=F, ftol.nonG=0.5, lambda_seq=seq(from=0, to=20, by=0.2)){

  #---TEST---
#    rep_num=5; NonGrass=F; B=50; n=1000; kfold=5; coef_normal=1; upper.uniform=3
#    dim.iPLS_IPE="BIC"; dim.SIMPLS="CV"; dim.PredEnv="BIC"
#    err_type<-c("Uniform")
#    estmrs<-c("PredEnv", "IPE", "LASSO", "Ridge")
#    SELECT_DIM=T
#    lambda_seq=seq(from=0, to=20, by=0.2)
#    u.SIMPLS=8; u.PredEnv=9
  #--- --- ---

  d<-ncol(dat$Gamma1); p<-nrow(dat$Gamma1); r<-nrow(dat$SiYvX)
  if(length(err_type)!=1){
    stop("Specify a single err_type!")
  }
  if(!err_type%in%c("Normal", "Uniform")){
    stop("Specify the right err_type!")
  }
  if(!is.null(u.PredEnv)){
    if(u.PredEnv>p|u.PredEnv<1)
      stop("Wrong u.SIMPLS!")
  }
  if(!is.null(u.SIMPLS)){
    if(u.SIMPLS>p|u.SIMPLS<1)
    stop("Wrong u.SIMPLS!")
  }
  if(length(dim.iPLS_IPE)!=1){
    stop("Specify a single dim.iPLS_IPE!")
  }
  if(!dim.iPLS_IPE%in%c("bootstrap", "CV", "BIC")){
    stop("Specify the right dim.iPLS_IPE!")
  }
  if(dim.iPLS_IPE=="CV"&!all(c("iPLS", "IPE")%in%estmrs)){
    stop("Specify both iPLS and IPE in estmrs if using CV in dim.iPLS_IPE!")
  }
  if(length(dim.SIMPLS)!=1){
    stop("Specify a single dim.SIMPLS!")
  }
  if(!dim.SIMPLS%in%c("bootstrap", "CV")){
    stop("Specify the right dim.SIMPLS!")
  }
  if(length(dim.PredEnv)!=1){
    stop("Specify a single dim.PredEnv!")
  }
  if(!dim.PredEnv%in%c("bootstrap", "CV", "BIC")){
    stop("Specify the right dim.PredEnv!")
  }

  esti_nm<-c("iPLS", "IPE", "PredEnv", "SIMPLS", "LASSO", "Ridge")
  esti_dim_nm<-c("iPLS", "IPE", "PredEnv", "SIMPLS")
  # Make sure estmrs contain right names
  if(!all(estmrs%in%esti_nm)){
    stop("Specify the right estimator in estmrs!")
  }

  dis<-matrix(NA, nrow=rep_num, ncol=length(estmrs)+1, dimnames=list(NULL, c("OLS", esti_nm[esti_nm%in%estmrs])))
  dimen<-matrix(NA, nrow=rep_num, ncol=sum(estmrs%in%esti_dim_nm), dimnames=list(NULL, estmrs[estmrs%in%esti_dim_nm]))

  output<-list(rep_num=rep_num, estmrs=estmrs, SELECT_DIM=SELECT_DIM, err_type=err_type, n=n, B=B, kfold=kfold, dat=dat, dim.PredEnv=dim.PredEnv, dim.iPLS_IPE=dim.iPLS_IPE, dim.SIMPLS=dim.SIMPLS)
  if(err_type=="Normal"){
    output$coef_normal=coef_normal
  }
  if(err_type=="Uniform"){
    output$upper.uniform=upper.uniform
  }
  X.list<-Y.list<-list()

  ptm<-proc.time()

  # Generate X and Y, and perform OLS estimation
  for(i in 1:rep_num){

    X.list[[i]]<-MASS::mvrnorm(n, rep(0, p), dat$SiX)
    eps<-switch(err_type,
                "Normal"=MASS::mvrnorm(n, rep(0, r), dat$SiYvX*coef_normal),
                "Uniform"=matrix(runif(n=n*r, max=upper.uniform)-upper.uniform/2, nrow=n, ncol=r))
    Y.list[[i]]<-X.list[[i]]%*%dat$beta+eps

    betaOLS<-chol2inv(chol(stats::cov(X.list[[i]])*(n-1)/n))%*%stats::cov(X.list[[i]], Y.list[[i]])*(n-1)/n
    dis[i,"OLS"]<-norm(betaOLS-dat$beta, type="F")
  }

  # Dimension Selection
  if(SELECT_DIM){
    estmrs_boot<-estmrs[estmrs%in%c("iPLS", "IPE", "PredEnv", "SIMPLS")]

    if(dim.SIMPLS=="CV"&"SIMPLS"%in%estmrs){
      cat("Dimension Selection by CV for SIMPLS\n")
      dimsn_SIMPLS.list<-mapply(FUN=function(X, Y) u.cv.SIMPLS(X=X, Y=Y, kfold=kfold, trace=-1),
                         X.list, Y.list, SIMPLIFY=F)
      estmrs_boot<-setdiff(estmrs_boot, "SIMPLS")
      cat("\n")
    }
    if(dim.PredEnv=="CV"&"PredEnv"%in%estmrs){
      cat("Dimension Selection by CV for Predictor Envelope\n")
      dimsn_PredEnv.list<-mapply(FUN=function(X, Y) u.cv.PredEnv(X=X, Y=Y, kfold=kfold, trace=-1),
                                X.list, Y.list, SIMPLIFY=F)
      estmrs_boot<-setdiff(estmrs_boot, "PredEnv")
      cat("\n")
    } else if(dim.PredEnv=="BIC"&"PredEnv"%in%estmrs){
      cat("Dimension Selection by BIC for Predictor Envelope\n")
      dimsn_PredEnv.list<-mapply(FUN=function(X, Y) Renvlp::u.xenv(X=X, Y=Y),
                                 X.list, Y.list, SIMPLIFY=F)
      estmrs_boot<-setdiff(estmrs_boot, "PredEnv")
      cat("\n")
    }
    if(dim.iPLS_IPE=="CV"&all(c("iPLS", "IPE")%in%estmrs)){
      cat("Dimension Selection by CV for inner PLS and IPE\n")
      dimsn_iPLS_IPE.list<-mapply(FUN=function(X, Y) u.cv.sim(X=X, Y=Y, kfold=kfold, trace=-1, NonGrass=NonGrass),
                                X.list, Y.list, SIMPLIFY=F)
      estmrs_boot<-setdiff(estmrs_boot, c("iPLS", "IPE"))
      cat("\n")
    }
    if(dim.iPLS_IPE=="BIC" & "IPE"%in%estmrs){
      cat("Dimension Selection by BIC for IPE\n")
      dimsn_IPE.list<-mapply(FUN=function(X, Y) aicbic.IPE(X=X, Y=Y, trace=-1, NonGrass=NonGrass),
                             X.list, Y.list, SIMPLIFY=F)
      estmrs_boot<-setdiff(estmrs_boot, c("IPE"))
      cat("\n")
    }

    cat("Dimension Selection by Bootstrap for", estmrs_boot, "\n")
    dimsn.list<-mapply(FUN=function(X, Y) d.boot.sim(X=X, Y=Y, B=B, estmrs=estmrs_boot, trace=-1),
                       X.list, Y.list, SIMPLIFY=F)
    cat("\n")
  } else{
    cat("True/Specified Dimensions Used\n")
    if(is.null(u.SIMPLS)){
      u.SIMPLS<-p
    } else{
      cat("Dimension of SIMPLS specified as:", u.SIMPLS, "\n")
    }
    if(is.null(u.PredEnv)){
      u.PredEnv<-p
    } else{
      cat("Dimension of Predictor Envelopes specified as:", u.PredEnv, "\n")
    }
    dimsn.list<-lapply(X=1:rep_num, FUN=function(x) list(d.boot.IPE=d, d.boot.iPLS=d, u.boot.PLS=u.SIMPLS, u.boot.PredEnv=u.PredEnv, dis.IPE=NA, dis.iPLS=NA, dis.PredEnv=NA, dis.PLS=NA))
  }

  # Estimation of estimators

  # innerPLS
  if("iPLS"%in%estmrs){
    if(SELECT_DIM&dim.iPLS_IPE=="CV"){
      dis[,"iPLS"]<-sapply(X=1:rep_num, FUN=function(x) norm(innerPLS(X=X.list[[x]], Y=Y.list[[x]], d=dimsn_iPLS_IPE.list[[x]]$d.cv.iPLS)$betaiPLS-dat$beta, type="F"))
      dimen[,"iPLS"]<-sapply(X=1:rep_num, FUN=function(x) dimsn_iPLS_IPE.list[[x]]$d.cv.iPLS)
    } else{
      dis[,"iPLS"]<-sapply(X=1:rep_num, FUN=function(x) norm(innerPLS(X=X.list[[x]], Y=Y.list[[x]], d=dimsn.list[[x]]$d.boot.iPLS)$betaiPLS-dat$beta, type="F"))
      output$dimen_dis.iPLS<-t(sapply(X=1:rep_num, FUN=function(x) dimsn.list[[x]]$dis.iPLS))
      output$mean_dimen_dis.iPLS<-colMeans(output$dimen_dis.iPLS)
      dimen[,"iPLS"]<-sapply(X=1:rep_num, FUN=function(x) dimsn.list[[x]]$d.boot.iPLS)
    }
  }

  # IPE
  if("IPE"%in%estmrs){
    if(SELECT_DIM&dim.iPLS_IPE=="CV"){
      dis[,"IPE"]<-sapply(X=1:rep_num, FUN=function(x) norm(IPE(X=X.list[[x]], Y=Y.list[[x]], d=dimsn_iPLS_IPE.list[[x]]$d.cv.IPE, NonGrass=NonGrass)$betaIPE-dat$beta, type="F"))
      dimen[,"IPE"]<-sapply(X=1:rep_num, FUN=function(x) dimsn_iPLS_IPE.list[[x]]$d.cv.IPE)
    } else if(SELECT_DIM&dim.iPLS_IPE=="BIC"){
      dis[,"IPE"]<-sapply(X=1:rep_num, FUN=function(x) norm(IPE(X=X.list[[x]], Y=Y.list[[x]], d=dimsn_IPE.list[[x]]$d.bic, NonGrass=NonGrass)$betaIPE-dat$beta, type="F"))
      dimen[,"IPE"]<-sapply(X=1:rep_num, FUN=function(x) dimsn_IPE.list[[x]]$d.bic)
    } else{
      dis[,"IPE"]<-sapply(X=1:rep_num, FUN=function(x) norm(IPE(X=X.list[[x]], Y=Y.list[[x]], d=dimsn.list[[x]]$d.boot.IPE, NonGrass=NonGrass)$betaIPE-dat$beta, type="F"))
      output$dimen_dis.IPE<-t(sapply(X=1:rep_num, FUN=function(x) dimsn.list[[x]]$dis.IPE))
      output$mean_dimen_dis.IPE<-colMeans(output$dimen_dis.IPE)
      dimen[,"IPE"]<-sapply(X=1:rep_num, FUN=function(x) dimsn.list[[x]]$d.boot.IPE)
    }
  }

  # Predictor Envelope
  if("PredEnv"%in%estmrs){
    if(SELECT_DIM&dim.PredEnv=="CV"){
      dis[,"PredEnv"]<-sapply(X=1:rep_num, FUN=function(x) norm(Renvlp::xenv(X=X.list[[x]], Y=Y.list[[x]], u=dimsn_PredEnv.list[[x]]$u.cv.PredEnv)$beta-dat$beta, type="F"))
      dimen[,"PredEnv"]<-sapply(X=1:rep_num, FUN=function(x) dimsn_PredEnv.list[[x]]$u.cv.PredEnv)
    } else if(SELECT_DIM&dim.PredEnv=="BIC"){
      dis[,"PredEnv"]<-sapply(X=1:rep_num, FUN=function(x) norm(Renvlp::xenv(X=X.list[[x]], Y=Y.list[[x]], u=dimsn_PredEnv.list[[x]]$u.bic)$beta-dat$beta, type="F"))
      dimen[,"PredEnv"]<-sapply(X=1:rep_num, FUN=function(x) dimsn_PredEnv.list[[x]]$u.bic)
      output$dimsn_PredEnv.list<-dimsn_PredEnv.list
    } else {
      dis[,"PredEnv"]<-sapply(X=1:rep_num, FUN=function(x) norm(Renvlp::xenv(X=X.list[[x]], Y=Y.list[[x]], u=dimsn.list[[x]]$u.boot.PredEnv, asy=F)$beta-dat$beta, type="F"))
      output$dimen_dis.PredEnv<-t(sapply(X=1:rep_num, FUN=function(x) dimsn.list[[x]]$dis.PredEnv))
      output$mean_dimen_dis.PredEnv<-colMeans(output$dimen_dis.PredEnv)
      dimen[,"PredEnv"]<-sapply(X=1:rep_num, FUN=function(x) dimsn.list[[x]]$u.boot.PredEnv)
    }
  }

  # SIMPLS
  if("SIMPLS"%in%estmrs){

    if(SELECT_DIM&dim.SIMPLS=="CV"){
      dis[,"SIMPLS"]<-sapply(X=1:rep_num, FUN=function(x) norm(SIMPLS_env(X=X.list[[x]], Y=Y.list[[x]], u=dimsn_SIMPLS.list[[x]]$u.cv)$beta.PLS-dat$beta, type="F"))
      dimen[,"SIMPLS"]<-sapply(X=1:rep_num, FUN=function(x) dimsn_SIMPLS.list[[x]]$u.cv)
    } else{
      dis[,"SIMPLS"]<-sapply(X=1:rep_num, FUN=function(x) norm(SIMPLS_env(X=X.list[[x]], Y=Y.list[[x]], u=dimsn.list[[x]]$u.boot.PLS)$beta.PLS-dat$beta, type="F"))
      output$dimen_dis.PLS<-t(sapply(X=1:rep_num, FUN=function(x) dimsn.list[[x]]$dis.PLS))
      output$mean_dimen_dis.PLS<-colMeans(output$dimen_dis.PLS)
      dimen[,"SIMPLS"]<-sapply(X=1:rep_num, FUN=function(x) dimsn.list[[x]]$u.boot.PLS)
    }
  }


  # LASSO
  if("LASSO"%in%estmrs){
    cat("LASSO\n")
    LASSO_sim<-function(X, Y){
      beta.lasso<-matrix(NA, nrow=p, ncol=r)
      for (j in 1:r) {
        cv.lasso<-glmnet::cv.glmnet(x=X, y=Y[,j], alpha=1, intercept=F, thresh=1e-14, lambda=lambda_seq)
        beta.lasso[,j]<-as.numeric(glmnet::glmnet(x=X, y=Y[,j], lambda=cv.lasso$lambda.min, thresh=1e-14)$beta)  # Will this be slow?
      }
      cat("-")  # To show the process
      return(beta.lasso)
    }
    dis[,"LASSO"]<-sapply(X=1:rep_num, FUN=function(x) norm(LASSO_sim(X=X.list[[x]], Y=Y.list[[x]])-dat$beta, type="F"))
    cat("\n")
  }

  # Ridge
  if("Ridge"%in%estmrs){
    cat("Ridge\n")
    Ridge_sim<-function(X, Y){
      SXY<-stats::cov(X, Y)*(n-1)/n
      SX<-stats::cov(X)*(n-1)/n
      beta.ridge<-matrix(NA, nrow=p, ncol=r)
      for (k in 1:r) {
        cv.ridge<-glmnet::cv.glmnet(x=X, y=Y[,k], alpha=0, lambda=lambda_seq)  # lambda is supplied
        beta.ridge[,k]<-chol2inv(chol(SX*n+cv.ridge$lambda.min*diag(p)))%*%SXY[,k]*n
      }
      cat("-")  # To show the process
      return(beta.ridge)
    }
    dis[,"Ridge"]<-sapply(X=1:rep_num, FUN=function(x) norm(Ridge_sim(X=X.list[[x]], Y=Y.list[[x]])-dat$beta, type="F"))
    cat("\n")
  }

  output$dimen<-dimen
  output$dis<-dis

  # Summary statistics
  output$time_used<-proc.time()-ptm
  output$mean_dimen<-colMeans(output$dimen)
  output$mean_dis<-colMeans(output$dis)
  output$MSE<-apply(X=output$dis^2, MARGIN=2, FUN=sum)/rep_num

  cat("\n")

  return(output)

}


#' Write Results to Files
#'
#' Creates or overwrites a folder named of the type of distribution used.
#'
#' @param result_simulation The output from \code{\link{Run_Simulation}}.
#' @param fdnm The name of the folder that the data is read from.
#' @param BOXPLOT Draw a boxplot of \code{result_simulation$dis}? Use
#' \code{\link[graphics]{boxplot}}. The arguments are
#' \code{outline=F, ylim=c(0,max(result_Simulation$dis)), ylab="Bias"}.
#' @param estinm The names of the estimators used.
#' @param work_dir0 The name of the father working directory of \code{fdnm}.
#' @param ... Optional arguments to \code{\link[graphics]{boxplot}}. For example,
#' \code{names=c(expression(beta[OLS]), expression(beta[iPLS]), expression(beta[IPE], expression(beta[Predenv])))}.
results_save<-function(result_Simulation, fdnm, BOXPLOT=T, estinm=c("OLS", result_Simulation$estmrs), dat=result_Simulation$dat, work_dir0="Y:/Research/Data Generated and Results/IPE", ylim_boxplot=c(0,max(result_Simulation$dis)), ylab_boxplot="Bias", ...){

  #---TEST---
#  result_Simulation<-resu_2
#  estinm=c("OLS", result_Simulation$estmrs)
#  work_dir0="Y:/Research/Data Generated and Results/IPE"
#  fdnm="d5r6p102021-12-16 19-46-19"
#  BOXPLOT=T
  #--- ---

  WD<-getwd()

  dis_write<-rbind(c(result_Simulation$time_used[3], rep(NA, length(estinm)-1)), colMeans(result_Simulation$dis), rep(NA, length(estinm)), result_Simulation$dis, deparse.level=0)
  colnames(dis_write)<-colnames(result_Simulation$dis)
  dimen_write<-rbind(colMeans(result_Simulation$dimen), rep(NA, length(result_Simulation$dimen[1,])), result_Simulation$dimen, deparse.level=1)

  work_dir1<-paste0(work_dir0, "/", fdnm)
  setwd(work_dir1)

  time_now<-gsub(":", "-", Sys.time())
  tpnm<-paste0(result_Simulation$err_type, " ", time_now)
  dir.create(tpnm)
  work_dir2<-paste0(work_dir1, "/", tpnm)
  setwd(work_dir2)

  write.table(x=round(dis_write, digits=5), file=paste0("0 distance", " ", time_now, ".txt"), row.names=F, col.names=T, sep=" | ")
  write.table(x=round(t(result_Simulation$MSE), digits=5), file=paste0("0 MSE", " ", time_now, ".txt"), row.names=F, col.names=T, sep=" | ")
  write.table(x=round(dimen_write, digits=5), file=paste0("0 dimension", " ", time_now, ".txt"), row.names=F, col.names=T, sep=" | ")

  if(result_Simulation$SELECT_DIM){
    if("iPLS"%in%estinm & result_Simulation$dim.iPLS_IPE=="bootstrap"){
      write.table(x=round(result_Simulation$dimen_dis.iPLS, digits=5), file=paste0("1 inp dimension selection", " ", time_now, ".txt"), row.names=F, col.names=as.character(1:nrow(dat$SiYvX)), sep=" | ")
    }
    if("IPE"%in%estinm & result_Simulation$dim.iPLS_IPE=="bootstrap"){
      write.table(x=round(result_Simulation$dimen_dis.IPE, digits=5), file=paste0("1 SV dimension selection", " ", time_now, ".txt"), row.names=F, col.names=as.character(1:nrow(dat$SiYvX)), sep=" | ")
    }
    if("PredEnv"%in%estinm & result_Simulation$dim.PredEnv=="bootstrap"){
      write.table(x=round(result_Simulation$dimen_dis.PredEnv, digits=5), file=paste0("1 predenv dimension selection", " ", time_now, ".txt"), row.names=F, col.names=as.character(1:nrow(dat$SiX)), sep=" | ")
    }
    if("SIMPLS"%in%estinm & result_Simulation$dim.SIMPLS=="bootstrap"){
      write.table(x=round(result_Simulation$dimen_dis.PLS, digits=5), file=paste0("1 PLS dimension selection", " ", time_now, ".txt"), row.names=F, col.names=as.character(1:nrow(dat$SiX)), sep=" | ")
    }

  } else{
    write.table(x=t(result_Simulation$mean_dimen), file="True Dimensions or Dimensions Specified.txt", row.names=F, col.names=T, sep="|")
  }

  write.table(x=NULL, file=paste0("rep.num=", result_Simulation$rep_num))
  write.table(x=NULL, file=paste0("B=", result_Simulation$B))
  write.table(x=NULL, file=paste0("kfold=", result_Simulation$kfold))
  write.table(x=NULL, file=paste0("n=", result_Simulation$n))


  if(result_Simulation$err_type=="Normal"){
    write.table(x=NULL, file=paste0("coef.normal=", result_Simulation$coef_normal))
  }
  if(result_Simulation$err_type=="Uniform"){
    write.table(x=NULL, file=paste0("upper.uniform=", result_Simulation$upper.uniform))
  }

  if(BOXPLOT){
    pdf(file=paste0("Boxplot_", tpnm, ".pdf"), height=6, width=9)
    boxplot(result_Simulation$dis, outline=F, ylim=ylim_boxplot, ylab=ylab_boxplot,...)
    #names=c(expression(beta[OLS]), expression(beta[iPLS]), expression(beta[IPE], expression(beta[Predenv])))
    dev.off()
  }

  setwd(WD)

}



