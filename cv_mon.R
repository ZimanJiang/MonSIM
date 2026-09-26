# cross validation for monotone
## prediction function
pred_PLSIM_Est<-function(est,xnew,znew=c()){
  # est: the output object from PLSIM_Est_mon()
  # xnew: new design matrix for nonparametric component
  # znew: new design matrix for parametric component 
  uold = est$x%*%est$alpha
  unew = xnew%*%est$alpha
  n_new = nrow(xnew)
  eta_new = c()
  for(i in 1:n_new){
    w = CalKernel.epan(unew[i]-uold,est$h)
    if(sum(w) == 0)
    {
      idx = which.min(abs(unew[i]-uold))
      if(length(idx)>1)idx==idx[1]
      eta_new[i] = est$eta[idx]
    }
    else
    {
      w = w/sum(w)
      eta_new[i] = sum(est$eta*w);
    }
  }
  if(length(znew)==0) ynew = eta_new
  else ynew = eta_new+znew%*%est$beta
  return(ynew)
}

#cross-validation to select the bandwidth with smallest mse
cv_bw_PLSIM_mon <- function(y,x,z=c(),alpha0,beta0=0,tol=1e-7,max.iter=10000,hlist=c(0.01,0.1,0.5),n.fold = 5,pen.lam=0,sval=1e-4){
  n = length(y)
  id = rep(1:n.fold,ceiling(n/n.fold))
  id = id[1:n]
  id = sample(id,n,replace = F)
  
  mse = matrix(nrow = length(hlist),ncol=n.fold)
  est = list()
  for (i in 1:length(hlist)){
    for(j in 1:n.fold){
      j_id = which(id==j)
      y_train = y[-j_id]
      x_train = x[-j_id,]
      y_test = y[j_id]
      x_test = x[j_id,]
      if(length(z)!=0){
        z_train = z[-j_id,]
        z_test = z[j_id,]
      }
      else{
        z_train = c()
        z_test = c()
      }
      est = PLSIM_Est_mon(y_train,x_train,z_train,hlist[i],alpha0,beta0=0,tol,max.iter,pen.lam = pen.lam,sval=sval)
      y_hat = pred_PLSIM_Est(est,x_test,z_test)
      mse[i,j] = mean(norm(y_test-y_hat,"2"))
    }
  }
  mmse = rowMeans(mse)
  best.id = which.min(mmse)
  best.bw = hlist[best.id]
  mmse.df = data.frame(bw = hlist,mse = mmse)
  return(list(best.bw = best.bw,mse = mmse.df))
}


best_lambda <- function(y,x,z=c(),h,alpha0,beta0=0,tol=1e-7,max.iter=10000,sval=1e-6,lam.list = c(0.01,0.1)){
  BIC.list = rep(Inf,length(lam.list))
  for (i in 1:length(lam.list)){
    res = PLSIM_Est_mon(y,x,z,h,alpha0,beta0,tol,max.iter,sval,pen.lam = lam.list[i])
    BIC.list[i] = res$BIC
    if(min(BIC.list)==BIC.list[i]) temp = res
  }
  best.lam = lam.list[which.min(BIC.list)]
  return(list(best.lambda = best.lam,lambda.BIC.list = data.frame(lambda=lam.list,BIC=BIC.list), best.res = temp))
}

