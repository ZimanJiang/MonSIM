

CalKernel.epan <- function(u,h){
  kernel = rep(0,length(u))
  w = u/h
  id = which(w<1&w>-1)
  kernel[id] = 0.75*(1-w[id]^2)/h
  return(kernel)
}

## estimate the nonparametric link 20241010
CalEta_mon <- function(y,x,z=c(),alpha,beta,h,sval=1e-6,inc = T){
  eps = 1e-16 #correction term
  lam = x%*%alpha #n*1
  ord = order(lam)
  lam = as.matrix(lam[ord]) #n*1
  #xstar = as.matrix(x[ord,])
  if(length(z)==0){
    zstar = z
  }else{
    zstar = as.matrix(z[ord,])
  }
  
  lam.mat = matrix(rep(lam,nrow(lam)),ncol = nrow(lam),byrow = T) #n*n
  if(length(z)==0)ystar = y[ord]
  else ystar = y[ord]-zstar%*%beta #n*1
  lam.u = lam.mat-t(lam.mat) #n*n
  #calculate Kernel estimated eta
  #K0 = apply(lam.u,1,CalKernel.epan,h) #n*n, weights of every (lambda-u)
  K0 = apply(lam.u,1,dnorm,mean=0,sd=h) #n*n, weights of every (lambda-u)
  K00 = rowSums(K0) # 1*n
  K10 = rowSums(K0*lam.u) # 1*n
  K01 = t(K0%*%(-ystar)) # 1*n
  K20 = rowSums(K0*(lam.u)^2) # 1*n
  K11 = t((K0*lam.u)%*%(-ystar)) # 1*n
  ###### bug?
  eta = (K20*K01-K10*K11)/(K10*K10-K00*K20+eps) # 1*n #by my formula
  slope = (K11*K00-K01*K10)/(K10*K10-K00*K20+eps)
  
  negid = which(slope<0)
  #negid = 1
  # if (length(negid)>0){for (i in 2:length(eta)){
  #   if(eta[i]<=eta[i-1]) eta[i] = eta[i-1]+sval*(lam[i]-lam[i-1])
  #   slope[i] = sval
  # }}
  
  # if(sum(diff(t(eta))<0)>0){for (i in 2:length(eta)){
  #   if(eta[i]<=eta[i-1]) eta[i] = eta[i-1]+sval*(lam[i]-lam[i-1])
  #   slope[i] = sval
  # }}
  
  for (i in 2:length(eta)){
    if(eta[i] < eta[i-1]){
      eta[i] = eta[i-1]+sval*(lam[i]-lam[i-1])
      slope[i] = sval
    }
    #print(eta[i])
  }
  
  #if(sum(diff(t(eta))<0)>0)print(sum(diff(t(eta))<0))
  
  eta[ord] = eta
  slope[ord] = slope
  return(list(a = t(eta), b = t(slope)))
}


## solve the parameters
SlovePara <- function(y,x,z=c(),a,b,alpha){
  # x: n*p
  # z: n*q
  # alpha: p*1
  # a: n*1
  # b: n*1
  p = ncol(x)
  lam = x%*%alpha
  ytilde = y-a+b*lam
  ############### bug
  x_b = matrix(rep(b,p),ncol=p)*x
  if(length(z)!=0)mat = cbind(x_b,z) #n*(p+q)
  else mat = x_b
  #fit <- lm(ytilde~0+mat)
  #ksi <- coef(fit)
  #print(mat)
  sing = try(solve(t(mat)%*%mat),silent=T)
  if(inherits(sing,"try-error")==1){
    #print(b)
    return(list(alpha = alpha, beta = beta))
  }
  ksi = solve(t(mat)%*%mat)%*%t(mat)%*%ytilde
  if(length(z)==0){
      alpha_new = ksi
      beta_new = 0
    }else{
      alpha_new = ksi[1:p,1,drop=F]
      beta_new = ksi[-(1:p),1,drop=F]
    }
  return(list(alpha = alpha_new, beta = beta_new))
}

## Solve the parameters with penalty
PenalPara <- function(y,x,z=c(),a,b,alpha, pen.lam){
  # x: n*p
  # z: n*q
  # alpha: p*1
  # a: n*1
  # b: n*1
  p = ncol(x)
  lam = x%*%alpha
  ytilde = y-a+b*lam
  x_b = matrix(rep(b,p),ncol=p)*x
  if(length(z)!=0){
    mat = cbind(x_b,z) 
    q = ncol(z)}
  else{
    mat = x_b
    q = 0
  }
  
  pen.mat = c(rep(1,p),rep(0,q))
  fit.lasso = glmnet::glmnet(mat,ytilde,family = "gaussian",alpha=1,lambda = pen.lam,penalty.factor = pen.mat,intercept = F)
  
  ksi = coef(fit.lasso)[-1]
  if(length(z)==0){
    alpha_new = matrix(ksi,ncol=1)
    beta_new = 0
  }else{
    alpha_new = matrix(ksi[1:p],ncol=1)
    beta_new = matrix(ksi[-(1:p)],ncol=1)
  }
  return(list(alpha = alpha_new, beta = beta_new))
}

## alternating estimation iterating steps
PLSIM_Est_mon <- function(y,x,z=c(),h,alpha0,beta0=0,tol=1e-7,max.iter=10000,sval=1e-4,pen.lam = 0){
  n = nrow(x)
  alpha = alpha0
  beta = beta0
  
  iter=1
  diff=Inf
  diff.rec=c()
  loss.rec=c()
  
  alpha_list = c()
  converge = TRUE
  while(iter<max.iter & diff>tol){
    if(iter==max.iter-1)converge = FALSE
    res = CalEta_mon(y,x,z,alpha,beta,h,sval) 
    a = res$a #n*1
    b = res$b #n*1
    
    if(pen.lam == 0 )para.res = SlovePara(y,x,z,a,b,alpha)
    else{
      para.res = PenalPara(y,x,z,a,b,alpha, pen.lam)
    }
    
    alpha_new = para.res$alpha
    #if(alpha_new[1]<0)alpha_new = -alpha_new
    alpha_new = alpha_new/norm(alpha_new,"2") # standardize
    beta_new = para.res$beta
    #plot(x%*%alpha,a,main=iter)
    #plot(beta,main=iter,type="o")
    
    
    #if(length(z)!=0)diff.rec[iter]=norm(alpha_new-alpha,"2")+norm(beta_new-beta,"2") #l2 norm
    #else diff.rec[iter]=norm(alpha_new-alpha,"2")
    
    if(length(z)!=0)diff.rec[iter]=norm(alpha_new-alpha)+norm(beta_new-beta) #l1 norm
    else diff.rec[iter]=norm(alpha_new-alpha)
    

    diff=abs(diff.rec[iter])
    #plot(x%*%alpha,a,main=iter,ylim=c(-max(abs(range(y))),0))
    iter=iter+1
    alpha = alpha_new
    beta = beta_new
    if(iter>(max.iter-10))alpha_list <- cbind(alpha_list,alpha)
    
    
  }
  
  ####last step
  res = CalEta_mon(y,x,z,alpha,beta,h,sval)
  a = res$a #n*1
  b = res$b #n*1
  if(length(z)==0){ 
    mse = mean((y-a)^2)
    df = sum(alpha!=0) 
    }
  else{
    mse = mean((y-a-z%*%beta)^2)
    df = sum(alpha!=0)+length(beta)
  }
  BIC = log(mse)+log(n)/n*df
  return(list(alpha=alpha,beta = beta,eta = a, mse = mse, BIC = BIC, diff=diff.rec,iter = iter,h=h,x=x,y=y,z=z,last_alpha = alpha_list,converge=converge))

}

