noise = 10
rep = 100

library(PLSiMCpp)
library(glmnet)
library(mnormt)
library(pls)
library(np)
library(scatterplot3d)
library(ks)
library(fields)
library(isotone)
library(Iso)
library(numDeriv)

library(locfit)
library(Rcpp)
library(MASS)

source("bwfun.R")
source("finalbetas_cvh_altconv_updated.R")
source("pe_cvh_altconv.R")
source("summary_measures_cvh_updated.R")

source("cv_mon.R")
source("MonSIM.R")

source("code_ispline.R")
source("predict_ispline.R")
source("predict_foster.R")

source('fpigsim.R')
source('fisher.R')
#source('predict_sse_ese.R')
#sourceCpp("SSE.cpp")
#sourceCpp("ESE.cpp")
#sourceCpp("LSE.cpp")

n = 200
p = 3


train_mse = matrix(nrow=rep,ncol=6)
y_mse_true = matrix(nrow=rep,ncol=6)
res_3_l = list()
res_3_l_foster = list()
res_3_l_liang = list()
res_3_l_wan = list()
#res_3_l_lm = list()
res_3_l_sse = list()
res_3_l_ese = list()
res_3_alpha = matrix(NA,nrow=rep,ncol=p)
res_3_alpha_foster = matrix(NA,nrow=rep,ncol=p)
res_3_alpha_liang = matrix(NA,nrow=rep,ncol=p)
res_3_alpha_wan = matrix(NA,nrow=rep,ncol=p)
res_3_alpha_sse = matrix(NA,nrow=rep,ncol=p)
res_3_alpha_ese = matrix(NA,nrow=rep,ncol=p)

time = matrix(NA,nrow=rep,ncol=6)
area.pp.sp = rep(NA,rep)
area.liang.sp = rep(NA,rep)
area.foster.sp = rep(NA,rep)
area.wan.sp = rep(NA,rep)
area.sse.sp = rep(NA,rep)
area.ese.sp = rep(NA,rep)

# preperation

## specific setting of Foster
beta.iter <- 500
approach <- "kern"
constrain <- 1
force1 <- 1
knowzero <- 0
l2loss <- 1
ss <- 100

##
raw_alpha = c(1,-1,1)
alpha = raw_alpha/norm(raw_alpha,"2")
set.seed(2236)
x = matrix(runif(n*p,-2,2),ncol = p)
train_id = sample(1:n,0.8*n)
test_id = setdiff(1:n, train_id)
x_train = x[train_id,]
x_test = x[test_id,]
lam_test = x_test%*%alpha
#y_test_true = 1/(exp(-lam_test)+1)
y_test_true = (lam_test)^3
lam_train = x_train%*%alpha


# area calculation
area.cal <- function(x.int,u){
  f <- approxfun(x.int,u)
  #area = integrate(function(x) abs(f(x)-1/(exp(-x)+1)),  quantile(lam_train,0.05),quantile(lam_train,0.95),stop.on.error = FALSE)$value
  area = integrate(function(x) abs(f(x)-x^3), -1.9,1.9,stop.on.error = FALSE)$value
  return(area)
}


for (i in 1:rep){
  flag = 0
  set.seed(2025+i)
  vec_noise = rnorm(0.8*n,sd=noise)
  #lam_train = x_train%*%alpha 
  #y_train = 1/(exp(-lam_train)+1)+vec_noise
  y_train = (lam_train)^3+vec_noise
  
  init = coef(lm(y_train~0+x_train))
  alpha0 = init[1:p]/norm(init[1:p],"2")
  
  # fit models
  t1 <- Sys.time()
  pp.bw = cv_bw_PLSIM_mon(y_train,x_train,alpha0=alpha0,tol=1e-5,max.iter=500,hlist=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),n.fold = 5,pen.lam = 0,sval=0.1)$best.bw
  res_3_pos = PLSIM_Est_mon(y_train,x_train,z=c(),pp.bw,alpha0,tol = 1e-5,sval =0.1,max.iter = 500)
  print(pp.bw)
  
  t2 <- Sys.time()
  liang.bw = plsim.bw(xdat=c(),ydat=data.frame(y_train),zdat=x_train,zeta_i = alpha0,bandwidthList=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7))$bandwidthBest
  res_3_liang = plsim.est(ydat=data.frame(y_train),zdat=x_train,h=liang.bw)
  print(liang.bw)
  
  t3 <- Sys.time()
  foster_alpha0 = init[1:p]/init[1]
  alw <- 1/(abs(foster_alpha0)^(3/5))
  res_3_foster = try(finalbetas_cvh_altconv(foster_alpha0,0,x_train,y_train,alw))
  if(inherits(res_3_foster, "try-error"))flag = 1
  if(flag==0)foster_beta = res_3_foster$original[-c(1:2)]
  
  t4 <- Sys.time()
  res_3_wan = Est.f.alpha(x_train, y_train, foster_alpha0, iter.lambda=100, iter.tot=200,
                          iter.alpha=30, crit.alpha=1e-4, d.knot=1, print=FALSE)

  t5 <- Sys.time()
  
  
  #save models
  res_3_l[[i]] = res_3_pos
  res_3_l_liang[[i]] = res_3_liang
  res_3_l_foster[[i]] = res_3_foster
  res_3_l_wan[[i]] = res_3_wan
  res_3_l_sse[[i]] = res_sse
  res_3_l_ese[[i]] = res_ese
  
  #save alpha
  res_3_alpha[i,] = res_3_pos$alpha
  res_3_alpha_liang[i,]=res_3_liang$zeta
  if(flag==0)res_3_alpha_foster[i,] = foster_beta/norm(foster_beta,"2")
  res_3_alpha_wan[i,] = res_3_wan$alpha/norm(res_3_wan$alpha,"2")

  
  area.pp.sp[i] = area.cal(res_3_pos$x%*%(res_3_pos$alpha),res_3_pos$eta)
  area.liang.sp[i] = area.cal(res_3_liang$Z_alpha,res_3_liang$eta)
  if(flag==0)area.foster.sp[i] = area.cal(res_3_foster$z/norm(foster_beta,"2"),res_3_foster$eta)
  area.wan.sp[i] = area.cal(sort(x_train%*%res_3_wan$alpha/norm(res_3_wan$alpha,"2")),res_3_wan$f)
  
  #predict and mse
  y_pred = pred_PLSIM_Est(res_3_pos,xnew = x_test)
  y_pred_liang = predict(res_3_liang,z_test = x_test)
  if(flag==0)y_pred_foster = predict_foster(z=res_3_foster$z,eta= res_3_foster$eta, beta=foster_beta,newx = x_test)
  y_pred_wan = predict_ispline(res_3_wan, Xnew = x_test)
  
  y_mse_true[i,1] = mean((y_pred-y_test_true)^2)
  y_mse_true[i,2] = mean((y_pred_liang-y_test_true)^2)
  if(flag==0)y_mse_true[i,3] = mean((y_pred_foster-y_test_true)^2)
  y_mse_true[i,4] = mean((y_pred_wan-y_test_true)^2)
  
  train_mse[i,1] <- res_3_pos$mse
  train_mse[i,2] <- res_3_liang$mse
  if(flag==0)train_mse[i,3] <- mean((sort(y_train)-res_3_foster$eta)^2)
  train_mse[i,4] <- res_3_wan$mse
  
  # time
  time[i,] = c(t2-t1,t3-t2,t4-t3,t5-t4,t6-t5,t7-t6)
  if(flag==1)time[i,3]=NA
  print(i)
  print(norm(res_3_alpha[i,]-alpha,"2"))
}


norm.cal <- function(est,alpha){
  norm.df <- rep(NA,nrow(est))
  for (i in 1:nrow(est)){
    norm.df[i] = norm(est[i,]-alpha,"2")
  }
  return(norm.df)
}

pp.alpha.norm <- norm.cal(res_3_alpha,alpha)
liang.alpha.norm <- norm.cal(res_3_alpha_liang,alpha)
foster.alpha.norm <- norm.cal(res_3_alpha_foster,alpha)
wan.alpha.norm <- norm.cal(res_3_alpha_wan,alpha)
#sse.alpha.norm <- norm.cal(res_3_alpha_sse,alpha)
#ese.alpha.norm <- norm.cal(res_3_alpha_ese,alpha)

sse.alpha.norm <- ese.alpha.norm <- NA

alpha.df = data.frame(proposed = pp.alpha.norm, liang = liang.alpha.norm, foster = foster.alpha.norm, wan = wan.alpha.norm,sse = sse.alpha.norm,ese = ese.alpha.norm, noise=noise)
area.sp.df = data.frame(proposed=area.pp.sp,liang=area.liang.sp,foster=area.foster.sp,wan=area.wan.sp,sse=area.sse.sp,ese=area.ese.sp,noise = noise)

y_mse_true <- data.frame(cbind(y_mse_true,rep(noise,rep)))
names(y_mse_true)<-c("proposed","liang","foster","wan","sse","ese","noise")

train_mse <- data.frame(cbind(train_mse,rep(noise,rep)))
names(train_mse)<-c("proposed","liang","foster","wan","sse","ese","noise")

time <- data.frame(cbind(time,rep(noise,rep)))
names(time)<-c("proposed","liang","foster","wan","sse","ese","noise")

write.csv(alpha.df,"output/alpha_cubic.csv",row.names = F)
write.csv(area.sp.df,"output/area_cubic.csv",row.names = F)
write.csv(y_mse_true,"output/test_mse_cubic.csv",row.names = F)
write.csv(train_mse,"output/train_mse_cubic.csv",row.names = F)
write.csv(time,"output/time_cubic.csv",row.names = F)
