if (getRversion() < "4.4.2") {
  stop("Please use R version 4.4.2.")
}

required_packages <- list(
  "tidyverse" = "2.0.0",
  "glmnet" = "4.1-8",
  "lars" = "1.3",
  "MASS" = "7.3-60",
  "latex2exp" = "0.9.6",
  "mnormt" = "2.1.1",
  "ggpubr" = "0.6.0",
  "viridis" = "0.6.4",
  "hrbrthemes" = "0.8.7",
  "readxl" = "1.4.3",
  "reshape2" = "1.4.4",
  "scales" = "1.3.0"
)

install_or_update <- function(pkg, min_version) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  } else {
    installed_version <- packageVersion(pkg)
    if (installed_version < min_version) {
      install.packages(pkg)
    }
  }
}

# Install or update each package
for (pkg in names(required_packages)) {
  install_or_update(pkg, required_packages[[pkg]])
}

# Load all required packages
lapply(names(required_packages), library, character.only = TRUE)

cat("All required packages are installed and loaded.\n")

get_xi <- function(w, h, n){
  # given omega, n and hkk, return xi
  ans = n*(1-w)/(n-1+w-n*(1-w)*h)
  abs(ans)
}

get_w <- function(xi, h, n){
  # reverse function of get_xi, 
  # given xi, hkk and n, return omega
  (n-xi*(n-1-n*h))/(n+(1+n*h)*xi)
}

centralize <- function(x, r = F){
  # centralize x with mean 0 and sd 1/sqrt(n-1)
  n = dim(x)[1]
  meanx <- drop(rep(1,n) %*% x)/n
  x = scale(x, meanx, FALSE)
  normx = sqrt(drop(rep(1,n) %*% (x^2)))
  if (r){
    return(list(m = meanx, d = scale(x,FALSE,normx), v = normx)) 
  }
  scale(x, FALSE, normx) 
}

# output the solution path of observation k when lambda is fixed 
beta_path <- function(X,y,k=1,lambda=50, plot=0, lb = 0){
  #' input:  X               matrix n by p      design matrix
  #'         y               vector n by 1      response vector
  #'         k               integer            observation of interest
  #'         lambda          float > 0          penalty parameter
  #'         plot            0                  no plot 
  #'                         1                  plot approximate 
  #'                         2                  plot exact solution path
  #'         lb              float < 1          the lower bound of omega 
  #' output: w_path          vector b by 1      a vector of breakpoints
  #'         hkk_path        vector b by 1      leverages at each breakpoint
  #'         beta_path       matrix b by p      beta estimate at each breakpoint
  #'         s_path          matrix b by p      beta hat's sign at each breakpoint
  #'         beta0_path      vector b by 1      beta0 at each breakpoint 
  #'         l1norm          vector b by 1      l1 norm of beta at each breakpoint  
  
  X = centralize(X)
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  # obtain lasso solution
  obj = lars(X,y,type='lasso')
  beta_hat = as.vector(predict.lars(obj,mode='lambda', s=lambda, type = 'coefficients')$coefficients)
  ybar = mean(y)
  # record the sign of each covariate 
  s = sign(beta_hat)
  # active set 
  A = s!=0
  
  # derivative of f in terms of covariates in nonactive set
  d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X%*%beta_hat
  if (sum(A)==0){
    XX = matrix(0,nrow=p,ncol=p)
    hk = as.matrix(rep(0,n),nrow=n)
  }
  else{
    XX = solve(t(X[,A])%*%X[,A])
    hk = X[,A]%*%XX%*%X[k,A]
  }
  
  # current predictor of yk 
  yk_hat = drop(ybar + X[k,] %*% beta_hat)
  
  # eazy case
  if (lambda == 0){
    X2 = cbind(1,X)
    beta1 = (solve(t(X2)%*%X2 - X2[k,]%*%t(X2[k,]))%*%(t(X2)%*%y - X2[k,]*y[k]))
    beta10 = beta1[1]
    beta1 = beta1[2:(p+1)]
    ind = c(which((beta_hat>0) & (beta1<0)),  which((beta_hat<0) & (beta1>0)))
    if(plot){
      coln = colnames(X)
      plot_helper0 <- function(w,ind){
        A = get_xi(w,hk[k],n)%*%t(XX%*%X[k,]*(yk_hat - y[k]))
        t(A[,ind]) + (XX%*%(t(X)%*%y))[ind]
      }
      if (length(ind)>0){
        if (length(ind) == 1){
          xx = data.frame(coef = coln[ind])
          print(ggplot()+xlim(lb,1) + xlab(TeX("$\\omega$")) + ylab(TeX("$\\beta$")) + 
                  geom_function(data = xx, fun = plot_helper0, args = list(ind),
                                aes(color = coef))+ 
                  ggtitle(paste("Solution path for Case",toString(k)))+
                  theme(plot.title = element_text(hjust = 0.5)))
        }
        else{
          fig = ggplot()+xlim(lb,1)
          for (j in 1:length(ind)){
            xx = data.frame(coef = coln[ind[j]])
            fig = fig + geom_function(data = xx, fun = plot_helper0, args = list(ind[j]),
                                      aes(color = coef)) 
          }
          print(fig + xlab(TeX("$\\omega$")) + ylab("$\\beta$") + ggtitle(paste("Solution path for Case",toString(k))) +
                  theme(plot.title = element_text(hjust = 0.5)))
        }
      }
      else{
        print(paste('there is no sign change for case',toString(k)))
      }
    }
    if(hk_){
      return(list(w_path = c(1,0),hkk_path = c(hk[k], hk[k]), 
                  beta_path = rbind(beta_hat, beta1), s_path = rbind(s,s = sign(beta1)), beta0_path = c(ybar, beta10),
                  hk_path = cbind(hk,hk), l1norm = sum(abs(beta_hat))))
    }else{
      return(list(w_path = c(1,0),hkk_path = c(hk[k], hk[k]), 
                  beta_path = rbind(beta_hat, beta1), s_path = rbind(s,s = sign(beta1)), 
                  beta0_path = c(ybar, beta10),l1norm = sum(abs(beta_hat))))
    }
  }
  
  hk_path = c()
  # beta path records beta hat's value at each breakpoint
  beta_path = c(beta_hat)
  # so does beta0_path records intercept
  beta0_path = c(ybar)
  # and sign change
  s_path = c(s)
  # and omega value at each breakpoint
  w_path = c()
  hkk_path = c()
  w = 1
  while (T){
    hk_path = cbind(hk_path,hk)
    w_path = c(w_path, w)
    hkk_path = c(hkk_path, hk[k])
    xi = get_xi(w, hk[k], n)
    bias = yk_hat - y[k]
    if (sum(A) == 0){
      xi_cand1 = c()
    }
    else{
      slope_beta = XX%*%X[k,A]*bias
      xi_cand1 = -beta_hat[A]/slope_beta
    }
    slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
    # xi candidates
    # xi_cand1 = -beta_hat[A]/slope_beta
    xi_cand2p = (-lambda-d_hat)/slope_d
    xi_cand2m = (lambda-d_hat)/slope_d
    xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
    xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
    ind = which(xi_cand == xi_cand0)
    
    # update beta
    if (sum(A) > 0){
      beta_hat[A] = beta_hat[A] + min(get_xi(lb,hk[k],n),xi_cand0)*slope_beta
    }
    beta_path = rbind(beta_path, beta_hat)
    beta_hat0 = ybar + min(get_xi(lb,hk[k],n),xi_cand0) * bias/n
    beta0_path = c(beta0_path, beta_hat0)
    
    # if the xi is off the bound, stop the algorithm
    if (xi_cand0 > get_xi(lb,hk[k],n)){
      w_path = c(w_path, lb)
      hkk_path = c(hkk_path, hk[k])
      s_path = rbind(s_path, s)
      hk_path = cbind(hk_path, hk)
      break
    }
    # if not, locate the covariate (go to or leave the active set)
    else if(ind<=sum(A)){
      coef = which(cumsum(A)==ind)[1]
      A[coef] = F
      s[coef] = 0
      # check if active set is empty 
      if (sum(A) == 0){
        hk_path = cbind(hk_path, 0,0)
        w_path = c(w_path, get_w(xi_cand0, hk[k], n), lb)
        hkk_path = c(hkk_path, 0, 0)
        s_path = rbind(s_path, s, s)
        beta_path = rbind(beta_path, beta_hat)
        beta0_path = c(beta0_path, ybar + get_xi(lb,0,n) * (ybar-y[k])/n)
        break
      }
    }
    else if(ind>p){
      coef = which(cumsum(1-A)==(ind-p))[1]
      A[coef] = T
      s[coef] = 1
    }
    else{
      coef = which(cumsum(1-A)==(ind-sum(A)))[1]
      A[coef] = T
      s[coef] = -1
    }
    
    s_path = rbind(s_path, s)
    # update omega, XX, beta_hat, d_hat
    w = get_w(xi_cand0, hk[k], n)
    XX = solve(t(X[,A])%*%X[,A])
    beta_hat[A] = XX%*%(t(X[,A])%*%y-lambda*s[A])
    beta_hat[!A] = 0
    d_hat = t(X[,!A])%*%y - t(X[,!A])%*%X%*%beta_hat
    # y_hat = ybar + X%*%beta_hat
    hk = X[,A]%*%XX%*%X[k,A]
    yk_hat = drop(ybar + X[k,] %*% beta_hat)
  }
  if(plot){
    plot_helper <-function(x, df){
      i = findInterval(-x, -w_path, rightmost.closed = T)
      beta1 = df[i]
      beta2 = df[i+1]
      w1 = w_path[i]
      w2 = w_path[i+1]
      hkk = hkk_path[i]
      beta1+(beta2-beta1)*(get_xi(x, hkk, n) - get_xi(w1, hkk, n))/(get_xi(w2, hkk, n) - get_xi(w1, hkk, n))
    }
    if (is.null(colnames(X))){
      coln = 1:p
    }
    else{
      coln = colnames(X)
    }
    num_z = apply(beta_path, 2, function(c) sum(abs(c)< 1e-10))
    ind = which(num_z>0 & num_z<length(w_path))
    # ind = which(num_z>=0)
    if (length(ind)>0){
      if(plot ==1){
        df = cbind(beta_path[,ind], w_path)
        colnames(df) = c(coln[ind], 'w')
        df = as_tibble(df) %>% gather('coef', 'val', -w)
        print(ggplot(df, aes(x = w, y = val, group=coef, color = coef))+
                geom_line()+ggtitle(paste("Approx Solution path for Case",toString(k)))+
                theme(plot.title = element_text(hjust = 0.5)) + ylab(TeX("$\\beta$")))
      }
      else{
        if (length(ind) == 1){
          xx = data.frame(coef = coln[ind])
          print(ggplot()+xlim(lb,1) + xlab(TeX("$\\omega$")) + ylab(TeX("$\\beta$")) + 
                  geom_function(data = xx, fun = plot_helper, args = list(beta_path[,ind]), aes(color = coef))+ 
                  ggtitle(paste("Solution path for Observation",toString(k)))+
                  theme(plot.title = element_text(hjust = 0.5)))
        }
        else{
          df = beta_path[,ind]
          colnames(df) = coln[ind]
          fig = ggplot()+xlim(lb,1)
          for (j in 1:length(ind)){
            xx = data.frame(coef = colnames(df)[j])
            fig = fig + geom_function(data = xx, fun = plot_helper, 
                                      args = list(df[,j]), aes(color = coef)) 
          }
          fig = fig + labs(color = "Feature")
          print(fig + xlab(TeX("Case weight $\\omega$")) + ylab(TeX("Coefficient estimate $\\hat{\\beta}$"))+ 
                  theme(plot.title = element_text(hjust = 0.5))+theme(panel.background = element_rect(fill = "white"),
                                                                      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)))
        }
      }
    }
    else{
      print(paste('there is no sign change for case',toString(k)))
    }
  }
  return(list(w_path = w_path, hkk_path = hkk_path, beta_path = beta_path, s_path = s_path, beta0_path = beta0_path,
              l1norm = sum(abs(beta_path[1,]))))
}


qr_delete = function(Q,R,col){
  # function test 
  # A = matrix(c(1, 2,3, 4,
  #              0, 4, 1, 5,
  #              0, 2, 2,9,
  #              1,  0,0,11), 
  #            nrow = 4, byrow = TRUE)
  # QR = qr(A)
  # PP = qr_delete(qr.Q(QR),qr.R(QR),2)
  # PP$Q%*%PP$R
  if (dim(R)[1]<col){
    warning("no such feature included", call. = FALSE)
  }
  p = dim(R)[1]
  R = R[,-col,drop=F]
  if (col==p){
    return(list(Q=Q[,-p,drop=F],R=R[-p,,drop=F]))
  }
  if (col <= (p-2)){
    for (i in col:(p-2)){
      a = R[i,i]
      b = R[i+1,i]
      rot = matrix(c(a/sqrt(a**2+b**2), -b/sqrt(a**2+b**2), 
                     b/sqrt(a**2+b**2), a/sqrt(a**2+b**2)), nrow=2)
      Q[,i:(i+1)] = Q[,i:(i+1)]%*%t(rot)
      R[i:(i+1),] = rot%*%R[i:(i+1),]
    }
  }
  temp = sqrt(R[p,p-1]**2 + R[p-1,p-1]**2)
  Q[,(p-1)] = (R[p-1,p-1] * Q[,(p-1)] + R[p,p-1] * Q[,p])/temp
  R[p-1,p-1] = temp
  return(list(Q=Q[,-p,drop=F],R=R[-p,,drop=F]))
}


qr_insert = function(Q,R,col,v){
  # test code
  # A = matrix(c(1, 2,3,
  #              0, 4, 1,
  #              0, 2, 2,
  #              1,  0,0),
  #            nrow = 4, byrow = TRUE)
  # QR = qr(A)
  # PP = qr_insert(qr.Q(QR),qr.R(QR),1,c(1,2,3,4))
  # PP$Q%*%PP$R
  p = dim(R)[1]
  if (is.null(p)){
    nor = sqrt(sum(v**2))
    return(list(Q= as.matrix(v/nor,length(v),1), R = as.matrix(nor,1,1)))
  }
  if (col>(1+p)){
    warning("there is not enough features before", call. = FALSE)
  }
  if (col == (1+p)){
    coor = t(Q)%*%v
    res = v - Q%*%coor
    nor = sqrt(sum(res**2))
    return(list(Q=cbind(Q, res/nor), R = cbind(rbind(R,0), c(coor, nor))))
  }
  new_coor = t(Q)%*%v
  new_q = v - Q%*%new_coor
  nor = sqrt(sum(new_q**2))
  new_q = new_q/nor
  Q = cbind(Q, new_q)
  if (col == 1){
    R = rbind(cbind(new_coor,R), c(nor,rep(0,p)))
  }
  else{
    R = rbind(cbind(R[,1:(col-1)], new_coor, R[,col:p]), c(rep(0, col-1), nor, rep(0,p-col+1)))
  }
  
  for (i in (p+1):(col+1)){
    a = R[i-1,col]
    b = R[i,col]
    rot = matrix(c(a/sqrt(a**2+b**2), -b/sqrt(a**2+b**2), 
                   b/sqrt(a**2+b**2), a/sqrt(a**2+b**2)), nrow=2)
    Q[,(i-1):i] = Q[,(i-1):i]%*%t(rot)
    R[(i-1):i,] = rot%*%R[(i-1):i,]
  }
  return(list(Q = Q, R = R))
}


#' calculate case influence for Lasso regression
#' @param X           input matrix, of dimension nobs x nvars; each row is an observation vector. 
#'                    Requirement: nvars >1; in other words, x should have 2 or more columns.
#' @param y           response variable. 
#' @param k           an integer, or vector of integers representing the index/indices 
#'                    of observations that are of interest. By default, `CookDisLasso()` considers
#'                    all observations. 
#' @param fineness    number of grid points for fractions being considered from 0 to 1. 
#'                    The default is set to be 100. Then we take fractions 0, 0.01, ..., 0.99
#' @param s           a value, or vector of values, indexing the penalty levels to consider.
#'                    Its values depends on the mode= argument. By default (mode="fraction"), 
#'                    s should take on values between 0 and 1. 
#' @param mode        Mode="fraction", then s should be a number between 0 and 1, and it refers 
#'                    to the ratio of the L1 norm of the coefficient vector, relative to the 
#'                    norm at the full LS solution if n>p or the minimum L1 norm linear 
#'                    interpolation solution if n<=p. Mode="norm" means s refers to the L1 norm 
#'                    of the coefficient vector. Mode="lambda" uses the lasso regularization 
#'                    parameter for s. Abbreviations allowed.
#' @param threshold   If TRUE, CookDisLasso prints the threshold for influential points. This 
#'                    can only set to be TRUE if k is 1:n (by default).
#' @return \item{CD_Mat}{matrix of cook's distance at each fraction}
#'         \item{Lambda_list}{a vector of lambdas at each fraction}
#'         \item{fraction}{a vector of fractions based on required fineness}
#'         \item{threshold_table}{threshold at each fraction if threshold = TRUE}
#'         \item{beta_hat_table}{a matrix of betahats at each fraction}
#'         \item{index}{samples of interest. If k is specified in input, index = k, otherwise 1:n}
#'         \item{lev_history}{record the leverages of samples of interest in each grid point specified in fineness}
#'         \item{denom}{estimated error variance}
#' @description
#' Using case-weight lasso model and the solution path algorithm, this function calculates 
#' the case influence for the lasso for all observations(for specified subset), all fractions(for specified lambdas) and also gives 
#' case influence graph. 
#' @details
#' 1. threshold needs case influence of all observations to calculate. So threshold should be set to FALSE if the functions only calculate 
#' case influence for a subset of observations,i.e. k is not NULL.
#' 2. fineness automatically chooses a list of lambda candidates from 0 to maximum lambda that penalizes all beta to 0. If lambda candidate 
#' is pre-specified, i.e. lambda is not NULL, fineness should be set to NULL.
#' @seealso plot method
#' @examples
#' library(lars)
#' data(diabetes)
#' attach(diabetes)
#' CookDisLasso(x,y)
#' detach(diabetes)
#' 
#' set.seed(10)
#' x = matrix(rnorm(200*50),nrow=200)
#' y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(200)
#' obj1 = CookDisLasso(x,y, k=1:10, s=c(0.1,0.2), mode = "fraction", threshold = FALSE)
#' plot(obj1, 'resid', 1)
#' obj2 = CookDisLasso(x,y, fineness=40, threshold = TRUE)
#' plot(obj2, 'case-influence-graph')
#' @export
CookDisLasso <- function(X, y, k, fineness, s, mode = c("fraction", "norm", "lambda"),
                         threshold = TRUE){  
  
  num_of_update = 0
  mode <- match.arg(mode)
  if (!missing(k) && threshold) {
    warning("Unable to calculate threshold without all observations.", call. = FALSE)
    threshold = FALSE
  }
  if (!missing(s) && !missing(fineness)) {
    warning("'lambda' and 'fineness' cannot be specified together.", call. = FALSE)
    fineness = NULL
  }
  if (missing(fineness) && missing(s)) {
    fineness = 100
  }
  
  X = centralize(X)
  n = dim(X)[1]
  p = dim(X)[2]
  
  # record all gram inverses we come across during the calculation 
  recorder = list()
  recorder[[paste(rep(0,p),collapse = '')]] = list(Q=0,R=0,Rinv=0)
  
  lambda_max = max(abs(t(X)%*%y))
  
  obj = lars(X,y,type='lasso',use.Gram = !(n<p | p>500))
  nbeta1 = drop(abs(obj$beta) %*% rep(1, p))
  ybar = mean(y)
  beta_hat_table = c()
  if (!missing(fineness)){
    fraction = seq(0.99,0,length.out=fineness)
    for (f in fraction){
      beta_hat = as.vector(predict.lars(obj,mode='fraction', s=f, type = 'coefficients')$coefficients)
      beta_hat_table = cbind(beta_hat_table, beta_hat)
    }
    nbeta2 = drop(rep(1, p) %*% abs(beta_hat_table))
    l_list = approx(nbeta1, c(obj$lambda,0), xout = nbeta2)$y
  } 
  else if (mode == "lambda"){    
    l_list = sort(s[s<lambda_max & s>=0])
    if (length(l_list) == 0){
      # warning("no valid lambda provided.", call. = FALSE)
      l_list = c(lambda_max/2)
    }
    for (l in l_list){
      beta_hat = as.vector(predict.lars(obj,mode='lambda', s=l, type = 'coefficients')$coefficients)
      beta_hat_table = cbind(beta_hat_table, beta_hat)
    }
    fraction = approx(c(obj$lambda,0), nbeta1/nbeta1[length(nbeta1)], xout = l_list)$y
  }
  else if (mode == "fraction"){
    fraction = sort(s[s<=1 & s>=0], decreasing = T)
    for (f in fraction){
      beta_hat = as.vector(predict.lars(obj,mode='fraction', s=f, type = 'coefficients')$coefficients)
      beta_hat_table = cbind(beta_hat_table, beta_hat)
    }
    nbeta2 = drop(rep(1, p) %*% abs(beta_hat_table))
    l_list = approx(nbeta1, c(obj$lambda,0), xout = nbeta2)$y
  }
  else {
    max_norm = sum(abs(as.vector(predict.lars(obj,mode='fraction', s=1, type = 'coefficients')$coefficients)))
    s = s/max_norm
    fraction = sort(s[s<=1 & s>=0], decreasing = T)
    for (f in fraction){
      beta_hat = as.vector(predict.lars(obj,mode='fraction', s=f, type = 'coefficients')$coefficients)
      beta_hat_table = cbind(beta_hat_table, beta_hat)
    }
    nbeta2 = drop(rep(1, p) %*% abs(beta_hat_table))
    l_list = approx(nbeta1, c(obj$lambda,0), xout = nbeta2)$y
  }
  
  CD_Mat = c()
  
  if (missing(k)){
    indices_vec = 1:n
  }
  else{
    indices_vec = k
  }
  
  if (n<=p){
    denom = 1
  } 
  else{
    denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)
  }
  
  if (0 %in% l_list){
    starter = 2
    if (n<=p){
      cd0 = rep(NA,length(indices_vec))
    } else{
      cd0 = cooks.distance(lm(y~X))[indices_vec]
    }
  } else{
    starter = 1
  }
  
  lev_history = c()
  
  if ((length(fraction)==1) & (starter==2)){
    return('do regular cook distance')
  }
  
  
  for (i in starter:length(fraction)){
    lambda = l_list[i]
    
    beta_hat_backup = beta_hat_table[,i]
    s_backup = sign(beta_hat_backup)
    A_backup = s_backup!=0
    # derivative of f in terms of covariates in nonactive set
    d_hat_backup = t(X[,!A_backup])%*%y - t(X[,!A_backup])%*%(X[,A_backup,drop=F]%*%beta_hat_backup[A_backup,drop=F])
    
    A_id = paste(A_backup*1,collapse = '')
    if (A_id %in% names(recorder)){
      QR_backup = recorder[[A_id]]
    }
    else{
      qr_result = qr(X[,A_backup])
      Rinv = backsolve(qr.R(qr_result), diag(dim(qr.R(qr_result))[1]))
      recorder[[A_id]] = list(Q=qr.Q(qr_result), R = qr.R(qr_result), Rinv = Rinv)
      QR_backup = recorder[[A_id]]
    }
    if (sum(A_backup)==0){
      h_backup = matrix(0,ncol=n,nrow=n)
    }
    else{
      h_backup = recorder[[A_id]]$Q%*%t(recorder[[A_id]]$Q)
      # O(n^2 p)--> thats for n samples 
    }
    lev_history = cbind(lev_history, diag(h_backup)+1/n)
    CD_list = c()
    
    y_tilde = X[,A_backup,drop=F]%*%beta_hat_backup[A_backup] + ybar
    
    for (k in indices_vec){
      beta_hat = beta_hat_backup
      s = s_backup
      A = A_backup
      
      d_hat = d_hat_backup
      QR = QR_backup
      hk = h_backup[,k]
      
      # current predictor of yk 
      yk_hat = drop(y_tilde[k])
      w = 1
      while (T){
        xi = get_xi(w, hk[k], n)
        bias = yk_hat - y[k]
        if (sum(A) == 0){
          xi_cand1 = c()
        }
        else{
          slope_beta = QR$Rinv%*%QR$Q[k,]*bias
          xi_cand1 = -beta_hat[A]/slope_beta
        }
        if (sum(A) == (n-1)){
          slope_d = rep(0, p-sum(A))
          xi_cand2p = rep(Inf, p-sum(A))
          xi_cand2m = rep(Inf, p-sum(A))
        }
        else{
          slope_d = (X[k,!A] - t(X[,!A])%*%hk)*bias
          xi_cand2p = (-lambda-d_hat)/slope_d
          xi_cand2m = (lambda-d_hat)/slope_d
        }
        # xi candidates
        xi_cand = c(xi_cand1, xi_cand2p, xi_cand2m)
        xi_cand0 = min(xi_cand[xi_cand>(xi+0.00001)],Inf)
        ind = which(xi_cand == xi_cand0)[1]
        
        # update beta
        if (sum(A) > 0){
          beta_hat[A] = beta_hat[A] + min(get_xi(0,hk[k],n),xi_cand0)*slope_beta
        }
  
        beta_hat0 = ybar + min(get_xi(0,hk[k],n),xi_cand0) * bias/n
        
        # if the xi is off the bound, stop the algorithm
        if (is.na(xi_cand0) | xi_cand0 > get_xi(0,hk[k],n)){
          break
        }
        # if not, locate the covariate (go to or leave the active set)
        else if(ind<=sum(A)){
          coef = which(cumsum(A)==ind)[1]
          # print(c('-',coef))
          A[coef] = F
          s[coef] = 0
          # check if active set is empty 
          if (sum(A) == 0){
            hk[k] = 0
            break
          }
        }
        else if(ind>p){
          coef = which(cumsum(1-A)==(ind-p))[1]
          
          A[coef] = T
          s[coef] = 1
        }
        else{
          coef = which(cumsum(1-A)==(ind-sum(A)))[1]
          
          A[coef] = T
          s[coef] = -1
        }
        num_of_update = num_of_update + 1
        
        # update omega, XX, beta_hat, d_hat
        w = get_w(xi_cand0, hk[k], n)
        A_id = paste(A*1,collapse = '')
        if (A_id %in% names(recorder)){
          QR = recorder[[A_id]]
        }
        else{
          if (A[coef]){
            QR = qr_insert(QR$Q, QR$R, sum(A[1:coef]), X[,coef])
          }
          else{
            QR = qr_delete(QR$Q, QR$R, sum(A[1:coef])+1)
          }
          
          QR$Rinv = backsolve(QR$R, diag(dim(QR$R)[1]))
          recorder[[A_id]] = QR

        }
  
        beta_hat[A] = QR$Rinv%*%(t(QR$Rinv)%*%( - lambda * s[A] + t(X[,A,drop=F])%*%y))
        beta_hat[!A] = 0
        
        d_hat = t(X[,!A])%*%y - t(X[,!A])%*%(X[,A,drop=F]%*%beta_hat[A,drop=F])
        hk =  QR$Q%*%t(QR$Q[k, ,drop=F])
        yk_hat = drop(ybar + X[k,A] %*% beta_hat[A])
      }
      
      #--------------------
      
      A = s!=0
      A_id = paste(A*1,collapse = '')
      if (sum(A) == 0){
        y_last = ybar + get_xi(0, hk[k], n)/n*(ybar - y[k])
      }
      else{
        mid = QR$Q%*%(t(QR$Rinv)%*%(t(X[,A])%*%y-lambda*s[A])) 
        y_hat = mid + ybar
        y_last = y_hat + get_xi(0, hk[k], n)*(QR$Q%*%t(QR$Q[k,,drop=F])+1/n)*(y_hat[k] - y[k])
      }
      
      CD_list = c(CD_list, sum((y_last - y_tilde)**2))
    }
    CD_Mat = cbind(CD_Mat, CD_list)
  }
  
  CD_Mat = CD_Mat/denom
  if (starter == 2){
    CD_Mat = cbind(cd0,CD_Mat)
  }
  
  colnames(CD_Mat) <- NULL
  colnames(beta_hat_table) <- NULL
  rownames(threshold) <- NULL
  ans = list(CD_Mat = CD_Mat, Lambda_list = l_list, fraction = fraction, beta_table = beta_hat_table,
             index = indices_vec, lev_history = lev_history, denom = denom/(p+1), X=X, y=y, 
             avg_update = num_of_update/length(fraction)/length(indices_vec))
  
  if (threshold){
    threshold_table <- matrix(0, n, dim(CD_Mat)[2])
    for (j in 1:dim(CD_Mat)[2]) {
      # Extract the current column
      col_vals <- CD_Mat[, j]
      
      # Pre-calculate the sum and sum of squares for the column
      col_sum <- sum(col_vals)
      col_sum_sq <- sum(col_vals^2)
      
      # Loop through each row to calculate leave-one-out variance
      for (i in 1:n) {
        # Calculate the sum and sum of squares without the i-th entry
        loo_sum <- col_sum - col_vals[i]
        loo_sum_sq <- col_sum_sq - col_vals[i]^2
        
        # Calculate the leave-one-out variance using the pre-calculated values
        threshold_table[i, j] <- sqrt((loo_sum_sq - (loo_sum^2 / (n - 1))) / (n - 1) / 2)*qchisq(0.95,1)
      }
    }
    ans$threshold = threshold_table
  }
  class(ans) <- "CookDisLasso"
  invisible(ans)
}
