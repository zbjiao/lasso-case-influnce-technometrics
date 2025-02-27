rm(list=ls())
source("caseweightlasso2024.R")

####################################################################################
#-------------------- figure 1 case influence graph mechanism ---------------------#
plot_example <- function(X,y,title){
  temp = CookDisLasso(X,y,k=10,fineness = 1000,threshold = F)  
  plot_table = as.data.frame(cbind(t(temp$CD_Mat), temp$fraction))
  colnames(plot_table) = c('val','fraction')
  p = ggplot(data=plot_table,aes(x=fraction,y=val)) +
    geom_line() + xlim(c(0,NA))+ xlab(TeX("$|coef|/max|coef|$"))+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none") +labs(y="Cook's distance")+ggtitle(title)+
    geom_line(data=df2,aes(x=fraction, y=threshold),linetype = "dotted", color = 'red1')+
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  return(p)
}

set.seed(1)
n = 10
root_x1 = rnorm(n,0,1)
root_x2 = rnorm(n,0,1)
root_error =  rnorm(n, 0, 1)

# scenario 2
x1 = root_x1
x1[10] = mean(root_x1[-10])
x2 = root_x2
x2[10] = 4
X = centralize(cbind(x1, x2))
y = X%*%c(4,1) + root_error
y[10] = mean(y[-10])

result = CookDisLasso(X,y,fineness = 1000)
df2 = data.frame(fraction = result$fraction, threshold = 
                   sqrt(apply((result$CD_Mat)[1:9,],2,var)/2)*qchisq(0.95,1))

p1 = plot_example(X,y,'Scenario II')


# scenario 3
x1 = root_x1
x1[10] = mean(root_x1[-10])
x2 = root_x2
x2[10] = 4
X = centralize(cbind(x1, x2))
y = X%*%c(4,1) + root_error
y[10] = lm(y[-10]~X[-10,])$coefficients%*%c(1,X[10,])
p2 = plot_example(X,y,'Scenario III')


# scenario 3
x1 = root_x1
x1[10] = mean(root_x1[-10])
x2 = root_x2
x2[10] = mean(root_x2[-10])
X = centralize(cbind(x1, x2))
y = X%*%c(4,1) + root_error
y[10] = 4
p3 =  plot_example(X,y,'Scenario I')

# scenario 4
x1 = root_x1
x1[10] = mean(root_x1[-10])+4
x2 = root_x2
x2[10] = mean(root_x2[-10])
X = centralize(cbind(x1, x2))
y = X%*%c(4,1) + root_error
y[10] = -4
p4 = plot_example(X,y,'Scenario IV')

ggarrange(p3+rremove('xlab'), 
          p1+rremove('ylab')+rremove('xlab'), 
          p2, 
          p4+rremove('ylab'), 
          ncol = 2, nrow = 2)

####################################################################################
#-------------------- figure 2 model selection example ---------------------------#

plot.lars <-
  function(x, xvar=c("norm","df","arc.length","step"), breaks = TRUE, 
           plottype = c("coefficients", "Cp"), 
           omit.zeros = TRUE, eps = 1e-10, ...)
  {
    object <- x
    plottype <- match.arg(plottype)
    xvar <- match.arg(xvar)
    coef1 <- object$beta	### Get rid of many zero coefficients
    if(x$type!="LASSO"&&xvar=="norm")# problem with discontinuity in norm
      coef1=betabreaker(x)
    stepid=trunc(as.numeric(dimnames(coef1)[[1]]))
    coef1 <- scale(coef1, FALSE, 1/object$normx)
    if(omit.zeros) {
      c1 <- drop(rep(1, nrow(coef1)) %*% abs(coef1))
      nonzeros <- c1 > eps
      cnums <- seq(nonzeros)[nonzeros]
      coef1 <- coef1[, nonzeros,drop=FALSE]
    }
    else cnums <- seq(ncol(coef1))
    s1<-switch(xvar,
               norm={
                 s1 <- apply(abs(coef1), 1, sum)
                 s1/max(s1)
               },
               df=object$df,
               arc.length=cumsum(c(0,object$arc.length)),
               step=seq(nrow(coef1))-1
    )
    xname<-switch(xvar,
                  norm="|coef|/max|coef|",
                  df="Df",
                  arc.length="Arc Length",
                  step="Step"
    )
    
    if(plottype == "Cp") {
      Cp <- object$Cp
      plot(s1, Cp, type = "b", xlab="|coef|/max|coef|", ...)
    }
    else {
      matplot(s1, coef1, xlab = '', ..., type = "b", lwd = 2,lty=1,
              pch = "*", ylab = "Coefficients")
      abline(h = 0, lty = 3)
      axis(4, at = coef1[nrow(coef1),  ], labels = paste(cnums
      ), cex.axis = 1, adj = 0)
      if(breaks) {
        axis(3, at = s1, labels = paste(stepid),cex.axis=1)
        abline(v = s1)
      }
      
    }
    invisible()
  }

plotCVLars <-
  function(cv.lars.object,se=TRUE){
    mode=cv.lars.object$mode
    xlab=switch(mode,
                fraction="|coef|/max|coef|",
                step="Number of steps"
    )
    index=cv.lars.object$index
    cv=cv.lars.object$cv
    cv.error=cv.lars.object$cv.error
    plot(index, cv, type = "l", ylim = 
           range(cv, cv + cv.error, cv - cv.error),xlab='|coef|/max|coef|',ylab="Cross-Validated MSE")
    if(se)
      error.bars(index, cv + cv.error, cv - cv.error, 
                 width = 1/length(index))
    invisible()
  }

cv.lars2 <-
  function(x, y, K = 10, index, 
           trace = FALSE, plot.it = TRUE, se = TRUE,
           type = c("lasso", "lar", "forward.stagewise", "stepwise"),
           mode=c("fraction", "step"),...)
  {
    type=match.arg(type)
    
    if(missing(mode)){
      mode=switch(type,
                  lasso="fraction",
                  lar="step",
                  forward.stagewise="fraction",
                  stepwise="step"
      )
    }
    else  mode=match.arg(mode)
    all.folds <- cv.folds(length(y), K)
    if(missing(index)){
      index=seq(from = 0, to = 1, length = 100)
      if(mode=="step"){
        fit=lars(x,y,type=type,...)
        nsteps=nrow(fit$beta)
        maxfold=max(sapply(all.folds,length))
        nsteps=min(nsteps,length(y)-maxfold)
        index=seq(nsteps)
      }
    }
    residmat <- matrix(0, length(index), K)
    for(i in seq(K)) {
      omit <- all.folds[[i]]
      fit <- lars(x[ - omit,,drop=FALSE  ], y[ - omit], trace = trace, type=type,...)
      fit <- predict(fit, x[omit,  ,drop=FALSE], mode = mode, s = index
      )$fit
      if(length(omit)==1)fit<-matrix(fit,nrow=1)
      residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
      if(trace)
        cat("\n CV Fold", i, "\n\n")
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object<-list(index = index, cv = cv, cv.error = cv.error,mode=mode)
    if(plot.it) plotCVLars(object,se=se)
    invisible(object)
  }

set.seed(73)
n <- 50                                     
mu <- rep(0,6)                                   
v <- diag(6)
rho = 0.5
v[1,2] = rho
v[2,1] = rho
v[3,4] = rho
v[4,3] = rho
v[5,6] = rho
v[6,5] = rho
betac = 10
X <- centralize(x = mvrnorm(n, mu, Sigma = v))
beta = c(betac,0,betac,0,betac,0)

y = X %*% beta + rnorm(n,0,2)
par(mfrow=c(1,3))
par(mar = c(4, 4, 2, 2) + 0.1)
result = CookDisLasso(X,y,fineness=100,threshold = F)

##fig 1
plot(lars(X,y))

##fig 2
cd = apply(result$CD_Mat,2,mean)
cd.error = apply(result$CD_Mat,2,sd)
plot(result$fraction, cd,type='l',
     xlab = '|coef|/max|coef|',ylab = "Mean Cook's distance",lwd=1)
abline(v = result$fraction[which.min(cd[0:50])],lty=2,col='gray',lwd=3)

##fig 3
cverror = cv.lars2(X,y,K=n, se=F)
abline(v = which.min(cverror$cv)/100 ,lty=2,col='gray',lwd=3)

####################################################################################
#------------- figure 3 case influence graph of 4 different measures. -------------#

case_influence_plotter = function(fraction, CD_Mat,title){
  plot_table = t(rbind(fraction, CD_Mat))
  # threshold_table = cbind(fraction, apply(CD_Mat,2,mean)*3)
  colnames(plot_table) = c('fraction',1:n)
  # colnames(threshold_table)=c('fraction','val')
  df = as_tibble(plot_table)%>%gather('obs','val',-fraction)
  # df2 = as_tibble(threshold_table)%>%mutate(obs='threshold')
  p = ggplot(data=df,aes(x=fraction,y=val,group=obs,color=obs)) +
    geom_line() + xlim(c(0,NA))+ labs(x = "|coef|/max|coef|")+ggtitle(title)+
    theme(plot.title = element_text(vjust = 0,size = 10,hjust = 0.5),legend.position = "none")+
    # geom_line(data=df2,aes(x=fraction, y=val),linetype = "dotted", color = 'red1')+
    labs(y="influence measure")+theme(panel.background = element_rect(fill = "white"),
                                      panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                                      panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                                      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  
}

df<-read.table("datasets/prostate.txt")
X = centralize(as.matrix(df[,1:8]))
y = as.matrix(df[,9])

p = dim(X)[2]
n = dim(X)[1]

obj = lars(X,y,type='lasso')
denom = sum(lm(y~X)$residual**2)/(n-p-1)*(p+1)

# CD exact
cd_exact = CookDisLasso(X,y,fineness = 500, threshold = F)
# CD approx JKSS
lambda_list = cd_exact$Lambda_list
cd_approx = c()
cd_approx2 = c()
cd_local = c()
for (l in lambda_list){
  # regular lasso estimate with lambda given
  beta_hat = predict.lars(obj,mode='lambda', s=l, type = 'coefficients')$coefficients
  index = which(beta_hat!=0)
  # W^-
  beta_hat_inv = 1/abs(beta_hat)
  beta_hat_inv[which(is.infinite(beta_hat_inv))] = 0
  
  resi = y - mean(y) - X%*%beta_hat 
  H = X%*%solve(t(X)%*%X + l*diag(beta_hat_inv))%*%t(X)
  H2 = cbind(1,X[,index])%*%solve(t(cbind(1,X[,index]))%*%cbind(1,X[,index]))%*%t(cbind(1,X[,index]))
  H_vec = diag(H)
  H_vec2 = diag(H2)
  cd_approx = cbind(cd_approx,H_vec/(1-H_vec)**2*resi**2/denom)
  cd_approx2 = cbind(cd_approx2, H_vec2/(1-H_vec2)**2*resi**2/denom)
  cd_local = cbind(cd_local, H_vec2*resi**2/denom)
}

p1 = case_influence_plotter(cd_exact$fraction,cd_exact$CD_Mat,'Exact')
p4 = case_influence_plotter(cd_exact$fraction,cd_approx,'Ridge')
p3 = case_influence_plotter(cd_exact$fraction,cd_approx2,'Approx')
p2 = case_influence_plotter(cd_exact$fraction,cd_local,'Local')

ggarrange(p1+rremove('xlab'), 
          p2+rremove('ylab')+rremove('xlab'), 
          p3, 
          p4+rremove('ylab'), 
          ncol = 2, nrow = 2)

