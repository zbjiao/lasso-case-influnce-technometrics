rm(list=ls())
source("caseweightlasso2024.R")

####################################################################################
#-------------------- figure 5 relationship between xi and omega ------------------#
plot_xi <- function(w, h){
  # given omega and hkk, return xi
  (1-w)/(1-(1-w)*h)
}
w = 0:100/100
lev = c(0.01,0.2,0.5)
df = data.frame(cbind(c(plot_xi(w,0.01), plot_xi(w,0.2), 
                        plot_xi(w,0.5)), rep(w,length(lev)), rep(lev,each = 101)))
colnames(df) = c('xi','w','leverage')
df$leverage = factor(df$leverage)
ggplot(data=df,aes(x=w,y=xi,linetype=leverage))+geom_line()+
  xlab(TeX("Case weight $\\omega$"))+ylab(TeX("$\\xi(\\omega)$")) + 
  guides(linetype=guide_legend(title=TeX("$h_{kk}$")))+
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))


####################################################################################
#----------- figure 6 a solution path demo using synthetic dataset ----------------#
set.seed(100)
x=matrix(rnorm(20*5),nrow=20)
y = x%*%c(10,0,0,0,0)
x[1,] = c(30,10,10,2,1)/5
y[1] = x[1,]%*%c(-10,10,5,-8,0)-2
x = centralize(x)
beta_path(x,y,1,lambda=1.5,plot=2)
rm(x); rm(y)

####################################################################################
#--------------------- figure 7 justify the chisq1 threshold ----------------------#

data <- c()
data(diabetes)
attach(diabetes)

# row 1 col 1
lambda_max = max(abs(t(x)%*%y))
result = CookDisLasso(x, y, s = 200, mode='lambda')
order = 1
data = rbind(data, cbind(order, result$CD_Mat[,1]))
# row 1 col 2
fit = cv.glmnet(x,y,standardize = F)
n = nrow(x)
result = CookDisLasso(x, y, s = fit$lambda.min*n, mode='lambda')
order = 2
data = rbind(data, cbind(order, result$CD_Mat[,1]))
# row 1 col 3
order = 3
data = rbind(data, cbind(order, cooks.distance(lm(y~x))))

detach(diabetes)

dat0=read.csv("datasets/GBM3600Set55and65andClinicalInformation.csv")
# the following data frame contains
# the gene expression data: columns are genes, rows are arrays (samples)
datExprdataOne =data.frame(t(dat0[-c(1:4),18:72]))
datExprdataTwo=data.frame( t(dat0[-c(1:4),73:137]))
names(datExprdataOne)=as.character( dat0$gbm133a[-c(1:4)])
names(datExprdataTwo)=as.character( dat0$gbm133a[-c(1:4)])
# these data frames contain the clinical data of the patients
datClinicaldataOne=dat0[1:4,18:72]
datClinicaldataTwo=dat0[1:4,73:137]
datClinicaldataOne =data.frame(t(dat0[c(1:4),18:72]))
names(datClinicaldataOne)=as.character(dat0[1:4,1])
datClinicaldataTwo=data.frame(t(dat0[c(1:4),73:137]))

#getting the gene names
gnames<-read.csv("datasets/GBM3600Set55and65andClinicalInformation.csv")
gnames<-as.vector(gnames[,6])[-(1:4)]
gnames<-c("intercept",gnames)

x1<-datExprdataOne[datClinicaldataOne[,1]==1,]
y1<-log(datClinicaldataOne[datClinicaldataOne[,1]==1,2])
x<-x1
x<-as.matrix(x)
x<-centralize(t(scale(t(log10(x)))))
y<-as.vector(y1)

# row 2 col 1
lambda_max = max(abs(t(x)%*%y))
result = CookDisLasso(x, y, s = 3, mode = 'l')
order = 4
data = rbind(data, cbind(order, result$CD_Mat[,1]))
# row 2 col 2
fit = cv.glmnet(x,y,standardize = F)
n = nrow(x)
result = CookDisLasso(x, y, s = fit$lambda.min*n, mode='l')
order = 5
data = rbind(data, cbind(order, result$CD_Mat[,1]))
# row 2 col 3
result = CookDisLasso(x, y, s = 0.01, mode='l')
order = 6
data = rbind(data, cbind(order, result$CD_Mat[,1]))


# row 3 col 1
df<-read.table("datasets/prostate.txt")
X = centralize(as.matrix(df[,1:8]))
y = as.matrix(df[,9])
lambda_max = max(abs(t(X)%*%y))
result = CookDisLasso(X, y, s = 8,mode='l')
order = 7
data = rbind(data, cbind(order, result$CD_Mat[,1]))

# row 3 col 2
fit = cv.glmnet(X,y,standardize = F)
n = nrow(X)
result = CookDisLasso(X, y, s = fit$lambda.min*n,mode='l')
order = 8
data = rbind(data, cbind(order, result$CD_Mat[,1]))

# row 3 col 3
order = 9
data = rbind(data, cbind(order, cooks.distance(lm(y~X))))



# row 4 col 1
set.seed(100)
x = matrix(rnorm(200*50),nrow=200)
x = centralize(x)
y = x[,1:5]%*%c(5,4,3,2,1) + rnorm(200)
lambda_max = max(abs(t(x)%*%y))

result = CookDisLasso(x, y, s = 0.2, mode='f')
order = 10
data = rbind(data, cbind(order, result$CD_Mat[,1]))

# row 4 col 2
fit = cv.glmnet(x,y,standardize = F)
n = nrow(x)
result = CookDisLasso(x, y, s = fit$lambda.min*n, mode='lambda')
order = 11
data = rbind(data, cbind(order, result$CD_Mat[,1]))

# row 4 col 3
order = 12
data = rbind(data, cbind(order, cooks.distance(lm(y~x))))

# col_names <- c("|A|=1", "0", "CV")
# Row names
# row_names <- c("Diabetes", "Glioblastoma", "Prostate", "Synethetic")
data = data[,1:2]
data = as.data.frame(data)
colnames(data) = c('order','value')
data[data$order<=3,'source'] = 'Diabetes'
data[data$order>3 & data$order<=6,'source'] = 'Glioblastoma'
data[data$order>6 & data$order<=9,'source'] = 'Prostate'
data[data$order>9 & data$order<=12,'source'] = 'Synthetic'
data[data$order %in% c(1,4,7,10),'pen'] = '|A|=1'
data[data$order %in% c(2,5,8,11),'pen'] = 'no penalty'
data[data$order %in% c(3,6,9,12),'pen'] = 'cross-validated penalty'
data$source = as.factor(data$source)
data$pen = as.factor(data$pen)

data%>%group_by(order)%>%mutate(value=value/sqrt(var(value)/2), value = ifelse(value<=9, value, NA))%>%
  ggplot( aes(x=value, color=source, fill=source)) +
  geom_histogram(alpha=0.6,aes(y=after_stat(density)),bins = 25,na.rm = T) +
  stat_function(fun = dchisq, args = list(df = 1),color='black')+ylim(c(0,1))+xlim(c(0,9))+
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.15, "lines"),
    axis.text.y = element_text(size=9),
    strip.text.x = element_text(size = 10, face = "bold.italic",hjust=0.5),
    strip.text.y = element_text(size = 10, face = "bold.italic",hjust=0.5)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )+
  xlab("") + ylab("Density")+
  facet_grid(source ~ pen) 

rm(x); rm(y)


####################################################################################
#---- figure 8 probability density estimates of studentized lasso resids  ---------#
get_std_residual = function(x,y,name){
  
  m = glmnet(x,y)
  selection = as.matrix(m$beta)>0
  n = dim(x)[1]
  h_summary = c()
  for (i in 1:dim(selection)[2]){
    if (sum(selection[,i])==0){
      h = rep(1/n, n)
    }
    else{
      x_ = cbind(1,x[,selection[,i]])
      h = diag(x_%*%solve(t(x_)%*%x_)%*%t(x_))
    }
    h_summary = cbind(h_summary, h)
  }
  
  
  residuals_df <- as.data.frame(t(apply(predict(m,x),2,function(x) x-y)/sqrt(1-h_summary)))
  # apply(m$beta,2,function(x) sum(abs(x)))/sum(abs(lm(y~x)$coefficients)[2:dim(x)[2]])
  # Add a column identifier
  residuals_df$fraction <- apply(m$beta,2,function(x) sum(abs(x)))/sum(abs(glmnet(x,y,lambda = 0)$beta))
  residuals_df$fraction = (residuals_df$fraction - min(residuals_df$fraction))/
    (max(residuals_df$fraction)- min(residuals_df$fraction))
  # sum(abs(lm(y~x)$coefficients)[2:dim(x)[2]])
  
  # Pivot longer to get a suitable format for ggplot
  residuals_long <- pivot_longer(residuals_df, -fraction, names_to = "Variable", values_to = "Residual")
  residuals_long = residuals_long %>% group_by(fraction) %>% mutate(Residual = Residual/sd(Residual))
  residuals_long$name = name
  return(residuals_long)
}

data(diabetes)
attach(diabetes)
residuals_long1 = get_std_residual(x,y,'Diabetes')
detach(diabetes)

dat0=read.csv("datasets/GBM3600Set55and65andClinicalInformation.csv")
# the following data frame contains
# the gene expression data: columns are genes, rows are arrays (samples)
datExprdataOne =data.frame(t(dat0[-c(1:4),18:72]))
datExprdataTwo=data.frame( t(dat0[-c(1:4),73:137]))
names(datExprdataOne)=as.character( dat0$gbm133a[-c(1:4)])
names(datExprdataTwo)=as.character( dat0$gbm133a[-c(1:4)])
# these data frames contain the clinical data of the patients
datClinicaldataOne=dat0[1:4,18:72]
datClinicaldataTwo=dat0[1:4,73:137]
datClinicaldataOne =data.frame(t(dat0[c(1:4),18:72]))
names(datClinicaldataOne)=as.character(dat0[1:4,1])
datClinicaldataTwo=data.frame(t(dat0[c(1:4),73:137]))

#getting the gene names
gnames<-read.csv("datasets/GBM3600Set55and65andClinicalInformation.csv")
gnames<-as.vector(gnames[,6])[-(1:4)]
gnames<-c("intercept",gnames)

x1<-datExprdataOne[datClinicaldataOne[,1]==1,]
y1<-log(datClinicaldataOne[datClinicaldataOne[,1]==1,2])
x<-x1
x<-as.matrix(x)
x<-centralize(t(scale(t(log10(x)))))
y<-as.vector(y1)
residuals_long2 = get_std_residual(x,y,'Glioblastoma')

df<-read.table("datasets/prostate.txt")
x = centralize(as.matrix(df[,1:8]))
y = as.matrix(df[,9])
residuals_long3 = get_std_residual(x,y,'Prostate')

x <- matrix(rnorm(200*500), nrow=500)
y <- x[,1:5] %*% c(5,4,3,2,1) + rnorm(500,0,4)
residuals_long4 = get_std_residual(x,y,'Synthetic')

residuals_long = rbind(residuals_long1,residuals_long2,residuals_long3,residuals_long4)


residuals_long$name = as.factor(residuals_long$name)

# plot

p = ggplot(residuals_long, aes(x = Residual, group=fraction,color = fraction)) +
  geom_density() + # Adjust alpha for fading effect
  scale_color_viridis_c(option = 'G',alpha = 0.5) + # Use a color scale that supports fading well
  theme_minimal() +
  stat_function(fun = dnorm,color='red',linetype="dashed",group=0,linewidth=1,alpha=.8)+
  labs(
    x = "Studentized Residuals",
    y = "Density") +
  theme(strip.text.x = element_text(size = 10))+
  guides(color = 'none') + facet_wrap(~name)
p

rm(x); rm(y)


####################################################################################
#------------------------ figure 9 model selection example ------------------------#
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
par(mfrow=c(2,2))
par(mar = c(4, 4, 2, 2) + 0.1)
result = CookDisLasso(X,y,fineness=100,threshold = F)

##fig 5-1
plot(lars(X,y))

##fig 5-2
cd = apply(result$CD_Mat,2,mean)
cd.error = apply(result$CD_Mat,2,sd)
plot(result$fraction, cd,type='l',
     xlab = '',ylab = "Mean Cook's distance",lwd=1)
abline(v = result$fraction[which.min(cd[0:50])],lty=2,col='gray',lwd=3)

##fig 5-3
matplot(result$fraction, t(result$CD_Mat),xlab = "|coef|/max|coef|",
        ylab = "Cook's distance under Lasso", lty = 1,type = "l")

##fig 5-4
cverror = cv.lars2(X,y,K=n, se=F)
abline(v = which.min(cverror$cv)/100 ,lty=2,col='gray',lwd=3)

####################################################################################
#--------------- figure 10 time complexity analysis simultion  ----------------#

npl_simu = read.csv('results/time_complexity_simualtion.csv')

npl_plotter = npl_simu %>% mutate(ratio = as.numeric(t2)/as.numeric(t1))

# 1
p121 = npl_plotter %>% filter(p<1,l==1) %>% ggplot(aes(x=n, y=ratio,group=p)) + 
  geom_line(aes(color=as.factor(p)))+  geom_point(aes(color=as.factor(p)))+
  theme_minimal() + labs(
    # title = "L = O(1), c<1",
    x = "n (with L = 1)",
    y = "Ratio of time consumed",
    color = "c"
  )+theme(plot.title = element_text(hjust = 0.5))
# +annotate("text", x = 40, y = 0.3, label = "L=1")

# 2
p122 = npl_plotter %>% filter(l==1, p>1) %>% ggplot(aes(x=p, y=ratio,group=n)) + 
  geom_line(aes(color=factor(n)))+  geom_point(aes(color=factor(n)))+
  theme_minimal() + labs(
    # title = "L = O(1), c>1",
    # x = TeX("c where $p = {n^c}$"),
    x = "c (with L = 1)",
    y = "Ratio of time consumed",
    color = "n"
  )+theme(plot.title = element_text(hjust = 0.5))
# +annotate("text", x = 1.1, y = 0.7, label = "L=1")

# 3
p123 = npl_plotter %>% filter(((p==0.7) & (l==3))|((p==0.8) & (l==4))|((p==0.9) & (l==5)))%>% 
  ggplot(aes(x=p, y=ratio,group=n)) + 
  geom_line(aes(color=factor(n)))+  geom_point(aes(color=factor(n)))+
  theme_minimal() + labs(
    # title = "c = 0.7",
    # x = "L",
    x = "c (with L = O(p))",
    y = "Ratio of time consumed",
    color = "n"
  )+theme(plot.title = element_text(hjust = 0.5))
# +annotate("text", x = 1.5, y = 0.8, label = "c=0.7")

# 4
p124 = npl_plotter %>% filter(n/l==10, p>1) %>% ggplot(aes(x=p, y=ratio,group=n)) + 
  geom_line(aes(color=factor(n)))+  geom_point(aes(color=factor(n)))+
  theme_minimal() + labs(
    # title = "L = O(n), c>1",
    # x = TeX("c where $p = {n^c}$"),
    x = "c (with L = n/10)",
    y = "Ratio of time consumed",
    color = "n"
  )+theme(plot.title = element_text(hjust = 0.5))
# +annotate("text", x = 1.3, y = 3, label = "L = n/10")

ggarrange(p121,
          p122+rremove('ylab'), 
          p123, 
          p124+rremove('ylab'), 
          ncol = 2, nrow = 2)

####################################################################################
#--------------------- figure 11 needle plot of 4 measures ------------------------#
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

new_slice_pot = function(a){
  df = cbind(1:dim(cd_approx)[1], cd_exact$CD_Mat[,a], cd_approx[,a], cd_approx2[,a], cd_local[,a])
  colnames(df) = c('index','Exact','JKSS','Approx','Local')
  beta_hat = predict.lars(obj,mode='lambda', s=lambda_list[a], type = 'coefficients')$coefficients
  # df = data.frame(df) %>% gather('method', 'cd',-index)
  df = data.frame(df)
  q1 = ggplot(data=df,aes(x=index, ymax=Exact,ymin=0)) + geom_linerange()+
    ggtitle('Exact')+theme(plot.title = element_text(hjust = 0.5,size=10))+
    ylim(0,0.09)+
    geom_segment(aes(x = 3 , y = 0, xend = 3, yend = df[3,2]),color='red',linewidth=1.5)+
    # geom_segment(aes(x = 58 , y = 0, xend = 58, yend = df[58,2]),color='red',linewidth=1.5)+
    geom_segment(aes(x = 69 , y = 0, xend = 69, yend = df[69,2]),color='blue',linewidth=1.5)+
    geom_segment(aes(x = 95 , y = 0, xend = 95, yend = df[95,2]),color='red',linewidth=1.5)+
    geom_segment(aes(x = 96 , y = 0, xend = 96, yend = df[96,2]),color='blue',linewidth=1.5)+
    geom_segment(aes(x = 39 , y = 0, xend = 39, yend = df[39,2]),color='blue',linewidth=1.5)
    # annotate("text", x = 25, y = max(df$Exact)*0.9, label = "Exact")+

  q2 = ggplot(data=df,aes(x=index, ymax=JKSS,ymin=0)) + geom_linerange()+
    ggtitle('Ridge')+theme(plot.title = element_text(hjust = 0.5,size=10))+
    ylim(0,0.09)+
    geom_segment(aes(x = 69 , y = 0, xend = 69, yend = df[69,3]),color='blue',linewidth=1.5)+
    geom_segment(aes(x = 96 , y = 0, xend = 96, yend = df[96,3]),color='blue',linewidth=1.5)+
    geom_segment(aes(x = 39 , y = 0, xend = 39, yend = df[39,3]),color='blue',linewidth=1.5)
    
  q3 = ggplot(data=df,aes(x=index, ymax=Approx,ymin=0)) + geom_linerange()+
    ggtitle('Approx')+theme(plot.title = element_text(hjust = 0.5,size=10))+
    ylim(0,0.09)+
    geom_segment(aes(x = 3 , y = 0, xend = 3, yend = df[3,4]),color='red',linewidth=1.5)+
    geom_segment(aes(x = 95 , y = 0, xend = 95, yend = df[95,4]),color='red',linewidth=1.5)
    
  q4 = ggplot(data=df,aes(x=index, ymax=Local,ymin=0)) + geom_linerange()+
    ggtitle('Local')+theme(plot.title = element_text(hjust = 0.5,size=10))+
    ylim(0,0.09)
    
  qq = ggarrange(q1+rremove('xlab')+rremove("x.grid")+rremove("x.text")+rremove('ylab')+
                   theme(panel.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                         panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)),
                 q3+rremove('xlab')+rremove("x.grid")+rremove("x.text")+rremove('ylab')+
                   theme(panel.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                         panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)),
                 q4+rremove('xlab')+rremove("x.grid")+rremove("x.text")+rremove('ylab')+
                   theme(panel.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                         panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)),
                 q2+rremove('xlab')+rremove("x.grid")+rremove("x.text")+rremove('ylab')+
                   theme(panel.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "grey",linewidth = 0.5),
                         panel.grid.minor = element_line(colour = "white",linewidth = 0.5),
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)),
                 ncol=1,nrow=4)
}
print(new_slice_pot(197))


####################################################################################
#--------------- figure 12 a solution path for diabetes data  ----------------#

if (exists("y")) rm(y)

par(mfrow=c(1,1))
data(diabetes)
attach(diabetes)

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
      axis(4, at = coef1[nrow(coef1),  ],
           labels = paste(cnums), cex.axis = 1, 
           las = 0,line=0)
      if(breaks) {
        axis(3, at = s1, labels = paste(stepid),cex.axis=1)
        abline(v = s1)
      }
      
    }
    invisible()
  }
colname=colnames(x)
plot(lars(x,y))
detach(diabetes)

####################################################################################
#--------------- figure 13 residual/leverage plot of diabetes data ----------------#
data(diabetes)
attach(diabetes)

fit = cv.glmnet(x,y,standardize = F)
est = glmnet(x,y,lambda=fit$lambda.min,standardize=FALSE,thresh=1e-16)
n = dim(x)[1]
p = dim(x)[2]
denom = sum(lm(y~x)$residual**2)/(n-p-1)
result = CookDisLasso(x,y,s=3,mode='l')
residy = y - mean(y) - (x%*%result$beta)
x_ = x[,result$beta!=0]
lev = diag(x_%*%solve(crossprod(x_))%*%t(x_))
df_ = data.frame(index = 1:442, lev = lev, res = residy/sqrt((1-lev)*denom), Case_Influence = result$CD_Mat)

ggplot(data=df_, aes(x=lev, y=res, size=Case_Influence)) +
  geom_point(aes(color = Case_Influence > 0.025)) +  
  scale_color_manual(values = c("TRUE" = "red", "F" = "black")) +
  xlab(TeX('leverage')) + 
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  ylab(TeX('studentized residual')) +
  theme(panel.background = element_rect(fill = "white"),
        legend.position = 'None',
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) 

detach(diabetes)
