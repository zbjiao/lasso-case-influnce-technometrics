rm(list=ls())
source("caseweightlasso2024.R")

data(diabetes)
attach(diabetes)

################################## Table 3 ##################################
set.seed(1)

result_before = c()
result_after = c()
lambda_before = c()
lambda_after = c()
x = centralize(x)

noe = 1000

for (i in 1:noe){
  fit = cv.glmnet(x,y,standardize = F)
  est = glmnet(x,y,lambda=fit$lambda.min,standardize=FALSE,thresh=1e-16)
  
  result = CookDisLasso(x,y,s=fit$lambda.min*nrow(x), mode = 'lambda')
  noninfluential_indices = which(result$CD_Mat<=result$threshold)
  x_trimmed = x[noninfluential_indices,]
  y_trimmed = y[noninfluential_indices]
  
  fit2 = cv.glmnet(x_trimmed,y_trimmed,standardize = F)
  est2 = glmnet(x_trimmed,y_trimmed,lambda=fit2$lambda.min,standardize=FALSE,thresh=1e-16)
  
  result_before = cbind(result_before,as.vector(coef(est)))
  result_after = cbind(result_after,as.vector(coef(est2)))
  lambda_before = c(lambda_before,fit$lambda.min*nrow(x))
  lambda_after = c(lambda_after,fit2$lambda.min*nrow(x_trimmed))
  
  if (i%%20==0) print(paste0(i,'/',noe))
}


rda_ans = data.frame(mean1 = apply(result_before,1,mean),
                     mean2 = apply(result_after,1,mean),
                     sd1 = apply(result_before,1,sd),
                     sd2 = apply(result_after,1,sd),
                     nonzero_rate1 = apply(result_before, 1,function(x) mean(x!=0)),
                     nonzero_rate2 = apply(result_after, 1,function(x) mean(x!=0)))

rda_ans = rda_ans[2:11,]
rownames(rda_ans) <- colnames(x)

lambda_row = c(mean(lambda_before), mean(lambda_after), sd(lambda_before), sd(lambda_after), NA, NA)
a_row = c(sum(rda_ans$nonzero_rate1), sum(rda_ans$nonzero_rate2), 
          sd(apply(result_before!=0,2,sum)), sd(apply(result_after!=0,2,sum)), NA, NA)

Table3 = rbind(rda_ans, lambda_row, a_row)
rownames(Table3)[11:12] = c('lambda','active_set')

write.csv(round(Table3,4), "results/Table3.csv")

# > mean(lambda_before)
# [1] 14.35132
# > mean(lambda_after)
# [1] 16.35667

# > sd(lambda_before)
# [1] 10.32509
# > sd(lambda_after)
# [1] 4.128719

# > sum(rda_ans$nonzero_rate1)
# [1] 7.983
# > sum(rda_ans$nonzero_rate2)
# [1] 7.666

# > sd(apply(result_before!=0,2,sum))
# [1] 1.084383
# > sd(apply(result_after!=0,2,sum))
# [1] 0.7422904

################################## mse comparison ##################################

set.seed(1)

get_mse <- function(est,x_train,y_train,x_test,y_test){
  beta0 = est$a0
  beta_lasso = as.vector(est$beta)
  active_indices = which(beta_lasso!=0)
  x_train_active = x_train[,active_indices]
  beta_ols = solve(crossprod(x_train_active))%*%t(x_train_active)%*%y_train
  x_test_active = x_test[,active_indices]
  mse_ols = sum((x_test_active %*% beta_ols + beta0 - y_test)**2)
  mse_lasso = sum((x_test %*% beta_lasso + beta0 - y_test)**2)
  return(c(mse_ols, mse_lasso))
}

ans_sum=c()
noe = 1000
for (iter in 1:noe){
  
  ## total number of samples 442
  test_indices = sample(442,440*0.2)
  # train_indices = setdiff(1:442, test_indices)
  # test_indices = 391:442
  # train_indices = 1:390
  
  x_train = x[-test_indices,]
  y_train = y[-test_indices]
  x_test = x[test_indices,]
  y_test = y[test_indices]
  
  ############# regular procedure #################
  central_pack = centralize(x_train,r=T)
  x_train_c = central_pack$d
  x_test_c =  scale(x_test, central_pack$m, central_pack$v)
  
  fit = cv.glmnet(x_train_c,y_train,standardize = F)
  est = glmnet(x_train_c,y_train,lambda=fit$lambda.min,standardize=FALSE,thresh=1e-16)
  
  ans1 = get_mse(est, x_train_c,y_train,x_test_c,y_test)
  
  ############# case influence procedure #################
  
  result = CookDisLasso(x_train_c, y_train, s = fit$lambda.min*nrow(x_train), mode='lambda')
  noninfluential_indices = which(result$CD_Mat<=result$threshold)
  x_train_trimmed = x_train[noninfluential_indices,]
  y_train_trimmed = y_train[noninfluential_indices]
  
  central_pack = centralize(x_train_trimmed,r=T)
  x_train_trimmed_c = central_pack$d
  x_test_c =  scale(x_test, central_pack$m, central_pack$v)
  
  
  fit2 = cv.glmnet(x_train_trimmed_c,y_train_trimmed,standardize = F)
  est2 = glmnet(x_train_trimmed_c,y_train_trimmed,lambda=fit2$lambda.min,standardize=FALSE,thresh=1e-16)
  
  ans2 = get_mse(est2, x_train_trimmed_c,y_train_trimmed,x_test_c,y_test)
  
  ans_sum = rbind(ans_sum, c(ans1,ans2))
  
  if (iter%%20==0) print(paste0(iter,'/',noe))
}

mean(ans_sum[,2]/ans_sum[,4])
mean(ans_sum[,1]/ans_sum[,3])
sd(ans_sum[,2]/ans_sum[,4])
sd(ans_sum[,1]/ans_sum[,3])

# > mean(ans_sum[,2]/ans_sum[,4])
# [1] 0.9989423
# > mean(ans_sum[,1]/ans_sum[,3])
# [1] 0.994879
# > sd(ans_sum[,2]/ans_sum[,4])
# [1] 0.02200895
# > sd(ans_sum[,1]/ans_sum[,3])
# [1] 0.02296662
