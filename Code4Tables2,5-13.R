rm(list=ls())
source("caseweightlasso2024.R")


set.seed(100)

cutoff_simu <- function(q, ev, iter = 1000, 
                        beta = cbind(c(1,1,1), c(1,-1,1)), 
                        a_v = c(0,2,3,5), b_v=c(0,2,3,5), 
                        r=0.2,
                        np_v = matrix(c(50,10, 50,100, 50,50, 100,10), nrow=2)
){
  
  results = c()
  total_configuration = 3 * (length(a_v)-1) * ncol(beta) * ncol(np_v)
  num_config = 0
  
  for (a in a_v){
    for (b in b_v){
      if ((a == 0) & (b == 0)){
        next
      }
      if ((a > 0) & (b > 0) & (a != b)){
        next
      }
      for (np in 1:dim(np_v)[2]){
        for (coef in 1:dim(beta)[2]){
          
          n = np_v[1, np]
          p = np_v[2, np]
          
          cur_beta = c(beta[, coef], rep(0, p-length(beta[, coef])))
          
          if (r!=0){
            cormat<-matrix(0,p,p)
            for(i in 1:p) {
              for(j in 1:p) {
                cormat[i,j] = r^(abs(i-j))
              }
            }
          }
          else{
            cormat = diag(p)
          }
          snr = beta[, coef] %*% cormat[1:length(beta[, coef]), 1:length(beta[, coef])] %*% beta[, coef]/ev**2
          
          
          for(i in 1:iter){

            result <- try(abpr_simu(a, b, np_v[1, np], np_v[2, np], q, beta[, coef], ev, cormat), silent = TRUE)

            if (inherits(result, "try-error")) {
              next
            }
            results = rbind(results, c(snr, a, b, np_v[1, np], np_v[2, np], coef, i, result))
          }
          num_config = num_config + 1
          print(c('one round', a, b, np_v[1, np], np_v[2, np], coef))
          print(paste0(num_config, ' out of ', total_configuration, ' configurations completed. '))
        }
      }
    }
  }
  result_df = data.frame(results)
  colnames(result_df) = c('snr','a', 'b', 'n', 'p', 'beta', 'iter', 'inclusion', 'tf0', 'tf1','avg_upd')
  summary_df = result_df %>% 
    group_by(a, b, n, p, beta,snr) %>% 
    summarise(inclusion_rate = mean(inclusion, na.rm = TRUE),
              outlier = mean(tf0, na.rm = TRUE), 
              others = mean(tf1, na.rm = TRUE),
              avg_upd = mean(avg_upd, na.rm = TRUE)) %>%
    ungroup()
  summary_df
}



abpr_simu <- function(a,b,n,p,q,beta,ev,cormat){
  x = mvrnorm(n=n,rep(0,p),cormat)
  ep = rnorm(n,0,ev)
  ep[1] = 0
  y = x[,1:length(beta)] %*% beta + ep
  x[1,q] = a
  
  y[1] = x[1,1:length(beta)] %*% beta + b
  x = centralize(x)
  fit = cv.glmnet(x,y,standardize = F)

  result = CookDisLasso(x,y,s=fit$lambda.min*n, mode = 'lambda')

  inclusion = result$beta_table[q,1]!=0
  
  tf0 = result$CD_Mat[1,1]>result$threshold[1,1]
  tf1 = sum(result$CD_Mat[2:n,1]>result$threshold[2:n,1])/(n-1)
  avg_upd = result$avg_update
  return(c(inclusion, tf0, tf1, avg_upd))
}


real_simu_result1 = cutoff_simu(1,1,1000,beta = cbind(c(1,1,1/2, 1/2), c(1,-1,1/2,-1/2), c(1,1,1,0), c(4,3,2,1)),
                                a_v = c(0,2,3,5), b_v=c(0,2,3,5), 
                                r=0.2,
                                np_v = matrix(c(50,10, 50,500, 200,200, 500,10), nrow=2))

real_simu_result2 = cutoff_simu(10,1,1000,beta = cbind(c(1,1,1/2, 1/2), c(1,-1,1/2,-1/2), c(1,1,1,0), c(4,3,2,1)),
                                a_v = c(0,2,3,5), b_v=c(0,2,3,5), 
                                r=0.2,
                                np_v = matrix(c(50,10, 50,500, 200,200, 500,10), nrow=2))

# write.csv(real_simu_result1 %>% arrange(beta, a,b), file = "results/simu_result1_q1.csv")
# write.csv(real_simu_result2 %>% arrange(beta, a,b), file = "results/simu_result1_q10.csv")

real_simu_result1 = round(real_simu_result1 %>% dplyr::select(a,b,n,p, beta, outlier, others),2)
real_simu_result2 = round(real_simu_result2 %>% dplyr::select(a,b,n,p, beta, outlier, others),2)


# Table 2
Table2a = real_simu_result2 %>% filter(beta == 2, a==b) %>% dplyr::select(-beta)
Table2b = real_simu_result1 %>% filter(beta == 2, a==b) %>% dplyr::select(-beta) 

# Table 5
Table5 = real_simu_result1 %>% filter(beta == 1) %>% dplyr::select(-beta)

# Table 6
Table6 = real_simu_result1 %>% filter(beta == 2) %>% dplyr::select(-beta)

# Table 7 
Table7 = real_simu_result1 %>% filter(beta == 3) %>% dplyr::select(-beta)

# Table 8
Table8 = real_simu_result1 %>% filter(beta == 4) %>% dplyr::select(-beta)

# Table 9
Table9 = real_simu_result2 %>% filter(beta == 1) %>% dplyr::select(-beta)

# Table 10
Table10 = real_simu_result2 %>% filter(beta == 2) %>% dplyr::select(-beta)

# Table 11
Table11 = real_simu_result2 %>% filter(beta == 3) %>% dplyr::select(-beta)

# Table 12
Table12 = real_simu_result2 %>% filter(beta == 4) %>% dplyr::select(-beta)

# Table 13
comparison_simu_result = cutoff_simu(100,1,1000,beta = cbind(c(1,2,3,4,5)),
                                     a_v = c(0,2,5,10), b_v=c(0,2,5,10), 
                                     r=0.5,
                                     np_v = matrix(c(50,1000), nrow=2))

Table13 = comparison_simu_result %>%  dplyr::select(a,b,outlier) %>% 
  rename('1' = outlier) %>%  
  # Rajaratnam et al's result:
  mutate('1df' = c(0.09,0.96,1.00,0.00,0.11,0.00,0.96,0.00,1.00)) %>% 
  arrange(-pmax(a,b),-a, -(a+b))

# # A tibble: 9 × 4
# a     b   `1` `1df`
# <dbl> <dbl> <dbl> <dbl>
#   1    10    10 0.999  1   
# 2    10     0 0.019  0   
# 3     0    10 0.998  1   
# 4     5     5 0.93   0.96
# 5     5     0 0.013  0   
# 6     0     5 0.905  0.96
# 7     2     2 0.206  0.11
# 8     2     0 0.016  0   
# 9     0     2 0.212  0.09

# save all generated tables to the results folder. 
for (i in c('2a','2b',5,6,7,8,9,10,11,12,13)) {
  df_name <- paste0("Table", i)
  
  if (exists(df_name)) {
    file_name <- paste0('results/', df_name, ".csv")
    
    write.csv(get(df_name), file = file_name)
  }
}


