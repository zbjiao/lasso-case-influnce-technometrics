rm(list=ls())
source("caseweightlasso2024.R")

set.seed(0)

# dataset
concrete = read_excel("datasets/Concrete_Data.xls")
colnames(concrete) = c('cement','bfs','fash','water','sp','ca','fagg','age','strength')
age = concrete$age
concrete = concrete %>% mutate(
  recip_water = 1/water,
  recip_cement = 1/cement
) %>% 
  dplyr::select(-age)
full_y = concrete$strength
concrete = concrete[,-8]

original_vars = names(concrete)

interaction_formula = as.formula(paste("~ (", paste(original_vars, collapse = " + "), ")^2"))
interaction_matrix = model.matrix(interaction_formula, data = concrete)
# remove intercept column
interaction_matrix = interaction_matrix[,2:dim(interaction_matrix)[2]]

interaction_df = as.data.frame(interaction_matrix)

# Remove interaction terms between covariates and their own reciprocals
interaction_cols_to_remove = sapply(names(interaction_df), function(colname) {
  # Split the column name by ':' to identify interaction terms
  terms = strsplit(colname, ":")[[1]]
  # Check if the interaction involves a covariate and its own reciprocal
  if (length(terms) == 2 && (paste0("recip_", terms[1]) == terms[2] || paste0("recip_", terms[2]) == terms[1])) {
    return(TRUE)  # Mark this column for removal
  }
  return(FALSE)
})

# Remove the unwanted interaction columns
interaction_df = interaction_df[, !interaction_cols_to_remove]
growing_age = ifelse(age>28, 28, age)
interaction_df2 = cbind(interaction_df, age, growing_age)

result_before = c()
result_after = c()
fraction_before = c()
fraction_after = c()
r_sq_before = c()
r_sq_after = c()

noe = 1000

finalized_df = cbind(log(1+interaction_df2), full_y)
colnames(finalized_df)[length(colnames(finalized_df))] = 'y'
test_index = sample(dim(finalized_df)[1], dim(finalized_df)[1]*0.25)
train_df = finalized_df[-test_index, ]
test_df = finalized_df[test_index,]

x = as.matrix(train_df[1:(dim(train_df)[2]-1)])
x_pack = centralize(x,r=T)
x = x_pack$d
y = train_df$y

for (i in 1:noe){
  fit = cv.lars(x,y,plot.it = F)
  obj = lars(x,y)
  frac = fit$index[which.min(fit$cv)]
  
  est = as.vector(predict.lars(obj,mode='fraction', s=frac, type = 'coefficients')$coefficients)

  result = CookDisLasso(x,y,s=c(min(frac,0.99)),mode = 'fraction')
  noninfluential_indices = which(result$CD_Mat < result$threshold)
  
  x_trimmed = x[noninfluential_indices,]
  y_trimmed = y[noninfluential_indices]
  
  x_trimmed_pack = centralize(x_trimmed,r=T)
  x_trimmed = x_trimmed_pack$d
  fit2 = cv.lars(x_trimmed,y_trimmed,plot.it = F)
  obj2 = lars(x_trimmed,y_trimmed)
  frac2 = fit2$index[which(fit2$cv == min(fit2$cv))]
  est2 = as.vector(predict.lars(obj2,mode='fraction', s=frac2, type = 'coefficients')$coefficients)
  
  if (i == 1){
    # grab the first model fitted as an example for demonstration in figure 9
    est_example_1 = est
    est_example_2 = est2
    m2 = mean(y_trimmed)
    x2 = x_trimmed_pack
  }
  
  x_test_1 = scale(test_df[,1:(dim(test_df)[2]-1)], x_pack$m, x_pack$v)
  x_test_2 = scale(x_test_1, x_trimmed_pack$m, x_trimmed_pack$v)
  
  r_sq_before = c(r_sq_before, 1-sum((test_df$y - (x_test_1 %*% est) - mean(y))**2)/length(test_df$y)/var(test_df$y))
  r_sq_after = c(r_sq_after, 1-sum((test_df$y - (x_test_2 %*% est2) - mean(y_trimmed))**2)/length(test_df$y)/var(test_df$y))
  
  result_before = cbind(result_before,est)
  result_after = cbind(result_after,est2)
  fraction_before = c(fraction_before,frac)
  fraction_after = c(fraction_after,frac2)
  if (i%%20 == 0) print(paste0(i,'/',1000))
  # print(c(frac, frac2, r_sq_before[i], r_sq_after[i]))
}



rda_ans = data.frame(mean1 = apply(result_before,1,mean),
                     mean2 = apply(result_after,1,mean),
                     sd1 = apply(result_before,1,sd),
                     sd2 = apply(result_after,1,sd),
                     nonzero_rate1 = apply(result_before, 1,function(x) mean(x!=0)),
                     nonzero_rate2 = apply(result_after, 1,function(x) mean(x!=0)),
                     nonzero_sd1 = apply(result_before, 1,function(x) sd(x!=0)),
                     nonzero_sd2 = apply(result_after, 1,function(x) sd(x!=0)))
rownames(rda_ans) = colnames(x)

sum(rda_ans$nonzero_sd2**2)/sum(rda_ans$nonzero_sd1**2)
# [1] 0.7864243
var(fraction_after)/var(fraction_before)
# [1] 0.6175607
mean(r_sq_after)
sd(r_sq_after)
mean(r_sq_before)
sd(r_sq_before)

# > mean(r_sq_after)
# [1] 0.8805996
# > sd(r_sq_after)
# [1] 0.00168687
# > mean(r_sq_before)
# [1] 0.8740616
# > sd(r_sq_before)
# [1] 0.001230681

# write.csv(round(rda_ans,2), "concrete_example.csv")


####################################### Figure 4 #######################################



process_new_observation = function(new_observation, original_data, reciprocal_vars, interaction_formula) {
  # Step 1: Add reciprocal variables to the new observation
  for (colname in names(new_observation)) {
    # If the value is non-zero, calculate its reciprocal
    if (new_observation[[colname]] != 0) {
      new_observation[[paste0("recip_", colname)]] = 1 / new_observation[[colname]]
    }
  }
  
  # Step 2: Create interaction terms for the new observation
  # Use the interaction formula created earlier for two-way interactions
  interaction_matrix_new = model.matrix(interaction_formula, data = as.data.frame(new_observation))
  
  # Remove the intercept column
  interaction_matrix_new = interaction_matrix_new[, 2:dim(interaction_matrix_new)[2]]
  
  
  # Convert to data frame
  interaction_df_new = t(as.data.frame(interaction_matrix_new))
  
  # Remove unwanted interaction columns
  interaction_df_new = log(1+interaction_df_new[,colnames(interaction_df)])
  
  return(interaction_df_new)
}


range = 0:100
df_wcr1 = (0.3+range/200)
df_wcr2 = (0.3+range/200)
age_range = c(1,3,7,14,28,60)
for (age in age_range){
  predicted_y1 = c()
  predicted_y2 = c()
  for (i in range){
    new_observation = data.frame(
      cement = 350, bfs = 10, fash = 10, water = 350*(0.3+i/200), sp = 5, ca = 1100, fagg = 750, age = age) %>%
      mutate(growing_age = ifelse(age<=28, age, 28),
             mature_age = ifelse(age>28, age-28, 0),
             agg_ratio = fagg/(ca+fagg)) 
    processed_new_observation = process_new_observation(new_observation, concrete, reciprocal_vars, interaction_formula)
    processed_new_observation['age'] = log(1+age)
    processed_new_observation['growing_age'] = log(1+ifelse(age>28, age, 28))
    x1 = (processed_new_observation - x_pack$m) / x_pack$v
    predicted_y1 = c(predicted_y1, x1 %*% 
                       as.vector(est_example_1) + mean(y))
    predicted_y2 = c(predicted_y2, ((x1- x2$m)/x2$v) %*% 
                       as.vector(est_example_2) + m2)
  }
  df_wcr1 = cbind(df_wcr1, predicted_y1)
  df_wcr2 = cbind(df_wcr2, predicted_y2)
}

colnames(df_wcr1) = c('ratio',age_range)
colnames(df_wcr2) = c('ratio',age_range)
df_wcr1 = data.frame(df_wcr1)
df_wcr2 = data.frame(df_wcr2)


df1_long = melt(df_wcr1, id.vars = "ratio", variable.name = "age", value.name = "value")
df1_long$approach = "Approach1"  

df2_long = melt(df_wcr2, id.vars = "ratio", variable.name = "age", value.name = "value")
df2_long$approach = "Approach2"  


df_combined = rbind(df1_long, df2_long)

# Rename age levels for proper legend labels
df_combined$age = factor(df_combined$age, levels = c("X1", "X3", "X7", "X14", "X28","X60"),
                          labels = c("1 day", "3 days", "7 days", "14 days", "28 days","60 days"))


p = ggplot(df_combined, aes(x = ratio, y = value, color = age, linetype = approach)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scales::seq_gradient_pal("black", "white")(seq(0, 0.8, length.out = 6))) + 
  labs( x = "Water to Cement ratio", y = "Compressive Strength (MPa)", color = "Age", linetype = "Approach") +
  theme_minimal()

p




