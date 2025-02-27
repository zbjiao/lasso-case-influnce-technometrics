rm(list=ls())
source("caseweightlasso2024.R")

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

set.seed(1)

result_before = c()
result_after = c()
lambda_before = c()
lambda_after = c()
x = centralize(x)

# number of experiment 
noe = 100

for (i in 1:noe){
  fit = cv.lars(x,y,plot.it = F, use.Gram = F)
  obj = lars(x,y,use.Gram = F)
  frac = fit$index[which(fit$cv == min(fit$cv))]
  est = as.vector(predict.lars(obj,mode='fraction', s=frac, type = 'coefficients')$coefficients)
  
  result = CookDisLasso(x,y,s=c(frac),mode = 'fraction')
  noninfluential_indices = which(result$CD_Mat < result$threshold)
  # print(which(result$CD_Mat > result$threshold))
  x_trimmed = x[noninfluential_indices,]
  y_trimmed = y[noninfluential_indices]
  
  x_trimmed_pack = centralize(x_trimmed,r=T)
  x_trimmed = x_trimmed_pack$d
  fit2 = cv.lars(x_trimmed,y_trimmed,plot.it = F, use.Gram = F)
  obj2 = lars(x_trimmed,y_trimmed, use.Gram = F)
  frac2 = fit2$index[which(fit2$cv == min(fit2$cv))]
  est2 = as.vector(predict.lars(obj2,mode='fraction', s=frac2, type = 'coefficients')$coefficients)
  
  result_before = cbind(result_before,est)
  result_after = cbind(result_after,est2)
  
  print(paste0(i,'/',noe))
}

before_index = which(apply(result_before!=0,1,mean)!=0)
after_index = which(apply(result_after!=0,1,mean)!=0)

A1 = data.frame(name = gnames[before_index+1], A1 = apply(result_before!=0,1,mean)[before_index])%>% arrange(-A1)
# name freq
# 1    PPAP2C 0.89
# 2    HS3ST2 0.89
# 3      IRF3 0.32
# 4     RNF44 0.32
# 5       JIK 0.22
# 6      PTEN 0.14
# 7    SNAP25 0.03
# 8     KCNC1 0.03
# 9     PTGDS 0.03
# 10 ARHGAP15 0.03
# 11    TRPM2 0.01
# 12    MAT2B 0.01

A2 = data.frame(name = gnames[after_index+1], A2 = apply(result_after!=0,1,mean)[after_index]) %>% arrange(-A2)
# name freq
# 1       CNN3 0.97
# 2   FLJ12443 0.97
# 3       PTEN 0.97
# 4      KCNC1 0.97
# 5      SYNJ2 0.97
# 6    CGI-115 0.97
# 7      PTGDS 0.91
# 8    ADIPOR1 0.91
# 9      GTSE1 0.90
# 10      IRF3 0.86
# 11     COPZ2 0.86
# 12   C11orf2 0.73
# 13   CSNK1A1 0.70
# 14    TM9SF2 0.57
# 15      TCF4 0.46
# 16      UROD 0.33
# 17     LRRN5 0.32
# 18     MARK4 0.32
# 19    CORO1A 0.26
# 20    CELSR2 0.24
# 21       JIK 0.24
# 22 GALNACT-2 0.19
# 23       JUN 0.17
# 24   SDCCAG8 0.12
# 25  FLJ22169 0.10
# 26  HLA-DPA1 0.08
# 27  HLA-DPB1 0.07
# 28    PPAP2C 0.05
# 29      ARF4 0.04
# 30     OSTF1 0.04
# 31       ME1 0.03
# 32    FAM38B 0.02
# 33   SLC31A2 0.01
# 34 NIF3L1BP1 0.01

top10_A1 = A1$name[1:10]
top10_A2 = A2$name[1:10]

# Rajaratnam et al. results:
IL = data.frame(name = gnames[2:3601], IL = apply(read.table('datasets/nonzeroNew.txt'),2,mean)) %>% 
  filter(IL>0) %>% arrange(-IL)
top10_IL = IL$name[1:10]

fullresult = merge(A2, IL, by='name',all=T) %>% merge(A1, by='name',all=T)%>% 
  replace(is.na(.), 0) %>% arrange(-A2, -IL, -A1)

Table4 = fullresult[fullresult$name %in% union(union(top10_A1, top10_A2),top10_IL),]

# > Table4
#         name   A2   IL   A1
# 1      KCNC1 0.97 0.85 0.03
# 2       PTEN 0.97 0.73 0.14
# 3       CNN3 0.97 0.73 0.00
# 4      SYNJ2 0.97 0.71 0.00
# 5    CGI-115 0.97 0.69 0.00
# 6   FLJ12443 0.97 0.68 0.00
# 7      PTGDS 0.91 0.59 0.03
# 8    ADIPOR1 0.91 0.54 0.00
# 9      GTSE1 0.90 0.57 0.00
# 10      IRF3 0.86 0.54 0.32
# 21       JIK 0.24 0.04 0.22
# 28    PPAP2C 0.05 0.02 0.89
# 43    HS3ST2 0.00 0.07 0.89
# 121    RNF44 0.00 0.00 0.32
# 122 ARHGAP15 0.00 0.00 0.03
# 123   SNAP25 0.00 0.00 0.03


write.csv(Table4, file = "results/Table4.csv")
