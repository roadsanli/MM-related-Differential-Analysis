
##### 01 UKB samples group-----------------------------------
res <- readRDS("hospital_ICD10_data.rds") ##UKB hospital ICD10
library(stringr)
data2 <- res[nchar(res$ p41270)==0,]
very_healthy <- data2$eid ##
result <- strsplit(res[,3], split = "[|]")
key_words <-"C90.0" 
c90_data_list <- lapply(result, function(item) {
  l <- c(1:length(item))
  c90 <- str_detect(item, key_words)
  return(l[c90])
})
data_res <- res[,-c(1:3)]

out <- data.frame("eid"=res$eid)
out$C90data <- "NA"
out$C90data <- as.Date(out$C90data,format="%Y-%m-%d")
for(i in 1:502173){
  pos <- c90_data_list[[i]]
  if(length(pos)>0){
    out$C90data[i] <- min(data_res[i,pos])
  }
}
out <- out[!is.na(out[,2]),]
data_time <- read.csv("RA_blood_follow__time.csv",header=T)
data_time <- data_time[match(out$eid,data_time$eid),]
differences <- as.numeric(difftime(out$C90data, data_time$ blood_select_time, units = "days"))/365

C90_time <- data.frame("eid"=out$eid,"c90data"=out$C90data,"differences"=differences) ##1579

plasmaid <- readRDS("protein+biomarker(CRP)+ESR+M05+M06+RA_pc(2).rds")
plasmaid <- plasmaid[,1] #52327
meta <- read.csv("251metobolic.csv",header=T)
meta_cleaned <- meta[!apply(meta[, 2:252], 1, function(x) all(is.na(x))), ]
metaid <- meta_cleaned[,1] #274257

length(intersect(very_healthy,plasmaid)) #5333  healthy and protein samples
length(intersect(very_healthy,metaid)) #29419   healthy and metabolites samples
C90_time_plasma <- C90_time[match(intersect(C90_time$eid,plasmaid),C90_time$eid),] #187
sum(C90_time_plasma$ differences < 0.5) #26
sum(C90_time_plasma$ differences > 0.5)# 161

C90_time_meta <- C90_time[match(intersect(C90_time$eid,metaid),C90_time$eid),] 
C90_time_meta <- C90_time_meta[!is.na(C90_time_meta $ differences),] ##880
sum(C90_time_meta $ differences < 0.5) #118  MM and metabolites samples
sum(C90_time_meta $ differences > 0.5)# 762  MM risk and metabolites samples(MM after metabolites)


##### 02 metabolics difference in MM risk --------------------------------------------
resd1 <- read.csv("251metobolic.csv",header=T) # include metabolites
resd <- resd1
library(glmnet)
trait <- "Multiple_myeloma"
###
case <- case_762$eid
cov <- read.csv("20pc_age_sex(2).csv",header=T) ##covariates
case_id <- intersect(resd$eid,case)
id <- case_id
resd <- resd[match(id,resd$eid),]
case_762 <- case_762[match(id,case_762$eid),]
resd$Multiple_myeloma <- case_762$ differences
###
inte <- intersect(resd$ eid,cov$eid)
res <- resd[match(inte,resd$ eid),]
pp <- match(res$ eid,cov$eid)
cov2 <- cov[pp,c(1:23,28)]     
cov2 <- cbind(res[,-1],cov2)  
##1:251 metabolites
#252 depress_TSH
#253:275 (pc1:20,sex.age.bmi)
############################################
X <- cov2
ncol <- ncol(X)
nrow <- nrow(X)
##delete 10% sample
del2 <- rep(0,nrow)
for(j in 1:nrow){
  del <- sum(is.na(X[j,1:251]))/251
  del2[j] <- ifelse(del > 0.1, '1', '0')
}
X <- X[which(del2==0),]

#####delete 10% protein
ncol <- ncol(X)
nrow <- nrow(X)
del2 <- rep(0,ncol)
for(j in 1:251){
  del <- sum(is.na(X[,j]))/nrow
  del2[j] <- ifelse(del > 0.1, '1', '0')
}
#print(sum(del2))
X <- X[,which(del2==0)]
np <- 251-sum(as.numeric(del2))
#############impute
for(i in 1:np){
  X[is.na(X[,i]),i] <- mean(X[!is.na(X[,i]),i])
  ####
  if(i%%100==0){
    print(i)}
}

# Scale protein data 
X[,1:np] <- apply(X[,1:np], 2, scale)

###
X <- X[!is.na(X$pca_1),]
######
###################
sexid <- match("sex",colnames(X))
ageid <- match("age",colnames(X))
age2id <- match("age2",colnames(X))
bmid <- match("bmi",colnames(X))
Multiple_myeloma_binaryid <- match("Multiple_myeloma",colnames(X))##


X2 <- X[,c(1:np,sexid,ageid,bmid, Multiple_myeloma_binaryid)]##

############################单变量
trait  <- "Multiple_myeloma"
newx <- rep(0,np)
P<- rep(0,np)
beta<- rep(0,np)
for(k in 1:np){
  pro = colnames(X2)[k]
  gg <- paste0(trait,"~",pro,"+sex+bmi+age")
  #gg <- paste0(pro,"~",trait,"+sex+bmi")
  resglm <- glm(gg,family="gaussian",data=X2)
  #resglm <- glm(gg,family="binomial",data=X2)  
  ss <- coef(summary(resglm)) 
  newx[k] <- pro
  P[k] <- ss[2,4] 
  beta[k] <-ss[2,1]  
  if(k%%300==0){
    print(k)}
}
#output
pro_name <- colnames(X2)[1:np]
log10PValue <- -log10(P)
data <- data.frame("gene_name"=pro_name, "PValue"=P, "log10PValue"=log10PValue,"effect"=beta)
write.csv(data,paste0(trait,"metobolic_risk.csv"),quote = FALSE, row.names = F)


#####forest plot
res <- read.csv("MM_metabolic_risk.csv",header=T)
res <- res[order(-res$PValue),]
p <- res$Adjust_PValue
beta <- res$Effect.size
se=sqrt(((beta)^2)/qchisq(p,1,lower.tail=F))
sq5 <- beta-1.96*se
xq5 <- beta+1.96*se
name <- res[,1]
####
library(forestplot)
# 创建数据框
study_data <- data.frame(
  study = factor(res[,1], levels = res[,1]),  # 研究名称
  mean = exp(beta),  # 效应值
  lower = exp(sq5),  # 置信区间下限
  upper = exp(xq5),  # 置信区间上限
  p_value = p  # p 值
)
###
study_data$mean <- round(study_data$mean, 2)  # 保留两位小数
study_data$lower <- round(study_data$lower, 2)
study_data$upper <- round(study_data$upper, 2)
study_data$p_value <- format(study_data$p_value, scientific = TRUE, digits = 2)  
####
set <- paste(study_data$lower, study_data$upper, sep = ",")
study_data$Effect95 <- paste(study_data$mean, set, sep = " (")
study_data$Effect95 <- paste(study_data$Effect95, ")", sep = "")
study_data$p_value <- as.numeric(study_data$p_value)
# 创建显著性标签
study_data$significance <- ifelse(study_data$p_value < 0.05, "Significant", "Not Significant")


# 绘制森林图
pl <- ggplot(study_data, aes(y = study, x = mean)) +
  # 添加置信区间
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, color = "black",size = 0.4) +
  # 添加圆形点，并映射填充颜色
  geom_point(aes(fill = significance), size = 2, shape = 21, color = "white") +
  # 添加参考线
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray", size = 0.45) +
  # 设置标题和轴标签
  labs(title = "biochemistry and metabolic Forest Plot in Multiple_Myeloma",
       x = "Odds Ratio (95% CI)",
       y = "") +
  # 设置填充颜色映射
  scale_fill_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  # 设置横坐标轴的区间
  #scale_x_continuous(limits = c(0, 5)) +
  # 移除颜色标签
  guides(fill = "none") +
  # 更换主题
  theme_classic() +
  theme(
    #aspect.ratio = 65 / 92,
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

ggsave(
  filename = "output.pdf",
  plot = pl,
  width = 93.459 / 25.4,  # 宽度从mm转换为英寸
  height = 65.242 / 25.4, # 高度从mm转换为英寸
  units = "in",           # 单位为英寸
  dpi = 300               # 分辨率为300 dpi
)




####03 metabolics difference in MM --------------------------------------------
library(glmnet)
trait <- "Multiple_myeloma"
###
resd <- read.csv("251metobolic.csv",header=T)
control <- read.csv("healthyid_29419.csv",header=T)
case <- read.csv("MMcase_118.csv",header=T)
case <- case[,1]
control <- control[,1]
cov <- read.csv("20pc_age_sex(2).csv",header=T)

####
case_id <- intersect(resd$eid,case)
control_id <- intersect(resd$eid,control)

id <- c(case_id,control_id)
resd <- resd[match(id,resd$eid),]
resd$Multiple_myeloma <- 0
resd$Multiple_myeloma[match(case_id,resd$eid)] <- 1
###

inte <- intersect(resd$ eid,cov$eid)
res <- resd[match(inte,resd$ eid),]
pp <- match(res$ eid,cov$eid)
cov2 <- cov[pp,c(1:23,28)]     
cov2 <- cbind(res[,-1],cov2)  
##1:251 代谢
#252 depress_TSH
#253:275 (pc1:20,sex.age.bmi)
############################################
X <- cov2
ncol <- ncol(X)
nrow <- nrow(X)
##delete 10% sample
del2 <- rep(0,nrow)
for(j in 1:nrow){
  del <- sum(is.na(X[j,1:251]))/251
  del2[j] <- ifelse(del > 0.1, '1', '0')
}
X <- X[which(del2==0),]

#####delete 10% protein
ncol <- ncol(X)
nrow <- nrow(X)
del2 <- rep(0,ncol)
for(j in 1:251){
  del <- sum(is.na(X[,j]))/nrow
  del2[j] <- ifelse(del > 0.1, '1', '0')
}
#print(sum(del2))
X <- X[,which(del2==0)]
np <- 251-sum(as.numeric(del2))
#############impute
for(i in 1:np){
  X[is.na(X[,i]),i] <- mean(X[!is.na(X[,i]),i])
  ####
  if(i%%100==0){
    print(i)}
}

# Scale protein data 
X[,1:np] <- apply(X[,1:np], 2, scale)

###
X <- X[!is.na(X$pca_1),]
######
###################
sexid <- match("sex",colnames(X))
ageid <- match("age",colnames(X))
age2id <- match("age2",colnames(X))
bmid <- match("bmi",colnames(X))
Multiple_myeloma_binaryid <- match("Multiple_myeloma",colnames(X))##
################################################多变量

X2 <- X[,c(1:np,sexid,ageid,bmid,Multiple_myeloma_binaryid)]##

############################单变量
newx <- rep(0,np)
log10PValue<- rep(0,np)
Beta<- rep(0,np)
for(k in 1:np){
  pro = colnames(X2)[k]
  gg <- paste0(pro,"~",trait,"+sex+bmi+age")
  #gg <- paste0(pro,"~",trait,"+sex+bmi")
  resglm <- glm(gg,family="gaussian",data=X2)
  #resglm <- glm(gg,family="binomial",data=X2)  
  ss <- coef(summary(resglm)) 
  newx[k] <- pro
  
  beta <- ss[match("Multiple_myeloma",rownames(ss)),1] ##
  se <- ss[match("Multiple_myeloma",rownames(ss)),2] # 
  res <- abs (beta/se)
  
  P <- ss[match("Multiple_myeloma",rownames(ss)),4]
  log10PValue[k] <- as.numeric(-log10(P))
  Beta[k] <-ss[match("Multiple_myeloma",rownames(ss)),1]  
  if(k%%300==0){
    print(k)}
}

#output
P <- 10^(-log10PValue)
pro_name <- newx
data <- data.frame("metabolic"=pro_name, "PValue"=P, "log10PValue"=log10PValue,"effect"=Beta)
trait  <- "MM"
write.csv(data,paste0(trait,"_metabolic.csv"),quote = FALSE, row.names = F)


