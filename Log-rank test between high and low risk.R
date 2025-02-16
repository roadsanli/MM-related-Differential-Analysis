
####input ( include test factors)
res <- read.csv("MM_test.csv",header=T) #row:sample ; col:factors

##high and low risk group
res[which(res[,1]<12),1] <- 1
res[which(res[,1]>12),1] <-0
res[which(res[,1]==12),1] <-1
n <- ncol(res)

output <- matrix(0,n-1,4)
data <- res
trait <- colnames(res)[2:n]


# 构建逻辑回归模型
for(i in 1:(n-1)){
  formu <- paste0("PFS ~ ",trait[i])
  model <- glm(formu, data = data, family = binomial())
  
  # 提取OR值、P值和置信区间
  summary(model)  # 显示详细结果
  or <- exp(coef(model)[2])  # OR值
  p_value <- summary(model)$coefficients[2,4]  # P值
  conf_int <- exp(confint(model)[2, ])  # OR值的置信区间
  
  # 输出结果
  cat("OR值:", round(or, 3), "\n")
  cat("P值:", signif(p_value, 3), "\n")
  cat("95%置信区间:", round(conf_int, 3), "\n")
  or <- as.numeric(round(or, 3))
  p <- signif(p_value, 3)
  z <- as.numeric(round(conf_int, 3))
  
  output[i,] <- c(p,or,z)
  print(i)
}

######
colnames(output) <- c("P value","OR","LCI","UCI")
rownames(output) <- trait


#OR = 1：变量与结果无关联。
#OR > 1：变量与结果正相关，值越大，关联越强。
#OR < 1：变量与结果负相关，值越小，关联越强（但在负方向）。

#################forest plot
library(ggplot2)
res <- data.frame(output)
res$trait <- trait
res <- res[order(-res$P.value),]
p <- res$P.value
or <- res$OR
sq5 <- res$LCI
xq5 <- res$UCI
name <- res$trait
####
library(forestplot)
# 创建数据框
study_data <- data.frame(
  study = name,  # 研究名称
  mean = or,  # 效应值
  lower = sq5,  # 置信区间下限
  upper = xq5,  # 置信区间上限
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

# 绘制森林图
#pdf("Biochemistry and metabolic Forest Plot in Multiple_Myeloma", width=3.64,height = 2.57)
# 确保 study 按照 study_data 中的顺序排列
plot_data <- data.frame(
  study = factor(study_data$study, levels = study_data$study),  # 指定因子顺序
  mean = study_data$mean,  # 效应值
  lower = study_data$lower,  # 置信区间下限
  upper = study_data$upper,  # 置信区间上限
  p_value = as.numeric(study_data$p_value)  # p 值
)
# 添加列来标记显著性
plot_data$significance <- ifelse(plot_data$p_value < 0.05, "Significant", "Not Significant")

# 绘制森林图
pl <- ggplot(plot_data, aes(y = study, x = mean)) +
  # 添加置信区间
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, color = "black",size = 0.5) +
  # 添加圆形点，并映射填充颜色
  geom_point(aes(fill = significance), size = 2, shape = 21, color = "white") +
  # 添加参考线
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray", size = 0.45) +
  # 设置标题和轴标签
  labs(title = "PFS Forest Plot in Multiple Myeloma",
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
  filename = "PFS.pdf",
  plot = pl,
  width = 93.459 / 25.4,  # 宽度从mm转换为英寸
  height = 80 / 25.4, # 高度从mm转换为英寸
  units = "in",           # 单位为英寸
  dpi = 300               # 分辨率为300 dpi
)


