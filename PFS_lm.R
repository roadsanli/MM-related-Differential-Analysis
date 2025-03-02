###############------------------------------------------------PFS相关因素的森林图
res <- read.csv("MM_test2.csv",header=T) #Include PFS factors of interest
res$PFS <- log(res$PFS+1) 
#res[,2:15] <- scale(res[,2:15])
n <- ncol(res)
output <- matrix(0,n-1,4)
data <- res
trait <- colnames(res)[2:n]
#-------------------------X_DN_cells_in_NDMM
# 构建逻辑回归模型
for(i in 1:(n-1)){
  formu <- paste0("PFS ~ ",trait[i])
  model <- lm(formu, data = data)
  beta <- summary(model)$coefficients[2,1]
  p_value <- summary(model)$coefficients[2,4]  # P值
  se <- summary(model)$coefficients[2,2]  # se值
  
  # 计算 95% 置信区间
  lower_CI <- beta -1.96*se  # 计算下限
  upper_CI <- beta +1.96*se   # 计算上限
  beta <- as.numeric(round(beta, 3))
  p <- signif(p_value, 3)
  output[i,] <- c(p,beta,lower_CI,upper_CI)
  print(i)
}

######
colnames(output) <- c("P value","Effect size","LCI","UCI")
rownames(output) <- trait

#################panne 的森林图
library(ggplot2)
res <- data.frame(output)
res$trait <- trait
res <- res[order(res$trait,decreasing = TRUE),]
p <- res$P.value
beta <- res$Effect.size
sq5 <- res$LCI
xq5 <- res$UCI
name <- res$trait
####
library(forestplot)
# 创建数据框
study_data <- data.frame(
  study = name,  # 研究名称
  mean = beta,  # 效应值
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
study_data <- study_data[order(study_data$mean),]
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

library(ggplot2)
library(cowplot)

# 重新整理数据，确保p值格式正确
plot_data$significance <- ifelse(plot_data$p_value < 0.05, "Significant", "Not Significant")

# 
forest_plot <- ggplot(plot_data, aes(y = study, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, color = "black", size = 0.5) +
  geom_point(aes(fill = significance), size = 2, shape = 21, color = "white") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray", size = 0.45) +
  labs(title = "PFS Forest Plot in Multiple Myeloma",
       x = "Effect size (95% CI)", y = "") +
  scale_fill_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  guides(fill = "none") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(), # 隐藏Y轴标签，避免重复
    axis.ticks.y = element_line(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

# 
p_value_table <- ggplot(plot_data, aes(y = study, x = 1, label = p_value)) +
  geom_text(size = 3) +
  labs(x = "P Value", y = "") +
  theme_void() +  # 去掉背景和坐标轴
  theme(
    axis.text.y = element_text(size = 10, hjust = 0.5), # 让Y轴文本对齐
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

# 
final_plot <- plot_grid( p_value_table, forest_plot,ncol = 2, rel_widths = c(2, 0.8))

# 
ggsave(
  filename = "PFS_with_Pvalue.pdf",
  plot = final_plot,
  width = 120 / 25.4,  # 宽度从mm转换为英寸
  height = 80 / 25.4, # 高度从mm转换为英寸
  units = "in",
  dpi = 300
)



