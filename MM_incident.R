################.    MM incident 火山图
library(survival)
library(survminer)
MM_cox_data <- readRDS("MM_cox_data.rds")
meta <- read.csv("251metobolic.csv",header=T)
cov <- read.csv("20pc_age_sex(2).csv",header=T)

id <- MM_cox_data$ eid
meta_id <- meta[match(id,meta$eid),]
MM_cox_data_id <- MM_cox_data[match(id,MM_cox_data$eid),]
cov_id <- cov[match(id,cov$eid),]

MM_data <- meta_id
MM_data$time <- MM_cox_data_id$time
MM_data$status <- MM_cox_data_id$status
MM_data$sex <- cov_id$sex
MM_data$age <- cov_id$age
MM_data$bmi <- cov_id$bmi


#######
MM_data <- readRDS("MM_data.rds")
Pro_name <- NULL
P <- NULL 
Beta <- NULL
SE <- NULL
Lower_CI <- NULL
Upper_CI <- NULL
# Scale protein data 
MM_data[,2:252] <- apply(MM_data[,2:252], 2, scale)
for(i in 2:252){
  fit1 <- coxph(Surv(time, status) ~ MM_data[,i]+sex+age+bmi,data = MM_data)
  OR <- summary(fit1)$ coefficients[1,2]
  SE_OR <- OR * summary(fit1)$ coefficients[1,3]
  p <- summary(fit1)$ coefficients[1,5]
  beta <- summary(fit1)$ coefficients[1,1]
  lower_CI <- OR -1.96*SE_OR  # 计算下限
  upper_CI <- OR +1.96*SE_OR   # 计算上限
  Pro_name <- c(Pro_name,colnames(MM_data)[i])
  P <- c(P,p)
  Beta <- c(Beta,beta)
  SE <- c(SE,summary(fit1)$ coefficients[1,3])
  Lower_CI <- c(Lower_CI,lower_CI)
  Upper_CI <- c(Upper_CI,upper_CI)
}
log10PValue <- -log10(P)

data <- data.frame("gene_name"=Pro_name, "PValue"=P, "log10PValue"=log10PValue,"OR"=exp(Beta),"Lower_CI"=Lower_CI,"Upper_CI"=Upper_CI)
write.csv(data,paste0("MM_incident_cox.csv"),quote = FALSE, row.names = F)
min(data$OR)
max(data$OR)
pdf(paste0("MM_incident_cox.pdf"),width = 5, height = 6)##
library(ggplot2)
library(ggrepel)
library(openxlsx)
##矫正后的显著-log10P值3.701
data$Sig = ifelse(data$log10PValue > 3.701 &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                    data$OR!=1,  #
                  ifelse(data$OR > 1 ,'Up','Down'),'no')



data = data.frame(data)
table(data$Sig) #查看数据统计情况

###
p1 <- ggplot(data, aes(x =OR, y=log10PValue, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=1, size=1) +  #点的透明度、大小
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(0.7, 1.3)) + #调整点的颜色和x轴的取值范围型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  labs(x="Hazard ratio", y="-log10[P]") +  #x、y轴标签
  ggtitle("Multiple Myeloma") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(size = 16,hjust = 0.5),#标题大小，水平居中
        legend.position="right", 
        legend.title = element_blank(),
        panel.grid = element_blank(),  # 去除背景中的方格线
        axis.text = element_text(size = 10),  # 调整轴标签大小为 12
        axis.title = element_text(size = 12)  # 调整轴标题大小为 14
  ) 

###添加基因名标记###
p2 <- p1 + geom_text_repel(
  data = subset(data, data$log10PValue > 3.701 & data$OR > 1),# 可以设置跟上面不同的阈值，用数值替换即可
  aes(label = gene_name), size = 3,
  box.padding = unit(0.5, "lines"),
  point.padding = unit(0.8, "lines"), max.overlaps = 5,segment.color = "black", show.legend = FALSE )


p3 <- p2 + geom_text_repel(
  data = subset(data, data$log10PValue > 3.701 & data$OR<  1),# 可以设置跟上面不同的阈值，用数值替换即可
  aes(label = gene_name), size = 3,
  box.padding = unit(0.5, "lines"),
  point.padding = unit(0.8, "lines"), max.overlaps = 5,segment.color = "black", show.legend = FALSE )
p3+geom_hline(yintercept = 3.701, linetype = "dashed", color = "grey60", size = 0.8)+
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", size = 0.8)



dev.off()


########### 森林图

data_up <- data[data$OR>1,]
data_down <- data[data$OR<1,]
data_up <- data_up[order(data_up$ PValue),]
data_down <- data_down[order(data_down$ PValue),]
res2 <- rbind(data_up[1:4,],data_down[1:4,])

output <- res2[,c(2,4,5,6)]
colnames(output) <- c("P value","OR","LCI","UCI")
rownames(output) <- res2[,1]
#
library(ggplot2)
res <- data.frame(output)
res$trait <- res2[,1]
#res <- res[order(res$trait,decreasing = TRUE),]
p <- res$P.value
beta <- res$OR
sq5 <- res$LCI
xq5 <- res$UCI
name <- res2[,1]
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
# 添加列来标记显著性
plot_data$significance <- ifelse(plot_data$p_value < 0.05, "Significant", "Not Significant")

# 绘制森林图
forest_plot <- ggplot(plot_data, aes(y = study, x = mean)) +
  # 添加置信区间
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, color = "black",size = 0.5) +
  # 添加圆形点，并映射填充颜色
  geom_point(aes(fill = significance), size = 2, shape = 21, color = "white") +
  # 添加参考线
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray", size = 0.45) +
  # 设置标题和轴标签
  labs(title = "PFS Forest Plot in Multiple Myeloma",
       x = "Hazard ratio (95% CI)",
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
    axis.text.y = element_blank(), # 隐藏Y轴标签，避免重复
    axis.ticks.y = element_line(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

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



ggsave(
  filename = "MM_incident_forest.pdf",
  plot = final_plot,
  width = 170 / 25.4,  # 宽度从mm转换为英寸
  height = 80 / 25.4, # 高度从mm转换为英寸
  units = "in",           # 单位为英寸
  dpi = 300               # 分辨率为300 dpi
)



