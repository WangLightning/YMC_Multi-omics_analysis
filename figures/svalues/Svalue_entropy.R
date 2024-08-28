rm(list = ls())
suppressMessages({
    library(tidyverse)
})

setwd("Z:/home/wanglinting/Yeast/CodeOcean/figures/svalues")
dir_result <- "../../results"

entropy <- function(x){
    s <- sum(x)
    p <- x/s 
    e <- -p*log(p)
    sum(e[!is.na(e)])
}

# 每层单独计算 -----------------------------------------------------------------

file_data <- file.path(dir_result, "temp", "H3K9ac", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
entropy(data$ratio) # 1.601024

file_data <- file.path(dir_result, "temp", "H3K18ac", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
entropy(data$ratio) # 1.116343

file_data <- file.path(dir_result, "temp", "concatenate", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
entropy(data$ratio) # 1.792348

file_data <- file.path(dir_result, "temp", "Metabolite_interpolated_lc", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
entropy(data$ratio) # 2.331486


# 第三层之后合并 -----------------------------------------------------------------

file_data <- file.path(dir_result, "temp", "H3K9ac", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
res <- c(data$ratio[1:2], 1-sum(data$ratio[1:2]))
entropy(res) # 0.9861619

file_data <- file.path(dir_result, "temp", "H3K18ac", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
res <- c(data$ratio[1:2], 1-sum(data$ratio[1:2]))
entropy(res) # 0.8479013

file_data <- file.path(dir_result, "temp", "concatenate", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
res <- c(data$ratio[1:2], 1-sum(data$ratio[1:2]))
entropy(res) # 1.043027

file_data <- file.path(dir_result, "temp", "Metabolite_interpolated_lc", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
res <- c(data$ratio[1:2], 1-sum(data$ratio[1:2]))
entropy(res) # 1.011868


# plot --------------------------------------------------------------------

rm(list = ls())

# 打开PDF设备，指定文件名和页面大小
# 页面大小单位为英寸，默认字号为12pt，使用cm除以1.48，再在AI中缩放58.33%即可
pdf("Svalue_entropy.pdf", width = 12 / 1.48, height = 5.5 / 1.48)  

# x <- 3:1
# y <- c(86.0649, 75.55277, 50.59236)
# y2 <- c(67.8563, 49.4773, 34.1819)

x <- 4:1
y <- c(86.0649, 70.24324, 75.55277, 50.59236)
y2 <- c(67.8563, 54.2053859133921, 49.4773, 34.1819)

# 设置边距 c(bottom, left, top, right)
par(mar = c(5, 12, 2, 2))  # 默认值通常是c(5, 4, 4, 2) + 0.1
# 画布
# plot(x, y, type = "n", ylim = c(0.5,3.5), xlim = c(-10,110), axes = F, ann = F)
plot(x, y, type = "n", ylim = c(0.5,5), xlim = c(-5,125), axes = F, ann = F)

# x轴
axis(1, at = seq(0, 100, 20))
# x轴标签
mtext("Contribution proportion (%)", side = 1, line = 2.5)
# y轴
# axis(2, at = 3:1, labels = c("Epigenome of H3K18ac", "Concatenated transcriptome", "Metabolome"), las = 1)
axis(2, at = 4:1, labels = c("Epigenome for H3K18ac", "Epigenome for H3K9ac", "Transcriptome", "Metabolome"), las = 1)
# axis(4, at = 4:1, labels = c("1.12", "1.60", "1.79", "2.33"), las = 1)
# mtext("Entropy", side = 4, line = 2.5)
box()

text(114, 4.6, labels = "Eigen-entropy")
text(115, 1:4, labels = c("2.331", "1.792", "1.601", "1.116"))

# 黑色外框
segments(0, x, 100, x, lwd = 20)
# 白色填充
segments(0, x, 100, x, lwd = 16, col = "white")
# 红色温度
segments(0, x, y, x, lwd = 16, col = "royalblue")
# 红色温度
segments(0, x, y2, x, lwd = 16, col = "#F08080")

legend("top",         # 图例位置
       legend = c("Level 1", "Level 2"),   # 图例文本
       col = c("#F08080", "royalblue"), # 对应的颜色
       lwd = 16, # 线宽
       horiz = TRUE, # 横向排列
       bty = "n")         


# 关闭PDF设备
dev.off()
