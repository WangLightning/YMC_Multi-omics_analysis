rm(list = ls())
suppressMessages({
    library(grid)
    library(ggrepel)
    library(tidyverse)
})

setwd("Z:/home/wanglinting/Yeast/CodeOcean/figures/sample_loading")



color_1A <- colorRampPalette(c("white", "orange2"))(20L)[c(8L)]
color_1B <- colorRampPalette(c("white", "royalblue"))(20L)[c(8L)]
color_2A <- colorRampPalette(c("white",  "#FF4848"))(20L)[c(8L)]
color_2B <- colorRampPalette(c("white",  "#36A336"))(20L)[c(8L)]
color_gradient_1A <- linearGradient(c("white", color_1A), x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
color_gradient_1A_2 <- linearGradient(c(color_1A, "white"), x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
color_gradient_1B <- linearGradient(c("white", color_1B), x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
color_gradient_1B_2 <- linearGradient(c(color_1B, "white"), x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
color_gradient_2A <- linearGradient(c("white", color_2A), x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
color_gradient_2A_2 <- linearGradient(c(color_2A, "white"), x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
color_gradient_2B <- linearGradient(c("white", color_2B), x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
color_gradient_2B_2 <- linearGradient(c(color_2B, "white"), x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)


p <- ggplot() +
    # 1A
    annotate("rect", xmin = 4, xmax = 5, ymin = 0, ymax = 1, fill = color_gradient_1A) +
    annotate("rect", xmin = 5, xmax = 8, ymin = 0, ymax = 1, fill = color_1A) +
    annotate("rect", xmin = 8, xmax = 9, ymin = 0, ymax = 1, fill = color_gradient_1A_2) +
    # 1B
    annotate("rect", xmin = 0, xmax = 2, ymin = 0, ymax = 1, fill = color_1B) +
    annotate("rect", xmin = 2, xmax = 4, ymin = 0, ymax = 1, fill = color_gradient_1B_2) +
    annotate("rect", xmin = 9, xmax = 11, ymin = 0, ymax = 1, fill = color_gradient_1B) +
    annotate("rect", xmin = 11, xmax = 13, ymin = 0, ymax = 1, fill = color_1B) +
    annotate("rect", xmin = 13, xmax = 15, ymin = 0, ymax = 1, fill = color_gradient_1B_2) +
    # 2A
    annotate("rect", xmin = 0, xmax = 2, ymin = -0.8, ymax = -0.2, fill = color_gradient_2A) +
    annotate("rect", xmin = 2, xmax = 5, ymin = -0.8, ymax = -0.2, fill = color_2A) +
    annotate("rect", xmin = 5, xmax = 7, ymin = -0.8, ymax = -0.2, fill = color_gradient_2A_2) +
    annotate("rect", xmin = 12, xmax = 14, ymin = -0.8, ymax = -0.2, fill = color_gradient_2A) +
    annotate("rect", xmin = 14, xmax = 15, ymin = -0.8, ymax = -0.2, fill = color_2A) +
    # 2B
    annotate("rect", xmin = 7, xmax = 8, ymin = -0.8, ymax = -0.2, fill = color_gradient_2B) +
    annotate("rect", xmin = 8, xmax = 11, ymin = -0.8, ymax = -0.2, fill = color_2B) +
    annotate("rect", xmin = 11, xmax = 12, ymin = -0.8, ymax = -0.2, fill = color_gradient_2B_2) +
    # text
    annotate("text", x = 16, y = -0.2, label = "Time", size = 2.3) +
    annotate("text", x = 6.5, y = 0.5, label = "1A", size = 3) +
    annotate("text", x = 2, y = 0.5, label = "1B", size = 3) +
    annotate("text", x = 12, y = 0.5, label = "1B", size = 3) +
    annotate("text", x = 3.5, y = -0.5, label = "2A", size = 3) +
    annotate("text", x = 9.5, y = -0.5, label = "2B", size = 3) +
    annotate("text", x = 14, y = -0.5, label = "2A", size = 3) +
    scale_y_continuous(breaks = c(-0.5, 0.5), labels = c("Level 2", "Level 1")) +
    theme_classic() +
    theme(
        axis.line.x = element_blank(),  # 去掉默认的 x 轴
        axis.line.y = element_line(linewidth = 0.25),  
        axis.text.x = element_blank(),  # 去掉 x 轴标签
        axis.ticks.x = element_blank(), # 去掉 x 轴刻度
        axis.title = element_blank(),  # 去掉轴标题
        axis.text.y = element_text(size = 7, hjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.ticks = element_line(linewidth = 0.25), # 坐标刻度线宽度，但是似乎比pt更大一点
        axis.ticks.length = unit(0.07, "cm"), # 坐标刻度长度
    ) +
    annotation_custom(
        grob = linesGrob(arrow = arrow(type = "open", length = unit(0.15, "cm")), gp = gpar(col = "black", lwd = 0.75)),
        xmin = -Inf, xmax = Inf, ymin = -0.1, ymax = -0.1
    )
ggsave("Sample_model.pdf", p, width = 9.5, height = 5.5, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸


