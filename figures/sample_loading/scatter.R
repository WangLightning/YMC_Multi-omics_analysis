rm(list = ls())
suppressMessages({
    library(ggrepel)
    library(tidyverse)
})


setwd("Z:/home/wanglinting/Yeast/CodeOcean/figures/sample_loading")

number_dim <- 6

file_sample <- "../../results/temp/rnaseq14/sample_name.txt"
file_v <- "../../results/temp/rnaseq14/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Sample loading at level ", 0:(number_dim-1)), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)

# # 450*450
# ggplot(data_v, aes(y = `Principal sample-eigenvector 1`, x = `Principal sample-eigenvector 2`)) +
#     geom_point() +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
#     geom_label_repel(
#         aes(label = Sample),
#         size = 3,
#         color = "black",
#         show.legend = FALSE,
#         min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
#         point.padding = unit(0.2, "lines") # 线和点之间的距离
#     ) +
#     geom_curve(
#         aes(
#             x = 0,
#             y = 0.05,
#             xend = -0.05,
#             yend = 0
#         ),
#         curvature = -2.5,
#         linewidth = 0.7,
#         arrow = arrow(length = unit(0.2, "cm"))
#     ) +
#     coord_fixed(ratio = 1) + 
#     xlim(c(-0.45, 0.45)) + 
#     ylim(c(-0.45, 0.45)) +
#     geom_text(aes(x = 0, y = 0.1, label = "1A"), color = "royalblue") +
#     geom_text(aes(x = 0, y = -0.1, label = "1B"), color = "orange2") +
#     geom_text(aes(x = 0.1, y = 0, label = "2A"), color = "#FF4848") +
#     geom_text(aes(x = -0.1, y = 0, label = "2B"), color = "#36A336") +
#     theme_classic()
    
p <- ggplot(data_v, aes(y = `Sample loading at level 2`, x = `Sample loading at level 1`)) +
    geom_point(size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.25) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.25) +
    geom_text_repel(
        aes(label = Sample),
        size = 2.3,
        color = "black",
        show.legend = FALSE,
        segment.size = 0.3,  # 单位为mm
        min.segment.length = 0, # 0始终为标签添加指引线段；若不想添加线段，则改为Inf
        point.padding = unit(0.1, "cm"), # 线和点之间的距离
    ) +
    geom_curve(
        aes(x = 0, y = 0.07, xend = -0.07, yend = 0),
        curvature = -2.5,
        linewidth = 0.25,
        arrow = arrow(length = unit(0.12, "cm"))
    ) +
    coord_fixed(ratio = 1) + 
    xlim(c(-0.5, 0.5)) + 
    ylim(c(-0.4, 0.4)) +
    annotate("rect", ymin = -0.03, ymax = 0.03, xmin = 0.11, xmax = 0.19, fill = "white") +
    annotate("rect", ymin = -0.03, ymax = 0.03, xmin = -0.19, xmax = -0.11, fill = "white") +
    annotate("rect", ymin = 0.12, ymax = 0.18, xmin = -0.04, xmax = 0.04, fill = "white") +
    annotate("rect", ymin = -0.18, ymax = -0.12, xmin = -0.04, xmax = 0.04, fill = "white") +
    annotate("text", y = 0, x = 0.15, label = "1A", color = "orange2", size = 3, fontface = "bold") +
    annotate("text", y = 0, x = -0.15, label = "1B", color = "royalblue", size = 3, fontface = "bold") +
    annotate("text", y = 0.15, x = 0, label = "2A", color = "#FF4848", size = 3, fontface = "bold") +
    annotate("text", y = -0.15, x = 0, label = "2B", color = "#36A336", size = 3, fontface = "bold") +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 6, vjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.text.y = element_text(size = 6, hjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.title = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.25), # 坐标刻度线宽度，但是似乎比pt更大一点
        axis.ticks.length = unit(0.07, "cm"), # 坐标刻度长度
        axis.line = element_line(linewidth = 0.25)
    )
ggsave("Sample_loading_scatter.pdf", p, width = 7, height = 5.5, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸



