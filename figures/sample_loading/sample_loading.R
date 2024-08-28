suppressMessages({
    library(tidyverse)
})


setwd("Z:/home/wanglinting/Yeast/CodeOcean/figures/sample_loading")

number_dim <- 6
# colorPalette <- c("orange2", "royalblue", "red2", "green4")
# scales::show_col(colorRampPalette(c("white", "green4"))(15L))
colorPalette <- c("orange2", "#5C7EE5", "#F45B5B", "#48AC48")


# function -----------------------------------------------------------

phase <- function(Dim, Loading){
    if (Dim == "Level 1" & Loading > 0) p <- "1A"
    else if (Dim == "Level 1" & Loading < 0) p <- "1B"
    else if (Dim == "Level 2" & Loading > 0) p <- "2A"
    else if (Dim == "Level 2" & Loading < 0) p <- "2B"
    p
}
phase <- Vectorize(phase)

# 自定义标签函数，隔5个显示标签
custom_labels_5 <- function(labels) {
    new_labels <- ifelse(seq_along(labels) %% 5 == 1, labels, "")
    return(new_labels)
}

# 自定义标签函数，隔3个显示标签
custom_labels_3 <- function(labels) {
    new_labels <- ifelse(seq_along(labels) %% 3 == 1, labels, "")
    return(new_labels)
}

# 自定义标签函数，隔4个显示标签
custom_labels_2 <- function(labels) {
    new_labels <- ifelse(seq_along(labels) %% 2 == 1, labels, "")
    return(new_labels)
}

# YMC2005 --------------------------------------------------------------

file_sample <- "../../results/temp/microarray/sample_name.txt"
file_v <- "../../results/temp/microarray/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Level ", 0:(number_dim-1)), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)


# eigen component 0
column_name <- c("Sample", "Level 0")
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") 
# 900*400
# sample_loading_0_YMC2005
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    geom_bar(stat = "identity", fill = "#fd999a") + 
    labs(title = "Transcriptome from Tu et al. 2005 (Level 0)") +
    theme_bw() 
ggsave("sample_loading_0_YMC2005.svg", width = 4.1667*900, height =  4.1667*400, units = "px" ) 


# eigen component 1&2
column_name <- c("Sample", paste0("Level ", 1:2))
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))
# 1000*500
# sample_loading_12_YMC2005
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Transcriptome from Tu et al. 2005") +
    theme_bw()
ggsave("sample_loading_12_YMC2005.svg", width = 4.1667*1000, height =  4.1667*500, units = "px" ) 


# YMC2014 --------------------------------------------------------------

file_sample <- "../../results/temp/rnaseq14/sample_name.txt"
file_v <- "../../results/temp/rnaseq14/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Level ", 0:(number_dim-1)), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)


# eigen component 0
column_name <- c("Sample", "Level 0")
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") 
# 600*400
# sample_loading_0_YMC2014
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    geom_bar(stat = "identity", fill = "#fd999a") + 
    labs(title = "Transcriptome from Zheng et al. 2014 (Level 0)") +
    theme_bw() 
ggsave("sample_loading_0_YMC2014.svg", width = 4.1667*600, height =  4.1667*400, units = "px" ) 



# eigen component 1&2
column_name <- c("Sample", paste0("Level ", 1:2))
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))
# 600*400
# sample_loading_12_YMC2014
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Transcriptome from Zheng et al. 2014") +
    theme_bw()
ggsave("sample_loading_12_YMC2014.svg", width = 4.1667*600, height =  4.1667*400, units = "px" ) 


# concatenation -----------------------------------------------------------

file_sample <- "../../results/temp/concatenate/sample_name.txt"
file_v <- "../../results/temp/concatenate/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Level ", 1:number_dim), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)

column_name <- c("Sample", paste0("Level ", 1:2))


# YMC2005
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    filter(str_detect(Sample, "^M")) %>% 
    mutate(Sample = factor(paste0("T",1:36), levels = paste0("T",1:36)))  %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_x_discrete(labels = custom_labels_5) +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Transcriptome (concatenated)", subtitle = "Samples from Tu et al. 2005 (3 cycles)") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 9, vjust = -2), 
        plot.subtitle = element_text(size = 7, vjust = -2),
        axis.text.x = element_text(size = 6, vjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.text.y = element_text(size = 5, hjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.title = element_text(size = 7),
        axis.line = element_blank(), # 去掉坐标轴，否则会和边框重复导致线宽加粗
        axis.ticks = element_line(linewidth = 0.25), # 坐标刻度线宽度，但是似乎比pt更大一点
        axis.ticks.length = unit(0.07, "cm"), # 坐标刻度长度
        panel.border = element_rect(linewidth = 0.25), # 分面的框线宽度
        panel.grid = element_line(linewidth = 0.25), # 分面内格线宽度
        panel.spacing.y = unit(0.01, "cm"), # 不同分面的距离
        strip.text = element_text(size = 7), # 分面标题大小
        strip.background = element_blank(), # 去除分面标题框背景
        strip.switch.pad.grid = unit(0.01, "cm"), # 分面标题距离分面的距离
        legend.position = c(0.8,1.15),
        legend.direction = "horizontal", # 图例水平放置
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6, margin = margin(0, 0, 0, 0.05, unit = "cm")), # 使图例更紧凑
        legend.key.size = unit(0.25, "cm"), # 调整图例符号的大小
    )
# ggsave("sample_loading_12_Concatenation_YMC2005.pdf", width = 11.5, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸
ggsave("sample_loading_12_Concatenation_YMC2005_ver2.pdf", width = 11, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸



# YMC2014
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    filter(str_detect(Sample, "^R")) %>% 
    mutate(Sample = factor(paste0("T",1:16), levels = paste0("T",1:16))) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_x_discrete(labels = custom_labels_3) +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Transcriptome (concatenated)", subtitle = "Samples from Zheng et al. 2014 (1 cycle)") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 9, vjust = -2), 
        plot.subtitle = element_text(size = 7, vjust = -2),
        axis.text.x = element_text(size = 6, vjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.text.y = element_text(size = 5, hjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.title = element_text(size = 7),
        axis.line = element_blank(), # 去掉坐标轴，否则会和边框重复导致线宽加粗
        axis.ticks = element_line(linewidth = 0.25), # 坐标刻度线宽度，但是似乎比pt更大一点
        axis.ticks.length = unit(0.07, "cm"), # 坐标刻度长度
        panel.border = element_rect(linewidth = 0.25), # 分面的框线宽度
        panel.grid = element_line(linewidth = 0.25), # 分面内格线宽度
        panel.spacing.y = unit(0.01, "cm"), # 不同分面的距离
        strip.text = element_text(size = 7), # 分面标题大小
        strip.background = element_blank(), # 去除分面标题框背景
        strip.switch.pad.grid = unit(0.01, "cm"), # 分面标题距离分面的距离
        legend.position = "none",
    )
# ggsave("sample_loading_12_Concatenation_YMC2014.pdf", width = 6.5, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸
ggsave("sample_loading_12_Concatenation_YMC2014_ver2.pdf", width = 6, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸


# H3K9ac ------------------------------------------------------------------

file_sample <- "../../results/temp/H3K9ac/sample_name.txt"
file_v <- "../../results/temp/H3K9ac/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Level ", 0:(number_dim-1)), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)

# eigen component 1&2
column_name <- c("Sample", paste0("Level ", 1:2))
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))

# # 600*400
# # sample_loading_12_H3K9ac
# ggplot(data_plot, aes(x = Sample, y = Loading)) +
#     facet_wrap(~Dim, scales = "free_y", ncol = 1) +
#     geom_bar(aes(fill = Phase), stat = "identity") +
#     scale_fill_manual(values = colorPalette) +
#     labs(title = "Epigenome of hisonte modification H3K9ac") +
#     theme_bw()

ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_x_discrete(labels = custom_labels_3) +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Epigenome", subtitle = "Samples of H3K9ac (1 cycle)") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 9, vjust = -2), 
        plot.subtitle = element_text(size = 7, vjust = -2),
        axis.text.x = element_text(size = 6, vjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.text.y = element_text(size = 5, hjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.title = element_text(size = 7),
        axis.line = element_blank(), # 去掉坐标轴，否则会和边框重复导致线宽加粗
        axis.ticks = element_line(linewidth = 0.25), # 坐标刻度线宽度，但是似乎比pt更大一点
        axis.ticks.length = unit(0.07, "cm"), # 坐标刻度长度
        panel.border = element_rect(linewidth = 0.25), # 分面的框线宽度
        panel.grid = element_line(linewidth = 0.25), # 分面内格线宽度
        panel.spacing.y = unit(0.01, "cm"), # 不同分面的距离
        strip.text = element_text(size = 7), # 分面标题大小
        strip.background = element_blank(), # 去除分面标题框背景
        strip.switch.pad.grid = unit(0.01, "cm"), # 分面标题距离分面的距离
        legend.position = "none",
    )
# ggsave("sample_loading_12_H3K9ac.pdf", width = 5.5, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸
ggsave("sample_loading_12_H3K9ac_ver2.pdf", width = 5, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸



# H3K18ac ------------------------------------------------------------------

file_sample <- "../../results/temp/H3K18ac/sample_name.txt"
file_v <- "../../results/temp/H3K18ac/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Level ", 0:(number_dim-1)), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)

# eigen component 1&2
column_name <- c("Sample", paste0("Level ", 1:2))
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))

# # 600*400
# # sample_loading_12_H3K18ac
# ggplot(data_plot, aes(x = Sample, y = Loading)) +
#     facet_wrap(~Dim, scales = "free_y", ncol = 1) +
#     geom_bar(aes(fill = Phase), stat = "identity") +
#     scale_fill_manual(values = colorPalette) +
#     labs(title = "Epigenome of hisonte modification H3K18ac") +
#     theme_bw() 

ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_x_discrete(labels = custom_labels_3) +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Epigenome", subtitle = "Samples of H3K18ac (1 cycle)") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 9, vjust = -2), 
        plot.subtitle = element_text(size = 7, vjust = -2),
        axis.text.x = element_text(size = 6, vjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.text.y = element_text(size = 5, hjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.title = element_text(size = 7),
        axis.line = element_blank(), # 去掉坐标轴，否则会和边框重复导致线宽加粗
        axis.ticks = element_line(linewidth = 0.25), # 坐标刻度线宽度，但是似乎比pt更大一点
        axis.ticks.length = unit(0.07, "cm"), # 坐标刻度长度
        panel.border = element_rect(linewidth = 0.25), # 分面的框线宽度
        panel.grid = element_line(linewidth = 0.25), # 分面内格线宽度
        panel.spacing.y = unit(0.01, "cm"), # 不同分面的距离
        strip.text = element_text(size = 7), # 分面标题大小
        strip.background = element_blank(), # 去除分面标题框背景
        strip.switch.pad.grid = unit(0.01, "cm"), # 分面标题距离分面的距离
        legend.position = "none",
    )
# ggsave("sample_loading_12_H3K18ac.pdf", width = 5.5, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸
ggsave("sample_loading_12_H3K18ac_ver2.pdf", width = 5, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸

# H4K5ac ------------------------------------------------------------------
# 
# file_sample <- "../../results/temp/H4K5ac/sample_name.txt"
# file_v <- "../../results/temp/H4K5ac/svd/V.txt"
# data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
# data_v <-  read_table(file_v, 
#                       col_names = paste0("Level ", 0:(number_dim-1)), 
#                       col_types = str_dup("d", number_dim)) 
# data_v <- data_sample %>% 
#     bind_cols(data_v)
# 
# # eigen component 1&2
# column_name <- c("Sample", paste0("Level ", 1:2))
# data_plot <- data_v %>% 
#     select(all_of(column_name)) %>% 
#     pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
#     mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))
# 
# # 600*400
# # sample_loading_12_H4K5ac
# ggplot(data_plot, aes(x = Sample, y = Loading)) +
#     facet_wrap(~Dim, scales = "free_y", ncol = 1) +
#     geom_bar(aes(fill = Phase), stat = "identity") +
#     scale_fill_manual(values = colorPalette) +
#     labs(title = "Epigenome of hisonte modification H4K5ac") +
#     theme_bw() 


# Metabolite_interpolated_lc ------------------------------------------------------------------

file_sample <- "../../results/temp/Metabolite_interpolated_lc/sample_name.txt"
file_v <- "../../results/temp/Metabolite_interpolated_lc/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Level ", 0:(number_dim-1)), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)

# eigen component 1&2
column_name <- c("Sample", paste0("Level ", 1:2))
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))


ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_x_discrete(labels = custom_labels_2) +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Metabolome (interpolated)", subtitle = "Samples from LC−MS method (2 cycles)") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 9, vjust = -2), 
        plot.subtitle = element_text(size = 7, vjust = -2),
        axis.text.x = element_text(size = 6, vjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.text.y = element_text(size = 5, hjust = 1.2), # 使y轴标签稍微靠近刻度线
        axis.title = element_text(size = 7),
        axis.line = element_blank(), # 去掉坐标轴，否则会和边框重复导致线宽加粗
        axis.ticks = element_line(linewidth = 0.25), # 坐标刻度线宽度，但是似乎比pt更大一点
        axis.ticks.length = unit(0.07, "cm"), # 坐标刻度长度
        panel.border = element_rect(linewidth = 0.25), # 分面的框线宽度
        panel.grid = element_line(linewidth = 0.25), # 分面内格线宽度
        panel.spacing.y = unit(0.01, "cm"), # 不同分面的距离
        strip.text = element_text(size = 7), # 分面标题大小
        strip.background = element_blank(), # 去除分面标题框背景
        strip.switch.pad.grid = unit(0.01, "cm"), # 分面标题距离分面的距离
        legend.position = "none",
    )
# ggsave("sample_loading_12_Metabolite_interpolated_lc.pdf", width = 7.2, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸
ggsave("sample_loading_12_Metabolite_interpolated_lc_ver2.pdf", width = 6.7, height = 7, units = "cm" ) # pdf会带有0.4cm左右的白边，可以适当加大尺寸


# Metabolite_interpolated_gc ------------------------------------------------------------------

file_sample <- "../../results/temp/Metabolite_interpolated_gc/sample_name.txt"
file_v <- "../../results/temp/Metabolite_interpolated_gc/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Level ", 0:(number_dim-1)), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)

# eigen component 1&2
column_name <- c("Sample", paste0("Level ", 1:2))
data_plot <- data_v %>% 
    select(all_of(column_name)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))
# 800*400
# sample_loading_12_Metabolite_interpolated_gc
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Metabolome (interpolated)", subtitle = "Samples from GC−TFOMS method" ) +
    theme_bw() 
ggsave("sample_loading_12_Metabolite_interpolated_gc.svg", width = 4.1667*800, height =  4.1667*400, units = "px" ) 



# Metabolite_lc ------------------------------------------------------------------

file_sample <- "../../results/temp/Metabolite_lc/sample_name.txt"
file_v <- "../../results/temp/Metabolite_lc/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v,
                      col_names = paste0("Level ", 0:(number_dim-1)),
                      col_types = str_dup("d", number_dim))
data_v <- data_sample %>%
    bind_cols(data_v)

# eigen component 1&2
column_name <- c("Sample", paste0("Level ", 1:2))
data_plot <- data_v %>%
    select(all_of(column_name)) %>%
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>%
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))

# 800*400
# sample_loading_12_Metabolite_lc
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Metabolome of samples from LC−MS") +
    theme_bw()



# Metabolite_gc ------------------------------------------------------------------

file_sample <- "../../results/temp/Metabolite_gc/sample_name.txt"
file_v <- "../../results/temp/Metabolite_gc/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v,
                      col_names = paste0("Level ", 0:(number_dim-1)),
                      col_types = str_dup("d", number_dim))
data_v <- data_sample %>%
    bind_cols(data_v)

# eigen component 1&2
column_name <- c("Sample", paste0("Level ", 1:2))
data_plot <- data_v %>%
    select(all_of(column_name)) %>%
    pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>%
    mutate(Dim = as.factor(Dim), Phase = phase(Dim, Loading), Phase = factor(Phase, levels = c("1A", "1B", "2A", "2B")))

# 800*400
# sample_loading_12_Metabolite_gc
ggplot(data_plot, aes(x = Sample, y = Loading)) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_bar(aes(fill = Phase), stat = "identity") +
    scale_fill_manual(values = colorPalette) +
    labs(title = "Metabolome of samples from GC−TFOMS") +
    theme_bw()



