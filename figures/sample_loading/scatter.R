rm(list = ls())
suppressMessages({
    library(ggrepel)
    library(tidyverse)
})



number_dim <- 6

file_sample <- "../../results/temp/rnaseq14/sample_name.txt"
file_v <- "../../results/temp/rnaseq14/svd/V.txt"
data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
data_v <-  read_table(file_v, 
                      col_names = paste0("Sample loading at level ", 0:(number_dim-1)), 
                      col_types = str_dup("d", number_dim)) 
data_v <- data_sample %>% 
    bind_cols(data_v)

    
p <- ggplot(data_v, aes(y = `Sample loading at level 2`, x = `Sample loading at level 1`)) +
    geom_point(size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.25) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.25) +
    geom_text_repel(
        aes(label = Sample),
        size = 2.3,
        color = "black",
        show.legend = FALSE,
        segment.size = 0.3,  
        min.segment.length = 0, 
        point.padding = unit(0.1, "cm"), 
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
        axis.text.x = element_text(size = 6, vjust = 1.2), 
        axis.text.y = element_text(size = 6, hjust = 1.2), 
        axis.title = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.25), 
        axis.ticks.length = unit(0.07, "cm"), 
        axis.line = element_line(linewidth = 0.25)
    )
ggsave("Sample_loading_scatter.pdf", p, width = 7, height = 5.5, units = "cm" ) 



