suppressMessages({
    library(tidyverse)
    library(RColorBrewer)
})

setwd("Z:/home/wanglinting/Yeast/CodeOcean/figures/svalues")

dir_result <- "../../results"

# expression --------------------------------------------------------------

data_list <- c("microarray",
               "rnaseq14",
               "concatenate")
name_list <- c("Tu et al. 2005",
               "Zheng et al. 2014",
               "Concatenated data")

data_res <- tibble(dim = NA)
for (data_name in data_list) {
    file_data <- file.path(dir_result, "temp", data_name, "loadings/s_values.csv")
    data <- read_csv(file_data, show_col_types = FALSE) %>% 
        select(-value)
    data_res <- data_res %>% 
        full_join(data, by = "dim")
}

data_plot <- data_res %>% 
    set_names(c("dim",name_list)) %>% 
    na.omit() %>% 
    filter(dim <= 5) %>% 
    pivot_longer(cols = -dim) %>% 
    mutate(name = factor(name, levels = name_list))

# 600*400
# Svalue_expr
ggplot(data_plot, aes(x = dim, y = value, fill = name)) +
    geom_col(position = position_dodge(0.9),width = 0.8, color = "grey") +
    scale_fill_manual(values = brewer.pal(4,"Set3")) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Principal eigen-component level", y = "Percentage", fill = "data") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          legend.position = c(0.8,0.8)
          )


# histone -----------------------------------------------------------------


data_list <- c(
    # "H3",
    "H3K4me3",
    "H3K9ac",
    "H3K14ac",
    "H3K18ac",
    "H3K36me3",
    "H3K56ac",
    "H4K5ac",
    "H4K16ac"
)
name_list <- c(
    # "H3",
    "H3K4me3",
    "H3K9ac",
    "H3K14ac",
    "H3K18ac",
    "H3K36me3",
    "H3K56ac",
    "H4K5ac",
    "H4K16ac"
)

data_res <- tibble(dim = NA)
for (data_name in data_list) {
    file_data <- file.path(dir_result, "temp", data_name, "loadings/s_values.csv")
    data <- read_csv(file_data, show_col_types = FALSE) %>% 
        select(-value)
    data_res <- data_res %>% 
        full_join(data, by = "dim")
}

data_plot <- data_res %>% 
    set_names(c("dim",name_list)) %>% 
    na.omit() %>% 
    filter(dim <= 5 & dim >= 1) %>% 
    pivot_longer(cols = -dim) %>% 
    mutate(name = factor(name, levels = name_list))

# 800*400
# Svalue_histone
ggplot(data_plot, aes(x = dim, y = value, fill = name)) +
    geom_col(position = position_dodge(0.9),width = 0.75, color = "grey") +
    scale_fill_manual(values = brewer.pal(8,"Set3")) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Principal eigen-component level", y = "Percentage", fill = "data") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          legend.position = c(0.75,0.8)
    ) +
    guides(fill = guide_legend(ncol = 2))


# metabolite --------------------------------------------------------------

data_list <- c(
    "Metabolite_interpolated_lc",
    "Metabolite_interpolated_gc",
    "Metabolite_lc",
    "Metabolite_gc"
    )
name_list <- c(
    "Metabolome (LC-MS, interpolated)",
    "Metabolome (GC−TFOMS, interpolated)",
    "Metabolome (LC-MS)",
    "Metabolome (GC−TFOMS)"
)


data_res <- tibble(dim = NA)
for (data_name in data_list) {
    file_data <- file.path(dir_result, "temp", data_name, "loadings/s_values.csv")
    data <- read_csv(file_data, show_col_types = FALSE) %>% 
        select(-value)
    data_res <- data_res %>% 
        full_join(data, by = "dim")
}


data_plot <- data_res %>% 
    set_names(c("dim",name_list)) %>% 
    na.omit() %>% 
    filter(dim <= 5 & dim >= 1) %>% 
    pivot_longer(cols = -dim) %>% 
    mutate(name = factor(name, levels = name_list))

# 600*400
# Svalue_metabolome
ggplot(data_plot, aes(x = dim, y = value, fill = name)) +
    geom_col(position = position_dodge(0.9),width = 0.8, color = "grey") +
    scale_fill_manual(values = brewer.pal(4,"Set2")) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Principal eigen-component level", y = "Percentage", fill = "data") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          legend.position = c(0.7,0.8)
    )
