suppressMessages({
    library(rlist)
    library(tidyverse)
})


# parameters --------------------------------------------------------------

dir_result <- "../results"

# read config file
config <- list.load("config.yaml")

data_list <- config$data_histone_list
number_dim <- config$number_dim

dir_data <- config$dir_data

dir_res <- file.path(dir_result, 'contribution')
if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)

# function -----------------------------------------------------------------

unitize = function(x) {   
    x / sqrt(sum((x) ^ 2)) 
}

sum_to_one <- function(x) {
    x / sum(x)
}

centralize = function(x) {   
    x - mean(x)
}

scalize <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)  
}

normalize <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

getRNASeqTime <- function(x) {
    switch(
        x,
        T1 = 0.356,
        T2 = 0.666,
        T3 = 0.782,
        T4 = 0.866,
        T5 = 1.017,
        T6 = 1.187,
        T7 = 1.357,
        T8 = 1.516,
        T9 = 1.686,
        T10 = 1.856,
        T11 = 2.031,
        T12 = 2.196,
        T13 = 2.366,
        T14 = 2.908,
        T15 = 3.427,
        T16 = 3.987
    )
}

getChipSeqTime <- function(x) {
    switch(
        x,
        T1 = 0.351,
        T2 = 0.618,
        T3 = 0.753,
        T4 = 0.836,
        T5 = 0.985,
        T6 = 1.160,
        T7 = 1.373,
        T8 = 1.600,
        T9 = 1.770,
        T10 = 1.996,
        T11 = 2.220,
        T12 = 2.414,
        T13 = 2.741,
        T14 = 3.183,
        T15 = 3.547,
        T16 = 4.033
    )
}


# plot -----------------------------------------------------------------

data_name <- "transcriptome_concatenate"
# data_name <- "MUREN_log_YMC2014"
dir_temp <- file.path(dir_result, "temp", data_name)
file_gene <- file.path(dir_temp, "gene_name.txt")
data_gene <- read_table(file_gene, col_names = "Symbol", col_types = "c")
file_u <- file.path(dir_temp, "svd/U.txt")
data_u <- read_table(
    file_u,
    col_names = paste0("Level ", 1:number_dim),
    # col_names = paste0("Level ", 0:(number_dim-1)),
    col_types = str_dup("d", number_dim)
    ) %>%
    bind_cols(data_gene)
data_ref <- data_u %>%
    distinct(Symbol, .keep_all = TRUE) %>%
    # 1 and 2 are divided into four vectors according to positive and negative, and sum_to_one
    mutate(
        `1A` = if_else(`Level 1` >= 0, `Level 1`, 0),
        `1B` = if_else(`Level 1` <= 0, abs(`Level 1`), 0),
        `2A` = if_else(`Level 2` >= 0, `Level 2`, 0),
        `2B` = if_else(`Level 2` <= 0, abs(`Level 2`), 0),
        # across(where(is.numeric), unitize)
        across(where(is.numeric), sum_to_one)
    ) %>% 
    select(c("Symbol", "1A", "1B", "2A", "2B"))

data_name <- "MUREN_log_YMC2014"
dir_temp <- file.path(dir_result, "temp", data_name)
file_sample <- file.path(dir_temp, "sample_name.txt")
data_sample <- read_table(file_sample, col_names = "Sample", col_types = "c")
file_v <- file.path(dir_temp, "svd/V.txt")
data_v <- read_table(
    file_v,
    col_names = paste0("Level ", 0:(number_dim-1)),
    col_types = str_dup("d", number_dim)
    ) %>%
    bind_cols(data_sample) %>%
    mutate(time = sapply(Sample, getRNASeqTime)) 
time_1 <- data_v$time
mat_1 <- data_v %>% 
    # 1 and 2 are divided into four vectors according to positive and negative, and sum_to_one
    mutate(
        `1A` = if_else(`Level 1` >= 0, `Level 1`, 0),
        `1B` = if_else(`Level 1` <= 0, abs(`Level 1`), 0),
        `2A` = if_else(`Level 2` >= 0, `Level 2`, 0),
        `2B` = if_else(`Level 2` <= 0, abs(`Level 2`), 0),
        # across(where(is.numeric), unitize)
        across(where(is.numeric), sum_to_one)
    ) %>% 
    select(c("1A", "1B", "2A", "2B")) %>% 
    as.matrix()

data_interp_res <- tibble(
    Data = character(),
    `1A` = numeric(),
    `1B` = numeric(),
    `2A` = numeric(),
    `2B` = numeric(),
    )
for (data_name in data_list) {

    # read data
    file_data <- file.path(dir_data, "histone_modification", str_glue("{data_name}.csv"))
    data <- read_csv(file_data, show_col_types = FALSE) %>% 
        distinct(Symbol, .keep_all = TRUE) %>% 
        # mutate(across(where(is.numeric), centralize)) %>%
        # mutate(across(where(is.numeric), normalize)) %>%
        # mutate(across(where(is.numeric), scalize)) %>%
        inner_join(data_ref, by = "Symbol")
    
    data1 <- data %>% 
        select(c("1A", "1B", "2A", "2B")) %>% 
        as.matrix()

    data2 <- data %>% 
        select(-c("Symbol", "1A", "1B", "2A", "2B")) %>% 
        as.matrix()

    mat_res <- crossprod(data1, data2) 
    mat_norm <- (mat_res - min(mat_res)) / (max(mat_res) - min(mat_res))
    data_res <- mat_norm %>% 
        as_tibble(rownames = 'Phase')

    write_csv(data_res, file.path(dir_res, str_glue("contribution_{data_name}.csv")))

    # barplot
    data_plot1 <- data_res %>%
        mutate(across(where(is.numeric), ~ . / sum(.) * 100)) %>%
        pivot_longer(cols = where(is.numeric), names_to = "Sample", values_to = "value") %>%
        mutate(Sample = factor(Sample, levels = unique(Sample)), Phase = as.factor(Phase))
    p1 <- ggplot(data_plot1, aes(x = Sample, y = value)) +
        geom_bar(stat = "identity", aes(fill = Phase)) +
        scale_fill_manual(values = c("#FFBF80", "#7D98EA", "#FF8080", "#6DBD6D")) +
        labs(y = "Precentage (%)") +
        theme_bw()
    ggsave(
        p1,
        file = file.path(dir_res, str_glue("contribution_{data_name}_1.pdf")),
        width = 7,
        height = 5
    )

    # barplot2
    data_plot2 <- data_res %>%
        pivot_longer(cols = where(is.numeric), names_to = "Sample", values_to = "value") %>%
        mutate(Sample = factor(Sample, levels = unique(Sample)), Phase = as.factor(Phase))
    p2 <- ggplot(data_plot2, aes(x = Sample, y = value, fill = Phase)) +
        geom_col(position = position_dodge(0.9), width = 0.8) +
        scale_fill_manual(values = c("#FFBF80", "#7D98EA", "#FF8080", "#6DBD6D")) +
        theme_bw()
    ggsave(
        p2,
        file = file.path(dir_res, str_glue("contribution_{data_name}_2.pdf")),
        width = 15,
        height = 5
    )

    # interpolation
    data_chip <- 
        # t(mat_res) %>%
        t(mat_norm) %>%
        as_tibble(rownames = 'Sample') %>%
        mutate(time = sapply(Sample, getChipSeqTime))
    time_2 <- data_chip$time
    mat_2 <- data_chip %>% 
        select(c("1A", "1B", "2A", "2B")) %>% 
        as.matrix()
    
    row_result <- numeric()
    for (i in c("1A", "1B", "2A", "2B")){
        # linear interpolation
        interp_fun <- approxfun(time_2, mat_2[,i], rule = 2)  
        interp_res <- interp_fun(time_1) %*% mat_1[,i]
        row_result <- c(row_result, interp_res)
    }
    
    data_interp_res <- data_interp_res %>% 
        add_row(
            Data = data_name,
            `1A` = row_result[1],
            `1B` = row_result[2],
            `2A` = row_result[3],
            `2B` = row_result[4]
        )
}

write_csv(data_interp_res, file.path(dir_res, "contribution_interp.csv"))






