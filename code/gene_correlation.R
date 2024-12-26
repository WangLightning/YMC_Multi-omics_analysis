suppressMessages({
    library(corrplot)
    library(rlist)
    library(openxlsx)
    library(tidyverse)
})


# parameters --------------------------------------------------------------

dir_result <- "../results"

# read config file
config <- list.load("config.yaml")

data_transcriptome_list <- config$data_transcriptome_list
data_histone_list <- config$data_histone_list
number_dim <- config$number_dim

dir_res <- file.path(dir_result, "gene_correlation")
if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)


# concatenated vs other ---------------------------------------------------------------

data_list <- c(data_transcriptome_list, data_histone_list)

# concatenated data as reference
data_name <- "transcriptome_concatenate"
file_ref_gene <- file.path(dir_result, "temp", data_name ,"gene_name.txt")
file_ref_u <- file.path(dir_result, "temp", data_name ,"svd/U.txt")
data_ref_gene <- read_table(file_ref_gene, col_names = "Symbol", col_types = "c")
data_ref_u <- read_table(file_ref_u, 
                         col_names = paste0(data_name, "_", as.character(1:number_dim)), 
                         col_types = str_dup("d", number_dim)) %>% 
    bind_cols(data_ref_gene)

wb <- createWorkbook()
for (data_name in data_list) {

    if(data_name == "transcriptome_concatenate") next

    dir_temp <- file.path(dir_result, "temp", data_name)
    file_gene <- file.path(dir_temp, "gene_name.txt")
    data_gene <- read_table(file_gene, col_names = "Symbol", col_types = "c")
    file_u <- file.path(dir_temp, "svd", "U.txt")
    data_u <- read_table(file_u, 
                         col_names = paste0(data_name, "_", as.character(0:(number_dim-1))), 
                         col_types = str_dup("d", number_dim)) %>% 
        bind_cols(data_gene)

    data_common_gene <- intersect(data_ref_gene, data_gene)

    data_1 <- data_common_gene %>% 
        left_join(data_ref_u, by = "Symbol") %>%
        arrange(Symbol) %>% 
        distinct(Symbol, .keep_all = TRUE) %>% 
        select(-1)

    data_2 <- data_common_gene %>% 
        left_join(data_u, by = "Symbol") %>%
        arrange(Symbol) %>% 
        distinct(Symbol, .keep_all = TRUE) %>% 
        select(-c(1,2))

    res <- data_1 %>%
        cor(data_2) %>%
        abs() %>%
        round(digits = 3) %>%
        as_tibble()

    addWorksheet(wb, data_name)
    writeData(wb, sheet = data_name, res)
}

saveWorkbook(wb, file = file.path(dir_res, "cor_transcriptome_concatenate.xlsx"), overwrite = T)

# pairwise  ---------------------------------------------------------------

data_name <- "transcriptome_concatenate"
file_gene <- file.path(dir_result, "temp", data_name ,"gene_name.txt")
file_u <- file.path(dir_result, "temp", data_name ,"svd/U.txt")
data_gene <- read_table(file_gene, col_names = "Symbol", col_types = "c")
data_u <- read_table(file_u, 
                     col_names = paste0(data_name, "_", as.character(1:number_dim)), 
                     col_types = str_dup("d", number_dim)) 
data <- data_gene %>% 
    bind_cols(data_u) %>% 
    distinct(Symbol, .keep_all = TRUE)

for (data_name in data_list) {
    
    if(data_name == "transcriptome_concatenate") next
    
    file_gene <- file.path(dir_result, "temp", data_name, "gene_name.txt")
    data_gene <- read_table(file_gene, col_names = "Symbol", col_types = "c")
    file_u <- file.path(dir_result, "temp", data_name, "svd/U.txt")
    data_u <- read_table(file_u, 
                         col_names = paste0(data_name, "_", as.character(0:(number_dim-1))), 
                         col_types = str_dup("d", number_dim)) 
    data_u <- data_gene %>% 
        bind_cols(data_u) %>% 
        distinct(Symbol, .keep_all = TRUE)
    data <- data %>% 
        inner_join(data_u, by = "Symbol")
}


wb <- createWorkbook()
for (i in 1:(number_dim-1)) {
    
    name <- as.character(i)
    data_res <- data %>% 
        select(ends_with(name))
    r <- colnames(data_res)
    res <- data_res %>%
        cor(data_res) %>%
        abs() %>%
        round(digits = 3) %>%
        as_tibble() %>% 
        mutate(data = r) %>% 
        select(data, everything())
    
    addWorksheet(wb, name)
    writeData(wb, sheet = name, res)
}
saveWorkbook(wb, file = file.path(dir_res, "cor_pairwise.xlsx"), overwrite = T)

# heat map -------------------------------------------------

datas <- c(
    data_transcriptome_list,
    data_histone_list
)
names <- c(
    "Transcriptome",
    "H3K4me3",
    "H3K9ac",
    "H3K14ac",
    "H3K18ac",
    "H3K36me3",
    "H3K56ac",
    "H4K5ac",
    "H4K16ac"
)
names <- c(
    "Tu et al. 2005 (Transcriptome)",
    "Zheng et al. 2014 (Transcriptome)",
    "Concatenated data (Transcriptome)",
    "H3K4me3 (Epigenome)",
    "H3K9ac (Epigenome)",
    "H3K14ac (Epigenome)",
    "H3K18ac (Epigenome)",
    "H3K36me3 (Epigenome)",
    "H3K56ac (Epigenome)",
    "H4K5ac (Epigenome)",
    "H4K16ac (Epigenome)"
)

# corr 1
name <- as.character(1)
data_res <- data %>% 
    select(ends_with(name)) %>% 
    select(starts_with(datas)) %>% 
    set_names(names)
res <- data_res %>%
    cor(data_res) %>%
    abs() %>%
    round(digits = 3) 


pdf(file.path(dir_res, "cor_1.pdf"), width = 8, height = 8)
corrplot(res, 
         col.lim = c(0,1),
         method = "number",
         tl.col = "black", 
         type = "lower"
)
dev.off()

# corr 2
name <- as.character(2)
data_res <- data %>% 
    select(ends_with(name)) %>% 
    select(starts_with(datas)) %>% 
    set_names(names)
res <- data_res %>%
    cor(data_res) %>%
    abs() %>%
    round(digits = 3) 


pdf(file.path(dir_res, "cor_2.pdf"), width = 8, height = 8)
corrplot(res, 
         col.lim = c(0,1),
         method = "number",
         tl.col = "black", 
         type = "lower"
)
dev.off()

