suppressMessages({
    library(ggrepel)
    library(openxlsx)
    library(tidyverse)
})

# ver2 --------------------------------------------------------------------

file_gene <- "gene 1 pvalue ver2.csv"
data_gene <- read_csv(file_gene, show_col_types = FALSE) %>%
    mutate(
        Loading = -Loading,
        Symbol = if_else(is.na(Symbol), ID, Symbol),
        Class = case_when(
            Loading > 0.025 & `p value` < 0.05 ~ "1A",
            Loading < -0.025 & `p value` < 0.05 ~ "1B",
            TRUE ~ "Normal"
        ),
        Class = factor(Class, levels = c("1A", "1B", "Normal"))
    ) 

wb <- createWorkbook()
addWorksheet(wb, "1A")
writeData(wb, sheet = "1A", data_gene %>% filter(Class == "1A"))
addWorksheet(wb, "1B")
writeData(wb, sheet = "1B", data_gene %>% filter(Class == "1B"))

gene_list <- bind_rows(
    data_gene %>% filter(Loading > 0.025 & `p value` < 0.05 ) %>% arrange(`p value`) %>% head(5),
    data_gene %>% filter(Loading < -0.025 & `p value` < 0.05 ) %>% arrange(`p value`) %>% head(5),
)  

p <- ggplot(data_gene, aes(Loading, -log10(`p value`))) +
    geom_point(size = 1, aes(color = Class)) +
    scale_color_manual(values = c("orange2","#5C7EE5",'gray')) +
    ylab('-log10 (P value)') +
    xlab('Loading') +
    geom_vline(xintercept = c(-0.025,0.025), lty = 2, col = "#303030") +
    geom_hline(yintercept = -log10(0.05), lty = 2, col = "#303030") +
    geom_label_repel(
        data = gene_list,
        aes(label = Symbol),
        size = 3,
        color = "black",
        show.legend = FALSE,
        min.segment.length = 0, 
        point.padding = unit(0.2, "lines") 
    ) +
    theme_classic() +
    theme(
        legend.title = element_blank(),
    )
ggsave("volcano_gene1_ver2.pdf", p, width = 7, height = 5)

file_gene <- "gene 2 pvalue ver2.csv"
data_gene <- read_csv(file_gene, show_col_types = FALSE) %>% 
    mutate(Class = case_when( 
        Loading > 0.025 & `p value` < 0.05 ~ "2A", 
        Loading < -0.025 & `p value` < 0.05 ~ "2B",  
        TRUE ~ "Normal"), 
        Class = factor(Class, levels = c("2A", "2B", "Normal"))) 

addWorksheet(wb, "2A")
writeData(wb, sheet = "2A", data_gene %>% filter(Class == "2A"))
addWorksheet(wb, "2B")
writeData(wb, sheet = "2B", data_gene %>% filter(Class == "2B"))
saveWorkbook(wb, file = "MarkGene_ver2.xlsx", overwrite = T)


gene_list <- bind_rows(
    data_gene %>% filter(Loading > 0.025 & `p value` < 0.05 ) %>% arrange(`p value`) %>% head(5),
    data_gene %>% filter(Loading < -0.025 & `p value` < 0.05 ) %>% arrange(`p value`) %>% head(5),
)  

p <- ggplot(data_gene, aes(Loading, -log10(`p value`))) +
    geom_point(size = 1, aes(color = Class)) +
    scale_color_manual(values = c("#F45B5B","#48AC48",'gray')) +
    ylab('-log10 (P value)') +
    xlab('Loading') +
    geom_vline(xintercept = c(-0.025,0.025), lty = 2, col = "#303030") +
    geom_hline(yintercept = -log10(0.05), lty = 2, col = "#303030") +
    geom_label_repel(
        data = gene_list,
        aes(label = Symbol),
        size = 3,
        color = "black",
        show.legend = FALSE,
        min.segment.length = 0, 
        point.padding = unit(0.2, "lines") 
    ) +
    theme_classic() +
    theme(
        legend.title = element_blank(),
    )
ggsave("volcano_gene2_ver2.pdf", p, width = 7, height = 5)

