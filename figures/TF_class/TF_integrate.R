suppressMessages({
    library(tidyverse)
    library(openxlsx)
    library(RColorBrewer)
})


colors <- c(
    "#8DD3C7",
    "#BEBADA",
    "#FDB462",
    "#FB8072",
    "#80B1D3",
    "#E3E36A",
    "#B3DE69",
    "#FCCDE5",
    "#4DAF4A",
    "#D9D9D9",
    "#377EB8",
    "#67D5B5",
    "#84B1ED",
    "#FFBC42",
    "#BC80BD",
    "#E41A1C",
    "#CCEBC5",
    "#FFED6F",
    "#984EA3",
    "#FF7F00"
)


data <- read.xlsx("../../results/concatenate/BASE/BASE_gene1.xlsx") %>%
    as_tibble() %>%
    rename(joint_1A = pvalue.up, joint_1B = pvalue.down) %>%
    left_join(read.xlsx("../../results/concatenate/BASE/BASE_gene2.xlsx") %>% as_tibble(), by = c("motif", "TF", "class")) %>%
    rename(joint_2A = pvalue.up, joint_2B = pvalue.down) %>% 
    mutate(class = as.factor(class))
data_class <- data %>% 
    distinct(TF, .keep_all = T) %>% 
    count(class) %>% 
    rename(all = n) %>% 
    arrange(desc(all)) %>% 
    mutate(color = colors)
sum(data_class$all) # 172


data_1A <- data %>% 
    filter(joint_1A <= 0.05) %>% 
    select(TF, class, p_value = joint_1A) %>% 
    arrange(p_value) %>% 
    distinct(TF, .keep_all = T) 
data_1A_class <- data_1A %>% 
    count(class)
data_1A_class_test <- data_1A_class %>% 
    left_join(data_class) %>% 
    mutate(n2 = all - n,
           n3 = nrow(data_1A) - n,
           n4 = sum(data_class$all) - n - n2 - n3)


data_1B <- data %>% 
    filter(joint_1B <= 0.05) %>% 
    select(TF, class, p_value = joint_1B) %>% 
    arrange(p_value) %>% 
    distinct(TF, .keep_all = T) 
data_1B_class <- data_1B %>% 
    count(class)
data_1B_class_test <- data_1B_class %>% 
    left_join(data_class) %>% 
    mutate(n2 = all - n,
           n3 = nrow(data_1B) - n,
           n4 = sum(data_class$all) - n - n2 - n3)


data_2A <- data %>% 
    filter(joint_2A <= 0.05) %>% 
    select(TF, class, p_value = joint_2A) %>% 
    arrange(p_value) %>% 
    distinct(TF, .keep_all = T) 
data_2A_class <- data_2A %>% 
    count(class) 
data_2A_class_test <- data_2A_class %>% 
    left_join(data_class) %>% 
    mutate(n2 = all - n,
           n3 = nrow(data_2A) - n,
           n4 = sum(data_class$all) - n - n2 - n3)


data_2B <- data %>% 
    filter(joint_2B <= 0.05) %>% 
    select(TF, class, p_value = joint_2B) %>% 
    arrange(p_value) %>% 
    distinct(TF, .keep_all = T) 
data_2B_class <- data_2B %>% 
    count(class)
data_2B_class_test <- data_2B_class %>% 
    left_join(data_class) %>% 
    mutate(n2 = all - n,
           n3 = nrow(data_2B) - n,
           n4 = sum(data_class$all) - n - n2 - n3)


wb <- createWorkbook()
addWorksheet(wb, "1A")
writeData(wb, "1A", data_1A)
addWorksheet(wb, "1B")
writeData(wb, "1B", data_1B)
addWorksheet(wb, "2A")
writeData(wb, "2A", data_2A)
addWorksheet(wb, "2B")
writeData(wb, "2B", data_2B)
saveWorkbook(wb, file = "TFs.xlsx", overwrite = T)


# test --------------------------------------------------------------------

C6_1A <- matrix(c(34, 13, 35, 90), ncol = 2)
chisq.test(C6_1A) # 3.177e-07

C2H2_1A <- matrix(c(23, 11, 46, 92), ncol = 2)
chisq.test(C2H2_1A) # 0.0005378

C6_1B <- matrix(c(20, 27, 14, 111), ncol = 2)
chisq.test(C6_1B) # 1.153e-05

Tr_1B <- matrix(c(4, 4, 30, 134), ncol = 2)
fisher.test(Tr_1B) # 0.04989

C2H2_2A <- matrix(c(11, 23, 20, 118), ncol = 2)
chisq.test(C2H2_2A) # 0.02942

C6_2A <- matrix(c(9, 38, 22, 103), ncol = 2)
chisq.test(C6_2A) # 0.9897

C2H2_2B <- matrix(c(6, 28, 13, 125), ncol = 2)
fisher.test(C2H2_2B) # 0.2182

C6_2B <- matrix(c(6, 41, 13, 112), ncol = 2)
fisher.test(C6_2B) # 0.7852

HE_2B <- matrix(c(3, 0, 16, 153), ncol = 2)
fisher.test(HE_2B) # 0.001163


# plot --------------------------------------------------------------------

data_1 <- data_1A_class %>%
    filter(class != 'Unknown') %>%
    arrange(desc(n)) %>%
    slice_head(n = 5)
data_2 <- data_1A_class %>%
    anti_join(data_1, by = 'class') %>%
    summarize(class = 'Unknown', n = sum(n))
data_plot <- bind_rows(data_1, data_2) %>% 
    left_join(data_class) %>% 
    mutate(class = ifelse(class == 'Unknown', 'Others & Unknown', class),
           class = as.factor(class)) 
# 500*200
# ggplot(data_plot, aes(x = class, y = n, fill = class)) +
#     geom_bar(stat = "identity", color = "grey") +
#     scale_fill_manual(values = (data_plot %>% arrange(class))$color) +
#     scale_x_discrete(limits = rev(data_plot$class)) +
#     theme_minimal() +
#     coord_flip() +
#     labs(
#         y = "Number",
#     ) +
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_text(size = 10, color = "black"),
#           axis.text.x = element_text(size = 8, color = "black"),
#           axis.title.x = element_text(size = 10, color = "black"),
#           legend.position = "none"
#     )
ggplot(data_plot, aes(x = class, y = n, fill = class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = (data_plot %>% arrange(class))$color) +
    scale_x_discrete(limits = rev(data_plot$class)) +
    theme_minimal() +
    coord_flip() +
    labs(y = "Number") +
    theme(
        axis.text.x = element_text(size = 6, color = "black"), 
        axis.text.y = element_text(size = 7, color = "black"), 
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
    )
ggsave("TF_class_1A.pdf", width = 8, height = 4, units = "cm" ) 


data_1 <- data_1B_class %>%
    filter(class != 'Unknown') %>%
    arrange(desc(n)) %>%
    slice_head(n = 5)
data_2 <- data_1B_class %>%
    anti_join(data_1, by = 'class') %>%
    summarize(class = 'Unknown', n = sum(n))
data_plot <- bind_rows(data_1, data_2) %>% 
    left_join(data_class) %>% 
    mutate(class = ifelse(class == 'Unknown', 'Others & Unknown', class),
           class = as.factor(class))
# # 500*200
# ggplot(data_plot, aes(x = class, y = n, fill = class)) +
#     geom_bar(stat = "identity", color = "grey") +
#     scale_fill_manual(values = (data_plot %>% arrange(class))$color) +
#     scale_x_discrete(limits = rev(data_plot$class)) +
#     theme_minimal() +
#     coord_flip() +
#     # ylim(c(0, 13)) +
#     labs(
#         y = "Number",
#     ) +
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_text(size = 10, color = "black"),
#           axis.text.x = element_text(size = 8, color = "black"),
#           axis.title.x = element_text(size = 10, color = "black"),
#           legend.position = "none"
#     )
ggplot(data_plot, aes(x = class, y = n, fill = class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = (data_plot %>% arrange(class))$color) +
    scale_x_discrete(limits = rev(data_plot$class)) +
    theme_minimal() +
    coord_flip() +
    labs(y = "Number") +
    theme(
        axis.text.x = element_text(size = 6, color = "black"), 
        axis.text.y = element_text(size = 7, color = "black"), 
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
    )
ggsave("TF_class_1B.pdf", width = 8, height = 4, units = "cm" ) 


data_1 <- data_2A_class %>%
    filter(class != 'Unknown') %>%
    arrange(desc(n)) %>%
    slice_head(n = 5)
data_2 <- data_2A_class %>%
    anti_join(data_1, by = 'class') %>%
    summarize(class = 'Unknown', n = sum(n))
data_plot <- bind_rows(data_1, data_2) %>% 
    left_join(data_class) %>% 
    mutate(class = ifelse(class == 'Unknown', 'Others & Unknown', class),
           class = as.factor(class))
# # 500*200
# ggplot(data_plot, aes(x = class, y = n, fill = class)) +
#     geom_bar(stat = "identity", color = "grey") +
#     scale_fill_manual(values = (data_plot %>% arrange(class))$color) +
#     scale_x_discrete(limits = rev(data_plot$class)) +
#     theme_minimal() +
#     coord_flip() +
#     # ylim(c(0, 11)) +
#     labs(
#         y = "Number",
#     ) +
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_text(size = 10, color = "black"),
#           axis.text.x = element_text(size = 8, color = "black"),
#           axis.title.x = element_text(size = 10, color = "black"),
#           legend.position = "none"
#     )
ggplot(data_plot, aes(x = class, y = n, fill = class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = (data_plot %>% arrange(class))$color) +
    scale_x_discrete(limits = rev(data_plot$class)) +
    theme_minimal() +
    coord_flip() +
    labs(y = "Number") +
    theme(
        axis.text.x = element_text(size = 6, color = "black"), 
        axis.text.y = element_text(size = 7, color = "black"), 
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
    )
ggsave("TF_class_2A.pdf", width = 8, height = 4, units = "cm" ) 

data_1 <- data_2B_class %>%
    filter(class != 'Unknown') %>%
    arrange(desc(n)) %>%
    slice_head(n = 5)
data_2 <- data_2B_class %>%
    anti_join(data_1, by = 'class') %>%
    summarize(class = 'Unknown', n = sum(n))
data_plot <- bind_rows(data_1, data_2) %>% 
    left_join(data_class) %>% 
    mutate(class = ifelse(class == 'Unknown', 'Others & Unknown', class),
           class = as.factor(class))
# # 500*200
# ggplot(data_plot, aes(x = class, y = n, fill = class)) +
#     geom_bar(stat = "identity", color = "grey") +
#     scale_fill_manual(values = (data_plot %>% arrange(class))$color) +
#     scale_x_discrete(limits = rev(data_plot$class)) +
#     theme_minimal() +
#     coord_flip() +
#     labs(
#         y = "Number",
#     ) +
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_text(size = 10, color = "black"),
#           axis.text.x = element_text(size = 8, color = "black"),
#           axis.title.x = element_text(size = 10, color = "black"),
#           legend.position = "none"
#     )
ggplot(data_plot, aes(x = class, y = n, fill = class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = (data_plot %>% arrange(class))$color) +
    scale_x_discrete(limits = rev(data_plot$class)) +
    theme_minimal() +
    coord_flip() +
    labs(y = "Number") +
    theme(
        axis.text.x = element_text(size = 6, color = "black"), 
        axis.text.y = element_text(size = 7, color = "black"), 
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
    )
ggsave("TF_class_2B.pdf", width = 8, height = 4, units = "cm" ) 





