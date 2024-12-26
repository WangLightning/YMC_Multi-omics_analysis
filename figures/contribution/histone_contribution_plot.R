suppressMessages({
    library(tidyverse)
})


# parameters --------------------------------------------------------------

dir_result <- "../../results"
dir_res <- file.path(dir_result, 'contribution')

names <- c(
    "H3K4me3",
    "H3K9ac",
    "H3K14ac",
    "H3K18ac",
    # "H3K36me3",
    "H3K56ac",
    "H4K5ac"
    # "H4K16ac"
)

data_interp_res <- read_csv(file.path(dir_res, "contribution_interp.csv")) %>% 
    filter(!str_detect(Data, "TSS")) %>% 
    mutate(data_name = names, CV = round(CV, digits = 3)) 

data_plot <- data_interp_res %>%
    select(-CV) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Phase", values_to = "value") %>%
    mutate(data_name = factor(data_name, levels = data_interp_res$data_name), Phase = as.factor(Phase))

p <- ggplot() +
    geom_col(data = data_plot, aes(x = data_name, y = value, fill = Phase), position = position_dodge(0.8), width = 0.7) +
    scale_fill_manual(values = c("#FFBF80", "#7D98EA", "#FF8080", "#6DBD6D")) +
    geom_label(
        data = data_interp_res, 
        aes(x = data_name, y = 0.81, label = CV)
        ) + 
    labs(y = "Relative contribution", x = NULL, fill = "Eigen-phase") +
    theme_bw() +

ggsave(
    p,
    file = "contribution_interp.pdf",
    width = 8.5,
    height = 4.5
)






