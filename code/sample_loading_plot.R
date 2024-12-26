suppressMessages({
    library(rlist)
    library(tidyverse)
})


# parameters --------------------------------------------------------------

dir_result <- "../results"

# read config file
config <- list.load("config.yaml")

number_dim <- config$number_dim
data_list <- c(config$data_transcriptome_list, config$data_histone_list, config$data_metabolome_list)

dir_res <- file.path(dir_result, "matrix_decomposition")
if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)


# plot -----------------------------------------------------------------

for (data_name in data_list) {

    dir_temp <- file.path(dir_result, "temp", data_name)
    file_sample <- file.path(dir_temp, "sample_name.txt")
    file_v <- file.path(dir_temp, "svd/V.txt")
    data_sample <- read_table(file_sample, col_names = "Sample", col_types = cols(Sample = col_factor()))
    
    if(data_name == "transcriptome_concatenate") {
        data_v <-  read_table(file_v, 
                            col_names = paste0("Level ", 1:number_dim), 
                            col_types = str_dup("d", number_dim))
    } else {
        data_v <-  read_table(file_v, 
                            col_names = paste0("Level ", 0:(number_dim-1)), 
                            col_types = str_dup("d", number_dim)) 
    }
    data_v <- data_sample %>% 
        bind_cols(data_v)
    
    number_sample <- nrow(data_sample)

    data_plot <- data_v %>% 
        pivot_longer(cols = where(is.numeric), names_to = "Dim", values_to = "Loading") %>% 
        mutate(Dim = as.factor(Dim), Sig = Loading > 0)
    p <- ggplot(data_plot, aes(x = Sample, y = Loading)) +
        facet_wrap(~Dim, scales = "free", ncol = 1) +
        geom_bar(aes(fill = Sig), stat = "identity") + 
        scale_fill_manual(values = c("#3d5eefe7", "#f44a3aef")) +
        theme_bw() +
        labs(x = "Loadings") +
        theme(axis.text.x = element_text(size = 7, color = "black", vjust = 1, hjust = 0.5),
                axis.text.y = element_text(size = 10, color = "black", vjust = 1, hjust = 0.5),
                axis.title.x = element_text(size = 15, color = "black", vjust = 0, hjust = 0.5),
                axis.title.y = element_text(size = 15, color = "black", vjust = 2, hjust = 0.5),
                legend.position = " none",
                strip.text = element_text(size = 15),
                strip.background = element_rect(fill = "white", colour = "black"),
        ) 
    ggsave(
        file.path(dir_res, str_glue("sample_loading_all_{data_name}.pdf")), 
        p, 
        width = number_sample * 0.4, 
        height = number_dim * 4
        )
}






