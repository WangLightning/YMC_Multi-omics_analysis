rm(list = ls())
suppressMessages({
  library(rlist)
  library(tidyverse)
})


# parameters --------------------------------------------------------------

# read config file
config <- list.load("../../code/config.yaml")
number_dim <- config$number_dim


# function ----------------------------------------------------------------

getMicroarrayTime <- function(x) {
    switch(
        x,
        T1 = 0.417,
        T2 = 0.833,
        T3 = 1.250,
        T4 = 1.667,
        T5 = 2.083,
        T6 = 2.500,
        T7 = 2.917,
        T8 = 3.333,
        T9 = 3.750,
        T10 = 4.167,
        T11 = 4.583,
        T12 = 5.000,
        T13 = 5.417,
        T14 = 5.833,
        T15 = 6.250,
        T16 = 6.667,
        T17 = 7.083,
        T18 = 7.500,
        T19 = 7.917,
        T20 = 8.333,
        T21 = 8.750,
        T22 = 9.167,
        T23 = 9.583,
        T24 = 10.000,
        T25 = 10.417,
        T26 = 10.833,
        T27 = 11.250,
        T28 = 11.667,
        T29 = 12.083,
        T30 = 12.500,
        T31 = 12.917,
        T32 = 13.333,
        T33 = 13.750,
        T34 = 14.167,
        T35 = 14.583,
        T36 = 15.000
    )
}

getlcTime <- function(x) {
    switch(
        x,
        T1 = 0.474,
        T2 = 0.830,
        T3 = 1.342,
        T4 = 1.727,
        T5 = 2.111,
        T6 = 2.539,
        T7 = 3.051,
        T8 = 3.521,
        T9 = 3.892,
        T10 = 4.319,
        T11 = 4.746,
        T12 = 5.202,
        T13 = 5.358,
        T14 = 5.757,
        T15 = 6.142,
        T16 = 6.569,
        T17 = 7.010,
        T18 = 7.395,
        T19 = 7.737,
        T20 = 8.036,
        T21 = 8.506,
        T22 = 8.962,
        T23 = 9.432,
        T24 = 9.902
    )
}

getgcTime <- function(x) {
    switch(
        x,
        T1 = 0.488,
        T2 = 0.801,
        T3 = 1.186,
        T4 = 1.656,
        T5 = 2.168,
        T6 = 2.681,
        T7 = 3.108,
        T8 = 3.322,
        T9 = 3.835,
        T10 = 4.347,
        T11 = 4.817,
        T12 = 5.159,
        T13 = 5.358,
        T14 = 5.800,
        T15 = 6.227,
        T16 = 6.654,
        T17 = 7.096,
        T18 = 7.594,
        T19 = 8.079,
        T20 = 8.349,
        T21 = 8.862,
        T22 = 9.346,
        T23 = 9.802,
        T24 = 10.101
    )
}


# read data ---------------------------------------------------------------

# microarray data
file_microarray_sample <-
    "../../results/temp/microarray/sample_name.txt"
file_microarray_v <- "../../results/temp/microarray/svd/V.txt"
data_microarray_sample <-
    read_table(file_microarray_sample,
               col_names = "Sample",
               show_col_types = F) %>%
    mutate(Time = sapply(Sample, getMicroarrayTime))
data_microarray_v <- read_table(file_microarray_v,
                                col_names = as.character(0:(number_dim - 1)),
                                show_col_types = F)
data_microarray <- data_microarray_sample %>%
    bind_cols(data_microarray_v)

# # Metabolite_gc
# file_metabolite_sample <-
#     "../../results/temp/Metabolite_gc/sample_name.txt"
# file_metabolite_v <- "../../results/temp/Metabolite_gc/svd/V.txt"
# data_metabolite_sample <-
#     read_table(file_metabolite_sample,
#                col_names = "Sample",
#                show_col_types = F) %>%
#     mutate(Time = sapply(Sample, getgcTime))
# data_metabolite_v <-  read_table(file_metabolite_v,
#                                  col_names = as.character(0:(number_dim - 1)),
#                                  show_col_types = F)
# data_metabolite <- data_metabolite_sample %>%
#     bind_cols(data_metabolite_v)
# 
# # Metabolite_lc
# file_metabolite_sample <-
#     "../../results/temp/Metabolite_lc/sample_name.txt"
# file_metabolite_v <- "../../results/temp/Metabolite_lc/svd/V.txt"
# data_metabolite_sample <-
#     read_table(file_metabolite_sample,
#                col_names = "Sample",
#                show_col_types = F) %>%
#     mutate(Time = sapply(Sample, getlcTime))
# data_metabolite_v <-  read_table(file_metabolite_v,
#                                  col_names = as.character(0:(number_dim - 1)),
#                                  show_col_types = F)
# data_metabolite <- data_metabolite_sample %>%
#     bind_cols(data_metabolite_v)

# # Metabolite_interpolated_gc
# file_metabolite_sample <-
#     "../../results/temp/Metabolite_interpolated_gc/sample_name.txt"
# file_metabolite_v <-
#     "../../results/temp/Metabolite_interpolated_gc/svd/V.txt"
# data_metabolite_sample <-
#     read_table(file_metabolite_sample,
#                col_names = "Sample",
#                show_col_types = F) %>%
#     mutate(Time = sapply(Sample, getgcTime))
# data_metabolite_v <-  read_table(file_metabolite_v,
#                                  col_names = as.character(0:(number_dim - 1)),
#                                  show_col_types = F)
# data_metabolite <- data_metabolite_sample %>%
#     bind_cols(data_metabolite_v)

# Metabolite_interpolated_lc
file_metabolite_sample <-
    "../../results/temp/Metabolite_interpolated_lc/sample_name.txt"
file_metabolite_v <-
    "../../results/temp/Metabolite_interpolated_lc/svd/V.txt"
data_metabolite_sample <-
    read_table(file_metabolite_sample,
               col_names = "Sample",
               show_col_types = F) %>%
    mutate(Time = sapply(Sample, getlcTime))
data_metabolite_v <-  read_table(file_metabolite_v,
                                 col_names = as.character(0:(number_dim - 1)),
                                 show_col_types = F)
data_metabolite <- data_metabolite_sample %>%
    bind_cols(data_metabolite_v)


# plot --------------------------------------------------------------------

# dim1
data_microarray_1 <- data_microarray %>%
    select(Time, Loading = `1`) %>%
    mutate(type = "Transcriptome from Tu et al. 2005",
           Dim = 'Level 1')
data_metabolite_1 <- data_metabolite %>%
    select(Time, Loading = `1`) %>%
    mutate(type = 'Metabolome',
           Dim = 'Level 1')

# dim2
data_microarray_2 <- data_microarray %>%
    select(Time, Loading = `2`) %>%
    mutate(type = "Transcriptome from Tu et al. 2005",
           Dim = 'Level 2')
data_metabolite_2 <- data_metabolite %>%
    select(Time, Loading = `2`) %>%
    mutate(type = 'Metabolome',
           Dim = 'Level 2')

data_res <- data_microarray_1 %>%
    bind_rows(data_metabolite_1) %>%
    bind_rows(data_microarray_2) %>% 
    bind_rows(data_metabolite_2) %>%
    mutate(data = factor(type, levels = c("Transcriptome from Tu et al. 2005", "Metabolome")))

# scale_factor <- (max(data_res$Loading) - min(data_res$Loading)) / 40
scale_factor <- 0.7 / 40
data_oxygen <- read_csv("oxygen_ref_YMC2005.csv", show_col_types = FALSE) %>%
    mutate(dO2_scaled = dO2 * scale_factor - 0.7)

colors <- c(
    'Dissolved oxygen' = 'grey',
    'Transcriptome from Tu et al. 2005' = '#FF7F50',
    'Metabolome' = '#00BFFF'
)
shapes <- c(
    'Dissolved oxygen' = NA,
    'Transcriptome from Tu et al. 2005' = 16,
    'Metabolome' = 15
) 

ggplot(data_res) +
    facet_wrap( ~ Dim, scales = "free_y", ncol = 1) +
    geom_line(data = data_oxygen, aes(x = Time, y = dO2_scaled, color = 'Dissolved oxygen'), linewidth = 0.25) +
    geom_line(aes(x = Time, y = Loading, color = type), linewidth = 0.35) +
    geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.2, color = 'grey45') +
    scale_y_continuous(sec.axis = sec_axis( ~ (. + 1.2) / scale_factor, name = expression(dO[2] * "(%)"))) +
    scale_color_manual(values = colors) +
    labs(x = "Time (h)", y = "Loading", color = "Data") +
    xlim(0.4, 10) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_bw() +
    theme(
        axis.text.x = element_text(size = 6, vjust = 1.2), 
        axis.text.y = element_text(size = 6, hjust = 1.2), 
        axis.title = element_text(size = 7),
        axis.line = element_blank(), 
        axis.ticks = element_line(linewidth = 0.25), 
        axis.ticks.length = unit(0.07, "cm"), 
        panel.border = element_rect(linewidth = 0.25), 
        panel.grid = element_line(linewidth = 0.25), 
        panel.spacing.y = unit(0.01, "cm"),
        strip.text = element_text(size = 7), 
        strip.background = element_blank(), 
        strip.switch.pad.grid = unit(0.01, "cm"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 6, margin = margin(0, 0, 0, -0.1, unit = "cm")),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "top"
    )
ggsave("sample_loading_metabolite&YMC2005_interpolated_lc.pdf", width = 9.3, height = 8.5, units = "cm" ) 
# ggsave("sample_loading_metabolite&YMC2005_interpolated_gc.pdf", width = 9.3, height = 8.5, units = "cm" ) 


