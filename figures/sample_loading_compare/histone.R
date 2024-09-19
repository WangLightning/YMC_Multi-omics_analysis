suppressMessages({
  library(rlist)
  library(tidyverse)
})


# parameters --------------------------------------------------------------

# read config file
config <- list.load("../../code/config.yaml")
number_dim <- config$number_dim


# functions ---------------------------------------------------------------

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

# getChipSeqTime <- function(x) {
#   switch(
#     x,
#     T1 = 0.345,
#     T2 = 0.610,
#     T3 = 0.745,
#     T4 = 0.888,
#     T5 = 1.106,
#     T6 = 1.362,
#     T7 = 1.616,
#     T8 = 1.780,
#     T9 = 1.912,
#     T10 = 2.031,
#     T11 = 2.258,
#     T12 = 2.400,
#     T13 = 2.886,
#     T14 = 3.251,
#     T15 = 3.647,
#     T16 = 4.017
#   )
# }
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


# read data ---------------------------------------------------------------

# rnaseq data
file_rnaseq_sample <- "../../results/temp/rnaseq14/sample_name.txt"
file_rnaseq_v <- "../../results/temp/rnaseq14/svd/V.txt"
data_rnaseq_sample <- read_table(file_rnaseq_sample, col_names = "Sample") %>%
    mutate(Time = sapply(Sample, getRNASeqTime))
data_rnaseq_v <- read_table(file_rnaseq_v, 
                            col_names = as.character(0:(number_dim-1)))
data_rnaseq <- data_rnaseq_sample %>%
    bind_cols(data_rnaseq_v)

# H3K9ac
file_H3K9ac_sample <- "../../results/temp/H3K9ac/sample_name.txt"
file_H3K9ac_v <- "../../results/temp/H3K9ac/svd/V.txt"
data_H3K9ac_sample <- read_table(file_H3K9ac_sample, col_names = "Sample") %>%
    mutate(Time = sapply(Sample, getChipSeqTime))
data_H3K9ac_v <-  read_table(file_H3K9ac_v, 
                              col_names = as.character(0:(number_dim-1)))
data_H3K9ac <- data_H3K9ac_sample %>%
    bind_cols(data_H3K9ac_v)


# H3K18ac
file_H3K18ac_sample <- "../../results/temp/H3K18ac/sample_name.txt"
file_H3K18ac_v <- "../../results/temp/H3K18ac/svd/V.txt"
data_H3K18ac_sample <- read_table(file_H3K18ac_sample, col_names = "Sample") %>%
    mutate(Time = sapply(Sample, getChipSeqTime))
data_H3K18ac_v <-  read_table(file_H3K18ac_v, 
                             col_names = as.character(0:(number_dim-1)))
data_H3K18ac <- data_H3K18ac_sample %>%
    bind_cols(data_H3K18ac_v)


# plot --------------------------------------------------------------------

# dim1
data_rnaseq_1 <- data_rnaseq %>%
    select(Time, Loading = `1`) %>%
    mutate(type = "Transcriptome from Zheng et al. 2014",
           Dim = 'Level 1')
data_H3K9ac_1 <- data_H3K9ac %>%
    select(Time, Loading = `1`) 
data_H3K9ac_1 <- data_H3K9ac_1 %>%
    mutate(type = 'Epigenome for H3K9ac',
           Loading = c(sign(crossprod(data_rnaseq_1$Loading, data_H3K9ac_1$Loading))) * Loading,
           Dim = 'Level 1')
data_H3K18ac_1 <- data_H3K18ac %>%
    select(Time, Loading = `1`) 
data_H3K18ac_1 <- data_H3K18ac_1 %>%
    mutate(type = 'Epigenome for H3K18ac',
           Loading = c(sign(crossprod(data_rnaseq_1$Loading, data_H3K18ac_1$Loading))) * Loading,
           Dim = 'Level 1')

# dim2
data_rnaseq_2 <- data_rnaseq %>%
    select(Time, Loading = `2`) %>%
    mutate(type = "Transcriptome from Zheng et al. 2014",
           Dim = 'Level 2')
data_H3K9ac_2 <- data_H3K9ac %>%
    select(Time, Loading = `2`) 
data_H3K9ac_2 <- data_H3K9ac_2 %>%
    mutate(type = 'Epigenome for H3K9ac',
           Loading = c(sign(crossprod(data_rnaseq_2$Loading, data_H3K9ac_2$Loading))) * Loading,
           Dim = 'Level 2')
data_H3K18ac_2 <- data_H3K18ac %>%
    select(Time, Loading = `2`) 
data_H3K18ac_2 <- data_H3K18ac_2 %>%
    mutate(type = 'Epigenome for H3K18ac',
           Loading = c(sign(crossprod(data_rnaseq_2$Loading, data_H3K18ac_2$Loading))) * Loading,
           Dim = 'Level 2')

data_res <- data_rnaseq_1 %>%
    bind_rows(data_H3K9ac_1) %>%
    bind_rows(data_H3K18ac_1) %>% 
    bind_rows(data_rnaseq_2) %>%
    bind_rows(data_H3K9ac_2) %>%
    bind_rows(data_H3K18ac_2) %>% 
    mutate(data = factor(type, levels=c("Transcriptome from Zheng et al. 2014", "Epigenome for H3K9ac", "Epigenome for H3K18ac")))

scale_factor <- (max(data_res$Loading) - min(data_res$Loading)) / 40
data_oxygen <- read_csv("oxygen_ref_YMC2014.csv") %>%
    mutate(dO2_scaled = dO2 * scale_factor - 1)

# 850*500
# ggplot(data_res) +
#     facet_wrap(~Dim, scales = "free_y", ncol = 1) +
#     geom_line(data = data_oxygen,
#               aes(x = Time, y = dO2_scaled, color = 'oxygen'),
#               color = "grey"
#               ) +
#     geom_line(aes(x = Time, y = Loading, color = data), linewidth = 0.5) +
#     geom_point(aes(x = Time, y = Loading, shape = data, color = data)) +
#     scale_y_continuous(sec.axis = sec_axis( ~ (. + 0.5) / scale_factor, name = "dO2(%)")) +
#     labs(x = "Time(h)", y = "Loading", color = "Data", shape = "Data") +
#     theme_bw() +
#     theme(
#         axis.text.x = element_text(size = 10, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         axis.title.x = element_text(size = 12, color = "black"),
#         axis.title.y = element_text(size = 12, color = "black"),
#         strip.text = element_text(size = 12),
#         strip.background = element_rect(fill = "white", colour = "black"),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.key.size = unit(0.3, "inches")
#         )

colors <- c(
    'Dissolved oxygen' = 'grey',
    'Transcriptome from Zheng et al. 2014' = '#FF7F50',
    'Epigenome for H3K9ac' = '#BC50F7',
    'Epigenome for H3K18ac' = '#32CD32'
)
shapes <- c(
    'Dissolved oxygen' = NA,
    'Transcriptome from Zheng et al. 2014' = 16,
    'Epigenome for H3K9ac' = 17,
    'Epigenome for H3K18ac' = 18
) 
breakss <- c(
    'Dissolved oxygen',
    'Epigenome for H3K9ac',
    'Transcriptome from Zheng et al. 2014',
    'Epigenome for H3K18ac'
) 

ggplot(data_res) +
    facet_wrap(~Dim, scales = "free_y", ncol = 1) +
    geom_line(data = data_oxygen, aes(x = Time, y = dO2_scaled, color = 'Dissolved oxygen'), linewidth = 0.25) +
    geom_line(aes(x = Time, y = Loading, color = data), linewidth = 0.35) +
    geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.2, color = 'grey45') +
    scale_y_continuous(sec.axis = sec_axis( ~ (. + 0.7) / scale_factor, name = expression(dO[2] * "(%)"))) +
    scale_color_manual(values = colors, breaks = breakss) +
    labs(x = "Time(h)", y = "Loading", color = "Data") +
    xlim(0.3, 4.05) +
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
ggsave("sample_loading_histone&YMC2014.pdf", width = 9, height = 8.5, units = "cm" ) 

