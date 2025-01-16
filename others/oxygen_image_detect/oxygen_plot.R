suppressMessages({
    library(tidyverse)
})

# LC-MS -------------------------------------------------------------------

data_raw <- read_csv('LC-MS_oxygen.csv') %>%  
    filter(Time > 42.7) %>%
    mutate(n = row_number())

# the time of each sample corresponds to the row number
times <- c(
    6,
    16,
    28,
    37,
    46,
    56,
    68,
    77,
    87,
    100,
    110,
    121,
    126,
    137,
    146,
    156,
    167,
    176,
    184,
    194,
    204,
    216,
    227,
    238
)

data_point <- data_raw %>% 
    filter(n %in% times) %>% 
    mutate(t = 1:length(times))

p <- ggplot() +
    geom_line(data = data_raw, aes(x = Time, y = dO2), color = "#FF7F50") +
    geom_point(data = data_point, aes(x = Time, y = dO2)) +
    geom_text(data = data_point, aes(x = Time, y = dO2, label = t), vjust = 1.5) +
    labs(title = "Oxygen concentrations of metabolome LC-MS", x = "Time (h)", y = "dO2 (%)") +
    theme_light()
ggsave("Oxygen_lc.pdf", p, height = 3.5, width = 7)


# TOFMS -------------------------------------------------------------------

data_raw <- read_csv('TOFMS_oxygen.csv') %>%  
    filter(Time > 65.5) %>% 
    mutate(n = row_number())

# the time of each sample corresponds to the row number
times <- c(
    4,
    14,
    24,
    36,
    48,
    60,
    72,
    83,
    95,
    107,
    118,
    127,
    136,
    147,
    157,
    167,
    178,
    190,
    202,
    212,
    224,
    236,
    248,
    255
)

data_point <- data_raw %>% 
    filter(n %in% times) %>% 
    mutate(t = 1:length(times))

p <- ggplot() +
    geom_line(data = data_raw, aes(x = Time, y = dO2), color = "#FF7F50") +
    geom_point(data = data_point, aes(x = Time, y = dO2)) +
    geom_text(data = data_point, aes(x = Time, y = dO2, label = t), vjust = 1.5) +
    labs(title = "Oxygen concentrations of metabolome GC-TFOMS", x = "Time (h)", y = "dO2 (%)") +
    theme_light()
ggsave("Oxygen_gc.pdf", p, height = 3.5, width = 7)


# YMC2014 -------------------------------------------------------------------

data_raw <- read_csv('YMC2014_oxygen_rnaseq.csv') %>%  
    mutate(n = row_number())

# the time of each sample corresponds to the row number
times <- c(
    133,
    248,
    291,
    322,
    378,
    441,
    504,
    563,
    626,
    689,
    754,
    815,
    878,
    1079,
    1275,
    1479
)

data_point <- data_raw %>% 
    filter(n %in% times) %>% 
    mutate(t = 1:length(times))

p <- ggplot() +
    geom_line(data = data_raw, aes(x = Time, y = dO2), color = "#FF7F50") +
    geom_point(data = data_point, aes(x = Time, y = dO2)) +
    geom_text(data = data_point, aes(x = Time, y = dO2, label = t), vjust = 1.5) +
    labs(title = "Oxygen concentrations of transcriptome of Zheng et al. 2014", x = "Time (h)", y = "dO2 (%)") +
    theme_light()
ggsave("Oxygen_YMC2014.pdf", p, height = 3.5, width = 6)

# histone -------------------------------------------------------------------

data_raw <- read_csv('YMC2014_oxygen_chipseq.csv') %>%  
    mutate(n = row_number())

# the time of each sample corresponds to the row number
times <- c(
    191,
    339,
    406,
    443,
    503,
    566,
    629,
    688,
    752,
    815,
    889,
    948,
    1108,
    1275,
    1428,
    1584
)

data_point <- data_raw %>% 
    filter(n %in% times) %>% 
    mutate(t = 1:length(times))

p <- ggplot() +
    geom_line(data = data_raw, aes(x = Time, y = dO2), color = "#FF7F50") +
    geom_point(data = data_point, aes(x = Time, y = dO2)) +
    geom_text(data = data_point, aes(x = Time, y = dO2, label = t), vjust = 1.5) +
    labs(title = "Oxygen concentrations of epigenome", x = "Time (h)", y = "dO2 (%)") +
    theme_light()
ggsave("Oxygen_histone.pdf", p, height = 3.5, width = 6)


# YMC2005 -------------------------------------------------------------------

data_raw <- read_csv('YMC2005_oxygen.csv') %>%  
    mutate(Time = Time * 5/12, n = row_number())

# the time of each sample corresponds to the row number
times <- c(
    1,
    30,
    60,
    89,
    118,
    147,
    177,
    206,
    235,
    264,
    294,
    323,
    352,
    381,
    411,
    440,
    469,
    498,
    528,
    557,
    586,
    615,
    645,
    674,
    703,
    732,
    762,
    791,
    820,
    849,
    879,
    908,
    937,
    966,
    996,
    1025
)

data_point <- data_raw %>% 
    filter(n %in% times) %>% 
    mutate(t = 1:length(times))

p <- ggplot() +
    geom_line(data = data_raw, aes(x = Time, y = dO2), color = "#FF7F50") +
    geom_point(data = data_point, aes(x = Time, y = dO2)) +
    geom_text(data = data_point, aes(x = Time, y = dO2, label = t), vjust = 1.5) +
    labs(title = "Oxygen concentrations of transcriptome of Tu et al. 2005", x = "Time (h)", y = "dO2 (%)") +
    theme_light()
ggsave("Oxygen_YMC2005.pdf", p, height = 3.5, width = 10)



