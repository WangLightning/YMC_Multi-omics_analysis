rm(list = ls())
suppressMessages({
    library(tidyverse)
})

dir_result <- "../../results"

entropy <- function(x){
    s <- sum(x)
    p <- x/s 
    e <- -p*log(p)
    sum(e[!is.na(e)])
}

# entropy  -----------------------------------------------------------------

file_data <- file.path(dir_result, "temp", "H3K9ac", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
entropy(data$ratio) # 1.601024

file_data <- file.path(dir_result, "temp", "H3K18ac", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
entropy(data$ratio) # 1.116343

file_data <- file.path(dir_result, "temp", "concatenate", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
entropy(data$ratio) # 1.792348

file_data <- file.path(dir_result, "temp", "Metabolite_interpolated_lc", "loadings/s_values.csv")
data <- read_csv(file_data, show_col_types = FALSE) %>% 
    select(-value) %>% 
    filter(dim >= 1) 
entropy(data$ratio) # 2.331486


# plot --------------------------------------------------------------------

rm(list = ls())


pdf("Svalue_entropy.pdf", width = 12 / 1.48, height = 5.5 / 1.48)  

# x <- 3:1
# y <- c(86.0649, 75.55277, 50.59236)
# y2 <- c(67.8563, 49.4773, 34.1819)

x <- 4:1
y <- c(86.0649, 70.24324, 75.55277, 50.59236)
y2 <- c(67.8563, 54.2053859133921, 49.4773, 34.1819)

par(mar = c(5, 12, 2, 2))  
# plot(x, y, type = "n", ylim = c(0.5,3.5), xlim = c(-10,110), axes = F, ann = F)
plot(x, y, type = "n", ylim = c(0.5,5), xlim = c(-5,125), axes = F, ann = F)


axis(1, at = seq(0, 100, 20))
mtext("Contribution proportion (%)", side = 1, line = 2.5)
# axis(2, at = 3:1, labels = c("Epigenome of H3K18ac", "Concatenated transcriptome", "Metabolome"), las = 1)
axis(2, at = 4:1, labels = c("Epigenome for H3K18ac", "Epigenome for H3K9ac", "Transcriptome", "Metabolome"), las = 1)
# axis(4, at = 4:1, labels = c("1.12", "1.60", "1.79", "2.33"), las = 1)
# mtext("Entropy", side = 4, line = 2.5)
box()

text(114, 4.6, labels = "Eigen-entropy")
text(115, 1:4, labels = c("2.331", "1.792", "1.601", "1.116"))


segments(0, x, 100, x, lwd = 20)
segments(0, x, 100, x, lwd = 16, col = "white")
segments(0, x, y, x, lwd = 16, col = "royalblue")
segments(0, x, y2, x, lwd = 16, col = "#F08080")

legend("top",         
       legend = c("Level 1", "Level 2"),   
       col = c("#F08080", "royalblue"), 
       lwd = 16, 
       horiz = TRUE, 
       bty = "n")         


dev.off()
