library(ggplot2)
library(tidyverse)
library(dplyr)
library(data.table)


args = commandArgs(T)

csv_assoc <- fread(
  args[2],
  sep="\t") 

#csv_assoc$Distance <- as.numeric(csv_assoc$Distance)
#csv_assoc %>% filter(Distance < 50) -> csv_assoc
if ("X" %in% csv_assoc$CHROM) {
  csv_assoc$CHROM <- as.numeric(recode(csv_assoc$CHROM, X = "23"))
} else {
  csv_assoc$CHROM <- as.numeric(csv_assoc$CHROM)
}

csv_assoc %>% arrange(CHROM, POS) -> csv_assoc
csv_assoc$order <- 1:nrow(csv_assoc)

    

g <- ggplot(csv_assoc, aes(x=POS, y=-log10(P))) +
  geom_point() +
  geom_hline(yintercept = 8) +
  geom_hline(yintercept = -log10(0.05/nrow(csv_assoc))) +
  ggtitle(args[4])

ggsave(filename = args[3], plot = g, width = 10)

