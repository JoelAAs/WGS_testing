library(tidyverse)
library(data.table)

args <- commandArgs(T)

gwas_data <- fread(
  args[1],
  sep="\t",
  )

gwas_data$test <- args[3]
gwas_data %>% filter(!is.na(P)) -> gwas_data
gwas_data %>% filter(P <= 5*10^-7) -> gwas_data

gwas_data %>%  {write.table(file = args[2], ., sep = "\t", row.names = F)}