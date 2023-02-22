	
library(data.table)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(qvalue)
library(ggrepel)

test = "Cns.Toxicity"
variants <- read.csv(paste0(
  "Resultat/burden/",test, ".tsv")
  , sep = "\t")

n_cases = 66
VF_t <- (1/n_cases)^(1/2)
VF_column <- paste0(test, "~rest_VF")

## OBS 2 är cases
case_freq <- fread(paste0("local_data/freq/",test,"~rest_2.afreq"), sep = "\t") 
control_freq <- fread(paste0("local_data/freq/",test,"~rest_1.afreq"), sep = "\t")

case <- merge(variants, case_freq[, c("ID", "ALT_FREQS", "OBS_CT")], by.x="plinkID", by.y="ID",suffixes = c("_test", "_case"), all.x=TRUE)
all <- merge(case, control_freq[, c("ID", "ALT_FREQS", "OBS_CT")], by.x="plinkID", by.y="ID", suffixes = c("_case", "_control"), all.x=TRUE)

keep <- all[, c("plinkID", "gene", "snp138", "P", "OR", "LOG.OR._SE", "ALT_FREQS_case", "ALT_FREQS_control")]
keep <- keep[order(keep$P),]

write.csv(keep[1:40,], "Resultat/CNS-tox/40_top_burden_variants.csv", row.names = FALSE)
