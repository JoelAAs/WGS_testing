library(ggplot2)
library(tidyverse)
library(dplyr)
library(data.table)

csv_assoc <- fread(
  "run_folder/selected_genes/variants/CNS.toxicitet~rest_excluded.logistic.hybrid",
  sep="\t") 

selected <-  fread(
  "run_folder/selected_genes/CNS.toxicitet~rest_excluded.tsv",
  sep = "\t")

csv_assoc %>%  filter(CHROM %in% unique(selected$CHROM)) -> csv_assoc
csv_assoc$CHROM <- as.numeric(csv_assoc$CHROM)
genes <- array("", nrow(csv_assoc))
for (i in 1:nrow(csv_assoc)){
  genes[i] <- filter(selected, CHROM == csv_assoc$CHROM[i] & START-20 <= csv_assoc$POS[i] & END+20 >= csv_assoc$POS[i]) %>%  {.$gene_name}
}
csv_assoc$gene <- genes

if ("X" %in% csv_assoc$CHROM) {
  csv_assoc$CHROM <- as.numeric(recode(csv_assoc$CHROM, X = "23"))
} else {
  csv_assoc$CHROM <- as.numeric(csv_assoc$CHROM)
}



csv_assoc %>% arrange(CHROM, POS) -> csv_assoc
csv_assoc$order <- 1:nrow(csv_assoc)
csv_assoc$logp <- -log10(csv_assoc$P)
mylow = "#1b8061"
myhigh = "#c44610"

for (mygene in unique(csv_assoc$gene)){
  g <- ggplot(csv_assoc %>%  filter(gene == mygene), aes(x=POS, y=-log10(P), colour = -log10(P))) +
    geom_point(aes()) +
    geom_hline(yintercept = 8) +
    geom_hline(yintercept = -log10(0.05/nrow(csv_assoc))) +
    ggtitle(mygene) + 
    theme_bw() +
    xlab(paste("Chromosome:", selected %>%  filter(gene_name == mygene) %>%  {.$CHROM[[1]]}, "Position")) +
    scale_colour_gradient(
      limits =c(0, 8),
      low = mylow,
      high = myhigh
    ) +
    theme(
      legend.position = "none"
      axis.
    )
  
  
  
  ggsave(filename = paste0("Resultat/CNS-tox/selected_genes/", mygene,".png"), plot = g, width = 5, height=5)
}
