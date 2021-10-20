library(ggplot2)
library(qvalue)
library(ggtext)
library(tidyverse)

df_burden <- read.csv("data/burden_results/burden_allow_all_results.tsv", sep = "\t", stringsAsFactors = F)
q <- qvalue(p=df_burden$p)


df_burden$fdr <- p.adjust(df_burden$p, method="fdr")
df_burden <- df_burden[order(df_burden$p, df_burden$gene, decreasing = T),]
df_burden$gene = fct_inorder(df_burden$gene)
g <- ggplot(df_burden, aes(x=gene, y=fdr, size = n, color = gene)) + 
  geom_point() + 
  geom_hline(aes(yintercept=0.05),color = "grey40", linetype = "dashed") +
  lims(y=c(0,1)) + 
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(df_burden$gene)))) + 
  coord_flip() +
  theme_bw() +
  labs(
    x="Gene",
    y= "FDR") +
  theme( 
    legend.position = "none",
    #panel.border = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.y = element_text(size = 8, vjust = 0.5)
  )

  ggsave("run_folder/Plots/CNS-tox_selected_burden.png",
         g, width = 86 , height=86, unit="mm", dpi=400)
