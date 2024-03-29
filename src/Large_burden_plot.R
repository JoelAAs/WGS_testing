library(data.table)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(qvalue)

df_genes <- read.csv("local_data/burden_results/burden_vf_threshold_full_results.tsv", sep = "\t", as.is=T)
df_genes_unfiltered <- df_genes
df_genes %>% filter(p != 1 & p != 0) -> df_genes

genes <- fread(
  "/home/joel/Documents/Resourses/ensembl/Homo_sapiens.GRCh37.87.chr.gtf", 
  sep = "\t", skip=5) %>%  filter(V3 == "gene") 
genes$gene <- sapply(genes$V9, function(x)
  gsub('[\ ]*"', "", str_split(str_split(x, ";")[[1]][3], "gene_name")[[1]][2]))
genes <- genes %>%  filter(gene %in% df_genes$gene)
genes <- genes[!duplicated(genes$gene), ]
genes <- genes[, c("V1", "V4", "V5", "gene")]
colnames(genes) <- c("CHROM", "POS", "END", "gene")
df_genes <- inner_join(df_genes, genes, by = "gene")


df_genes$CHROM <- as.numeric(recode(df_genes$CHROM, X = "23", Y="23"))
df_genes$POS

data_cum <- df_genes %>% 
  group_by(CHROM) %>% 
  summarize(max_pos = max(POS)) %>% 
  mutate(add_pos = lag(cumsum(as.numeric(max_pos)), default = 0)) %>% 
  select(CHROM, add_pos)

df_genes <- df_genes %>% 
  inner_join(data_cum, by = "CHROM") %>% 
  mutate(pos_cum = POS + add_pos)


axis_set <- df_genes %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(pos_cum))

ylim <- df_genes %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

sig <- 0.05/length(genes$gene)

df_genes$fdr <- p.adjust(df_genes$p, "fdr")

manhattan_plot2 <- ggplot(
  df_genes, aes(x=pos_cum, y = fdr,
                 color = as_factor(CHROM), size = -log10(p))) +
  geom_hline(yintercept = 0.05, color = "grey40", linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHROM)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL,
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
#    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  ) + facet_wrap(test ~ .)

hist(df_genes_unfiltered$p)
df_genes$fdr <- p.adjust(df_genes$p, "fdr")
df_genes$q <- qvalue(df_genes$p)
ggplot(df_genes, aes(x=p, fill = test)) + geom_histogram(bins = 40) + facet_wrap(test ~.)

df_genes %>%filter(p == 0)
manhattan_plot2

hist(df_genes$p)


