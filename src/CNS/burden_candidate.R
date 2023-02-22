library(data.table)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(qvalue)
library(ggrepel)

args <- commandArgs(T)

mytest <-"Cns.Toxicity~rest"

df_genes <- read.csv("local_data/burden_results/burden_vf_threshold_full_results.tsv", sep = "\t", as.is=T)
df_genes %>% filter(p != 1 & p != 0) -> df_genes
df_genes %>% filter(test == mytest) -> df_genes

selected_genes <- read.csv("local_data/selected_genes/CNS.toxicitet~rest_excluded.tsv", stringsAsFactors = F)


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

df_genes %>% filter(gene %in% selected_genes$Gene) -> df_genes

df_genes$CHROM <- as.numeric(recode(df_genes$CHROM, X = "23", Y="23"))

data_cum <- df_genes %>% 
  group_by(gene) %>% 
  summarize(max_pos = max(POS)) %>% 
  mutate(add_pos = lag(cumsum(as.numeric(max_pos)), default = 0)) %>% 
  select(gene, add_pos)

df_genes <- df_genes %>% 
  inner_join(data_cum, by = "gene") %>% 
  mutate(pos_cum = POS + add_pos)

axis_set <- df_genes %>% 
  group_by(gene) %>% 
  summarize(center = mean(pos_cum))

ylim <- df_genes %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 1) %>% 
  pull(ylim)

sig <- -log10(0.05/nrow(df_genes))


lab_line=paste0("Adjusted \u03b1: e-", round(sig, 2))

g <- ggplot(
  df_genes, aes(x=pos_cum, y = -log10(p),
                color = as_factor(gene), size = -log10(p))) +
  geom_hline(yintercept = sig, color = "grey40", linetype = "dashed") +
  annotate("text", label=lab_line, x=axis_set$center[4], y = sig + 0.2) +
  geom_label_repel(
    aes(label =
          ifelse(-log10(p)> 4, gene, "")
    ), point.padding = 0.5, force=10, nudge_y =0.4) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$gene, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$gene)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL,
       y = "-log<sub>10</sub>(p)") + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

size = 178/2
ggsave(paste0("paper-prep/cns-tox/", mytest, "_selected.png"), g, width = size*((1+sqrt(5))/2), height = size, units="mm", dpi = 400)



