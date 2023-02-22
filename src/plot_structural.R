library(ggplot2)
library(tidyverse)
library(ggtext)


args <- commandArgs(T)
assoc = args[1]

df_struc <- read.csv(
  "/home/joel/Documents/swedegene/WGS_pipeline/master/structural_analysis/data/gene_struct.csv", sep = "\t", as.is=T, check.names = F, stringsAsFactors = F)
#assocs <- read.table("/home/joel/Documents/swedegene/WGS_pipeline/master/structural_analysis/data/assoc_comparisons")$V1
df_pheno <- read.csv("/home/joel/Documents/swedegene/WGS_pipeline/master/structural_analysis/data/phenofile_new", sep ="\t", as.is=T, check.names = F, stringsAsFactors = F)
df_struc <- t(df_struc)
colnames(df_struc) <- df_struc["ID", ]
rownames(df_pheno) <- df_pheno$IID


DF_all <- read.csv(
  "/home/joel/Documents/swedegene/WGS_pipeline/master/structural_analysis/Results/DF_all.tsv", sep = "\t", as.is=T, check.names = F, stringsAsFactors = F)
DF_all$CHROM <- DF_all[, "#CHROM"]
DF_all[, "#CHROM"] <- NULL
DF_all %>% filter(!status == "Low variants") -> DF_all

DF_all$gene <- gsub(" ", "", DF_all$gene)

cns <- DF_all %>% filter(test == assoc)

sig =-log10(0.05/nrow(cns))

cns <- cns %>% filter(status =="Success")


data_cum <- cns %>%
  group_by(CHROM) %>%
  summarize(max_pos = max(pos)) %>%
  mutate(add_pos = lag(cumsum(as.numeric(max_pos)), default = 0)) %>%
  select(CHROM, add_pos)

cns <- cns %>%
  inner_join(data_cum, by ="CHROM") %>%
  mutate(pos_cum = pos + add_pos)

axis_set <- cns %>%
  group_by(CHROM) %>%
  summarize(center = mean(pos_cum))

ylim <- cns %>%
  filter(p == min(p)) %>%
  mutate(ylim=abs(floor(log10(p)))+2) %>%
  pull(ylim)


lab = paste0("Adjusted \u237a: e-", round(sig, digits = 2))

cns_plot <- ggplot(
    cns, aes(x=pos_cum, y =-log10(p),
             color=as.factor(CHROM), size = -log10(p))) +
    geom_point(alpha= 0.75) +
    geom_hline(aes(yintercept=sig), linetype="dashed") +
    annotate(
      "text",
      label=lab ,
      x=axis_set$center[3], y=sig+0.2) +
    scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0,ylim)) +
    scale_color_manual(values = rep(c("#276fbf", "#183059"), length(axis_set$CHROM))) +
    scale_size_continuous(range = c(0.5, 3)) +
    labs(x = NULL,
         y="-log<sub>10</sub>(p)") +
    theme_bw() +
    theme(
      legend.position ="none",
      panel.border = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )


mes <- 178/2
ggsave(args[2], cns_plot, width = mes*((1+sqrt(5))/2), height = mes, units="mm", dpi = 400)

