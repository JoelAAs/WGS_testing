library(data.table)
library(ggplot2)
library(ggtext)
library(tidyverse)

## https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

args <- commandArgs(T)

gwas_data <- fread(
  "run_folder/exon_distance/concated/Cns.Toxicity~rest_annotated.logistic_distance",
  sep="\t",
) 
selected_genes <- read.csv("local_data/selected_genes/CNS.toxicitet~rest_excluded.tsv", stringsAsFactors = F)

gwas_data <- gwas_data %>% filter(ClosestElement != "None")
gwas_data$gene <- gwas_data$ClosestElement %>% sapply(function(x) str_split(x, "-")[[1]][2])

gwas_data %>% filter(gene %in% selected_genes$Gene) -> gwas_data

gwas_data$CHROM <- as.numeric(recode(gwas_data$CHROM, X = "23"))

gwas_data %>% filter(!is.na(P)) -> gwas_data
gwas_data %>% filter(Distance < 15000) -> gwas_data

data_cum <- gwas_data %>% 
  group_by(gene) %>% 
  summarize(max_pos = max(POS)) %>% 
  mutate(add_pos = lag(cumsum(as.numeric(max_pos)), default = 0)) %>% 
  select(gene, add_pos)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "gene") %>% 
  mutate(pos_cum = POS + add_pos)


axis_set <- gwas_data %>% 
  group_by(gene) %>% 
  summarize(center = mean(pos_cum))

ylim <- gwas_data %>% 
  filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)

sig <- -log10(0.05/nrow(gwas_data))

lab_line <- paste0("Adjusted \u03b1: ", round(sig, 2))

manhattan_plot <- ggplot(
  gwas_data, aes(x=pos_cum, y = -log10(P),
                 color = as_factor(gene), size = -log10(P))) +
  geom_hline(yintercept = sig, color = "grey40", linetype = "dashed") +
  annotate("text", label=lab_line, x=axis_set$center[4], y = sig + 0.2) +
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
size=178/2
ggsave("paper-prep/cns-tox/candidate_genes.png", manhattan_plot, width = size*((1+sqrt(5))/2), height = size, units="mm", dpi = 400)
