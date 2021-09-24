library(data.table)
library(ggplot2)
library(ggtext)
library(tidyverse)

## https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

args <- commandArgs(T)

gwas_data <- fread(
  args[2],
  sep="\t",
  ) 

unique(gwas_data$CHROM)

gwas_data$CHROM <- as.numeric(recode(gwas_data$CHROM, X = "23"))

gwas_data %>% filter(!is.na(P)) -> gwas_data
gwas_data %>% filter(Distance < 4000) -> gwas_data

data_cum <- gwas_data %>% 
  group_by(CHROM) %>% 
  summarize(max_pos = max(POS)) %>% 
  mutate(add_pos = lag(cumsum(as.numeric(max_pos)), default = 0)) %>% 
  select(CHROM, add_pos)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "CHROM") %>% 
  mutate(pos_cum = POS + add_pos)


axis_set <- gwas_data %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(pos_cum))

ylim <- gwas_data %>% 
  filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)

sig <- 5e-8

manhattan_plot <- ggplot(
  gwas_data, aes(x=pos_cum, y = -log10(P),
                 color = as_factor(CHROM), size = -log10(P))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHROM)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL,
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

ggsave(args[3], manhattan_plot, width = 7, height = 5)
