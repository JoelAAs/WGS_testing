library(ggplot2)
library(plotly)
library(ggtext)
library(tidyverse)
library(ggrepel)


args <- commandArgs(T)

selected_subset = args[2]

hla_assoc_df <- read.csv(args[1],
  sep="\t", as.is = T, check.names = F)
hla_assoc_df$PVALUE <- as.numeric(hla_assoc_df$PVALUE)
hla_assoc_df <- hla_assoc_df[!is.na(hla_assoc_df$PVALUE),]

#hla_assoc_df_sub <- hla_assoc_df[hla_assoc_df$TEST == selected_subset, ]

# Associations
#gg1 <- ggplot(hla_assoc_df, aes(x = ALLELE, text=ALLELE, label=TEST, y=-log10(PVALUE), color=HLA_GROUP)) +
#  geom_point() +
#  geom_hline(yintercept = -log10(0.05/nrow(hla_assoc_df))) + # Bon ferroni on everything
#  ylim(0, 60) +
#  facet_grid(GENETIC_MODEL ~ .)
#
#ggplotly(gg1)


### Good CNS

hla_assoc_df %>%  filter(TEST == selected_subset & GENETIC_MODEL=="dominant") -> test
test %>%  filter(PVALUE != 1) -> test
test %>% filter(PVALUE != 0) -> test
test <- test[order(test$ALLELE),]
test$POS <- 1:nrow(test)

data_cum <- test %>%
  group_by(HLA_GROUP) %>%
  summarize(max_pos = max(POS)) %>%
  mutate(add_pos = lag(cumsum(as.numeric(max_pos)), default = 0)) %>%
  select(HLA_GROUP, add_pos)

test <- test %>%
  inner_join(data_cum, by = "HLA_GROUP") %>%
  mutate(pos_cum = POS + add_pos)


axis_set <- test %>%
  group_by(HLA_GROUP) %>%
  summarize(center = mean(pos_cum))

ylim <- test %>%
  filter(PVALUE == min(PVALUE)) %>%
  mutate(ylim = abs(floor(log10(PVALUE))) + 1) %>%
  pull(ylim)


bon_sig <- 0.05/nrow(test)
sig_label = paste0("Adjusted \u03b1: e-", round(-log10(bon_sig), 2))

if (ylim < -log10(bon_sig)){
  ylim <- round(-log10(bon_sig) + 1)
}

g <- ggplot(test, aes(x = pos_cum, y=-log10(PVALUE), color=as_factor(HLA_GROUP))) +
  geom_point() +
  geom_hline(yintercept = -log10(bon_sig),  color = "grey40", linetype = "dashed") +
  annotate("text",
           y=-log10(bon_sig) + 0.2,
           x= axis_set$center[5],
           label = sig_label) +
  theme_bw() +
  geom_label_repel(
    aes(label =
          ifelse(-log10(PVALUE)> 3, ALLELE, "")
    ), point.padding = 0.5, force=10, nudge_y =-0.2, size = 3) +
  scale_x_continuous(label = axis_set$HLA_GROUP, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$HLA_GROUP)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL,
       y = "-log<sub>10</sub>(p)") +
  theme(
    legend.position = "None",
    panel.border = element_blank(),

    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

size=178/2
ggsave(args[3], g, width = size*((1+sqrt(5))/2), height = size, units="mm", dpi = 400)
