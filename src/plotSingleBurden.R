  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(tidyverse)
  library(qvalue)
  library(ggrepel)
  
  args <- commandArgs(T)
  
  mytest <- args[2]
  
  df_genes <- read.csv("local_data/burden_results/burden_dynamic_threshold_full_results.tsv", sep = "\t", as.is=T)
  df_genes %>% filter(test == mytest) -> df_genes

  sig <- 0.95#-log10(0.05/nrow(df_genes))
  df_genes$qvalue <- qvalue(df_genes$p)$qvalues
  
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
  
  #ylim <- df_genes %>% 
  #  filter(p == min(p)) %>% 
  #  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  #  pull(ylim)
  

  
  lab_line=paste0("Qvalue : 0.05")
  
  g <- ggplot(
    df_genes, aes(x=pos_cum, y = 1- qvalue,
                  color = as_factor(CHROM), size = -log10(p))) +
    geom_hline(yintercept = sig, color = "grey40", linetype = "dashed") +
    annotate("text", label=lab_line, x=axis_set$center[2], y = sig + 0.2) +
    geom_label_repel(
      aes(label =
            ifelse(qvalue < 0.15, gene, "")
          ), point.padding = 0.5, force=10, nudge_y =0.4) +
    geom_point(alpha = 0.75) +
    scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHROM)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL,
         y = "1 - qvalue") +
    theme_bw() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )
  
  size = 178/2
  ggsave(paste0("run_folder/Plots/Burden/", mytest, "_all.png"), g, width = size*((1+sqrt(5))/2), height = size, units="mm", dpi = 400)
  
  
  
    