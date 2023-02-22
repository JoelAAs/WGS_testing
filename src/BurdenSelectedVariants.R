library(data.table)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(qvalue)
library(ggrepel)


tests_genes <- read.delim("local_data/burden_results/burden_selected/tests_genes.tsv", sep="\t")
pheno_table_count <-  read.delim("run_folder/phenofile_counts.table", sep = "\t")
pheno_table_count$test_name <- sapply(pheno_table_count$assoc, function(x) gsub("~rest", "", x))
burden_tests <- read.csv("local_data/burden_results/burden_dynamic_threshold_full_results.tsv", sep = "\t")

for (test in unique(tests_genes$test)){
  n_cases <- pheno_table_count %>% filter(test_name == test) %>%  {.$n_cases} 
  VF_t <- (1/n_cases)^(1/2)
  VF_column <- paste0(test, "~rest_VF")
  
  test_name <- paste0(test, "~rest") 
  n_tests <- burden_tests %>% filter(test == test_name) %>% nrow()
  
  print(paste("Testing:", test))
  selected_vars <- tests_genes[tests_genes$test == test, "gene"] %>% 
    lapply(function(x) read.csv(
      paste0("local_data/burden_results/burden_selected/", x, "_var_final.tsv"), sep=";", as.is=T, check.names = F)) %>% bind_rows()
  
  selected_vars$plinkID <- apply(
    selected_vars, 1, function(x)
      gsub(" ", "", paste0(x[3], "_", x[4], ":", x[6], ":", x[7])))
  
  plink_res <- fread(
      paste0("run_folder/no_pc/", test, "~rest.no_pc.glm.logistic.hybrid"),
      sep = "\t") %>% filter(`#CHROM` %in% unique(selected_vars$CHROM)) %>% filter(ID %in% selected_vars$plinkID) 

  plink_res$plinkID <- plink_res$ID
  plink_res %>%  filter(TEST == "ADD") -> plink_res
  selected_vars <- merge(selected_vars, plink_res, by="plinkID")
  selected_vars %>%  filter(VF != 0 & Weight != 0) -> current_all
  
  current_all[, VF_column] <- as.numeric(gsub(",", ".", current_all[, VF_column]))
  current_all$gene <- sapply(current_all$ClosestElement, function(x) strsplit(x, "-")[[1]][2])
  current_all <- current_all[current_all[, VF_column] <= VF_t, ]
  
  sig_table <- data.frame(gene=unique(current_all$gene))
  sig_table$sig <-  sapply(sig_table$gene, function(x) -log10(0.05/(nrow(current_all[current_all$gene == x,]) + n_tests)))
  sig_table$P <-sig_table$sig
  sig_table$label <- sapply(sig_table$sig, function(x) paste0("Adjusted \u03b1: e-", round(x, 2)))
  sig_table$x <- current_all %>% 
    group_by(gene) %>% 
    summarize(x = floor((max(POS.x) - min(POS.x))/5)) %>% 
    {.$x}
  
  g <- ggplot(
    current_all, aes(x=POS.x, y = -log10(P),
                  color =-log10(P) , size = -log10(P))) +
    geom_point(alpha = 0.75) +
    geom_hline(data=sig_table, aes(yintercept = sig), color = "grey40", linetype = "dashed") +
    geom_text(data=sig_table,
                     aes(label = label, 
                         y = sig
                     ), 
              x = -Inf, hjust = -0.1, nudge_y = 0.4,  size = 4) +
    geom_label_repel(
      aes(label =
            ifelse(-log10(P)> 3, snp138, "")
      ), point.padding = 0.5, force=10, nudge_y =-0.6, size = 3) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Genomic position",
         y = "-log<sub>10</sub>(p)") + 
    theme_bw() +
    theme( 
      legend.position = "none",
      #panel.border = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(size = 8, vjust = 0.5)
    ) + facet_wrap(gene ~ ., scales="free")
  
  size = 178
  write.table(current_all, paste0("Resultat/burden/",test,".tsv"), sep = "\t")
  ggsave(paste0("Resultat/burden/", test, "_selected_variants.png"), g, width = size*((1+sqrt(5))/2), height = size, units="mm", dpi = 400)
}
