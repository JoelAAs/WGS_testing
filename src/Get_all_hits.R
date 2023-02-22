library(tidyverse)
library(data.table)

tests = c(
  "Hypersensitivity.All~rest",
  "Hypersensitivity.Narrow~rest",
  "Fototoxicitet~rest",
  "SJS.TEN~rest",
  "Cytopenia.All~rest",
  "White.Blood.Cytopenia~rest",
  "Trombocytopeni~rest",
  "Blödning~rest",
  "CNS.toxicitet~rest",
  "Hypo.Siadh~rest",
  "Levertoxicitet~rest",
  "Metabol.rubbning~rest",
  "Njurtoxicitet~rest",
  "Patologisk.fraktur~rest",
  "Pankreatit~rest",
  "Narkolepsi~rest"
)



variants_list <- list()
i=1
for (test in tests) {
  variants <- paste0("run_folder/exon_distance/concated/", test, "_annotated.logistic_distance")
  variants %>% 
    fread(sep="\t", stringsAsFactors = F) %>%  
    filter(P < 10e-8) -> variants
  variants$group <- test
  variants_list[[i]] <- variants
  i = i + 1
}

burden <- fread("data/burden_results/burden_vf_threshold_full_results.tsv", sep = "\t", stringsAsFactors = F) %>% 
  filter(p != 1 & p != 0) %>% 
  group_by(test) %>%  
  mutate(fdr = p.adjust(p, "fdr")) %>% 
  filter(fdr < 0.15) 


all <- do.call(rbind, variants_list)
