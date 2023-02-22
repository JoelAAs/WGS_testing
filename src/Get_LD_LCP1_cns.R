library(data.table)
library(tidyverse)

phenofile <- read.csv(
  "local_data/phenofile",
  sep = "\t",
  stringsAsFactors = F,
  check.names = F)

excluded <- c("RI-1933-4606",
              "RI-1933-1457",
              "RI-1933-1334",
              "RI-1933-1315",
              "RI-1933-1096",
              "RB-1745-651",
              "RB-1745-5754",
              "RB-1745-2961",
              "RB-1745-2816",
              "RI-1933-1728",
              "RI-1933-1594",
              "RI-1933-1521",
              "RB-1745-1009",
              "RB-1745-701",
              "RI-1933-697",
              "RB-1745-693",
              "RB-1745-580",
              "RB-1745-575",
              "RI-1933-572",
              "RB-1745-577")


samples <- phenofile[phenofile[ ,"Cns.Toxicity~rest"] == 2, "IID"]
controls <- phenofile[phenofile[ ,"Cns.Toxicity~rest"] == 1, "IID"] %>% 
  {.[!grepl("SweGen",.)]}  
samples <- samples[!samples %in% excluded] %>% unique()
#samples <- c(samples, controls)

lcp1 <- fread("local_data/data_cns/lcp1_subset.vcf", sep = "\t")
lcp1$ID <- paste0(lcp1$CHROM,"_", paste(lcp1$POS, lcp1$REF, lcp1$ALT, sep = ":")) 
lcp1 <- as.data.frame(lcp1)
lcp1_cases <- lcp1[,  c("ID", samples)]

for (sam in samples) {
  lcp1_cases[, sam] <-  sapply(lcp1_cases[, sam], function(x) sum(as.numeric(strsplit(strsplit(x, ":")[[1]][1], "/")[[1]])))
}

target <- "13_46718722:T:C"
lcp1_cases$p <- rowSums(lcp1_cases[, 2:ncol(lcp1_cases)],na.rm = T)/(nrow(lcp1_cases)*2)
lcp1_cases %>%  filter(ID == target) %>%  {.$p} -> pb
lcp1_cases_t <- as.data.frame(t(lcp1_cases))
colnames(lcp1_cases_t) <- lcp1_cases_t[1,]
lcp1_cases_t <- lcp1_cases_t[2:nrow(lcp1_cases_t),]
lcp1_cases_t <- sapply(lcp1_cases_t, as.numeric) %>%  as.data.frame()
pab = rep(0, length(colnames(lcp1_cases_t)))
j = 0 
for (variant in colnames(lcp1_cases_t)) {
  j = j + 1
  pab[j] = sum(apply(lcp1_cases_t[1:(nrow(lcp1_cases_t)-1), c(variant, target)], 1,FUN = function(x) min(x, na.rm=T)))
}

lcp1_cases$pab <- pab/(nrow(lcp1_cases)*2)
lcp1_cases$D   <- lcp1_cases$pab - lcp1_cases$p*pb
lcp1_cases$D_prim <- apply(lcp1_cases[, c("D", "pab", "p")], 1, FUN = function(x){
  if (x["D"] <0) {
    return(x["D"]/max(c(-x["pab"], -(1-x["p"])*(1-pb))))
  } else {
    return(x["D"]/max(c(x["p"]*(1-pb), pb*(1-x["p"]))))
  }
})

lcp1_cases$r2 <- lcp1_cases$D^2/(lcp1_cases$p*(1-lcp1_cases$p)*pb*(1-pb))

lcp1_cases %>%  filter(!is.infinite(D_prim) &  !is.infinite(r2)) -> lcp1_cases
top_ld <- lcp1_cases %>% filter(D_prim > 0.95 | r2 > 0.95) %>%  {.[, c("ID", "pab", "p", "D_prim", "r2")]}
view(top_ld)

write.table(top_ld, "top_ld_vars.csv")

