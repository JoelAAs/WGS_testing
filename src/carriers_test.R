library(tidyverse)
library(VennDiagram)

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
              "RB-1745-577",
              "RB-1745-5400")


drug_df <- read.csv("../Phenotype_mapping/data/all_medication.csv",
                    stringsAsFactors = F,
                    sep = "\t")

adr_df <- read.csv("../Phenotype_mapping/data/adr_diagnosis.csv",
                   stringsAsFactors = F,
                   sep = "\t")

variants <- read.csv("Resultat/burden/Cns.Toxicity.tsv",
                     stringsAsFactors = F,
                     sep = "\t")

pheno <- read.csv("local_data/phenofile",
                  stringsAsFactors = F,
                  sep = "\t")

cases <- pheno[pheno$Cns.Toxicity.rest == 2, "IID"]  %>% {.[!. %in% excluded]}
controls <- pheno[pheno$Cns.Toxicity.rest == 1, "IID"] %>% 
  {.[!grepl("SweGen",.)]} %>% 
  {.[!. %in% excluded]}


gmat_retsat <- read.csv("local_data/burden_results/burden_selected/RETSAT_gmat.tsv",
                      stringsAsFactors = F,
                      sep = "\t", check.names = F)
gmat_lcp1 <- read.csv("local_data/burden_results/burden_selected/LCP1_gmat.tsv",
                      stringsAsFactors = F,
                      sep = "\t", check.names = F)
gmat_SFMBT2 <- read.csv("local_data/burden_results/burden_selected/SFMBT2_gmat.tsv",
                        stringsAsFactors = F,
                        sep = "\t", check.names = F)
  

drug_df[drug_df$patient == "RB-1745-5400", "patient"] = "RB-1745-merged"
drug_df %>% 
  filter(patient %in% cases) %>% 
  filter(suspect_med == "True") -> case_med

adr_df[adr_df$patient == "RB-1745-5400", "patient"] = "RB-1745-merged"
adr_df %>% 
  filter(patient %in% cases) -> case_adr

pheno_df <- merge(case_med, case_adr, by="patient")


retsat_burden_rsid <- c("rs143283662", "rs139043592")
lcp1_burden_rsid <- c("rs6561297", "rs10492451")

get_geno <- function(gmat, variants, rsids, ids){
  gmat %>% 
    filter(PK %in% variants[variants$snp138 %in% rsids, "PK"]) %>% 
    {t(.[, c("PK", ids)])} -> out
  
  colnames(out) <- rsids
  out <- as.data.frame(out)
  out$patient <- row.names(out)
  out %>% 
    filter(patient != "PK") %>%
    return()  
}

get_geno(gmat_retsat, variants, retsat_burden_rsid, cases) -> retsat_cases_geno
get_geno(gmat_retsat, variants, retsat_burden_rsid, controls) -> retsat_controls_geno

get_geno(gmat_lcp1, variants, lcp1_burden_rsid, cases) -> lcp1_cases_geno
get_geno(gmat_lcp1, variants, lcp1_burden_rsid, controls) -> lcp1_controls_geno

#get_geno(gmat_tgoln2, variants, tgoln2_burden_rsid, cases) -> tgoln2_cases_geno
#get_geno(gmat_tgoln2, variants, tgoln2_burden_rsid, controls) -> tgoln2_controls_geno

pheno_df <- merge(pheno_df, retsat_cases_geno, by="patient") %>% 
  merge(lcp1_cases_geno, by = "patient") 

print_group <- function(pheno_df, column){
  pheno_df %>% 
    {.[, c("patient", "medication")]} %>% 
    {.[!duplicated(.),]} %>% 
    group_by(medication) %>% 
    count() %>% 
    mutate(n_p = n/nrow(pheno_df)) -> tot_med
  
  
  pheno_df %>% 
    filter((!!sym(column)) != 0) %>%
    {.[, c("patient", "medication")]} %>% 
    {.[!duplicated(.),]} %>% 
    group_by(medication) %>% 
    count() %>%  
    merge(tot_med, by="medication", suffixes = c("", "_total")) %>% 
    {.[order(.$n,decreasing = T),]} -> sum_med
    
  sum_med$n_p_sub <- sum_med$n/sum(sum_med$n)
  sum_med$delta_p <- (sum_med$n_p_sub - sum_med$n_p)*100
  
  print("####### Medication delta #######")
  print(paste("--------", column, "--------"))
  print(sum_med[order(sum_med$delta_p, decreasing = T),])
  print("----------------------")
  print("")
  
  pheno_df %>% 
    {.[, c("patient", "diagnosis")]} %>% 
    {.[!duplicated(.),]} %>% 
    group_by(diagnosis) %>% 
    count() %>% 
    mutate(n_p = n/nrow(pheno_df)) -> tot_diag
    
  pheno_df %>% 
    filter((!!sym(column)) != 0) %>% 
    {.[, c("patient", "diagnosis")]} %>% 
    {.[!duplicated(.),]} %>% 
    group_by(diagnosis) %>%
    count() %>% 
    merge(tot_diag, by="diagnosis", suffixes = c("", "_total")) %>% 
    {.[order(.$n, decreasing = T),]} -> sum_diag
  
  sum_diag$n_p_sub <- sum_diag$n/sum(sum_diag$n)
  sum_diag$delta_p <- (sum_diag$n_p_sub - sum_diag$n_p)*100 
  
  print("####### Diagnosis delta #######")
  print(paste("--------", column, "--------"))
  print(sum_diag[order(sum_diag$delta_p, decreasing = T),])
  print("----------------------")
  print("")
  
}

print(retsat_burden_rsid)
print_group(pheno_df, retsat_burden_rsid[1])
print_group(pheno_df, retsat_burden_rsid[2])

print(lcp1_burden_rsid)
print_group(pheno_df, lcp1_burden_rsid[1])
print_group(pheno_df, lcp1_burden_rsid[2])


case_carriers  <- merge(retsat_cases_geno, lcp1_cases_geno, by = "patient")
cor(case_carriers[2:5])


retsat_carriers <- case_carriers %>% filter(rs143283662 == 1) %>% {.$patient}
lcp1_carriers <- case_carriers %>% filter(rs10492451 == 1) %>% {.$patient}


intersect(retsat_carriers, lcp1_carriers)
venn.diagram(
  x = list(retsat_carriers, lcp1_carriers),
  category.names = c("RETSAT" , "LCP1 "),
  filename = 'LCP1vsRETSAT.png',
  output=TRUE)




retsat_carriers <- pheno_df %>% filter(rs143283662 == 1) 
