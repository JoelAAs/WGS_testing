library(ggplot2)
library(plotly)
library(tidyverse)

cases <-  c("RB-1745-2486", "RB-1745-2911", "RI-1933-1142", "RI-1933-2449", "RI-1933-2462",
            "RI-1933-2463", "RI-1933-2470", "RI-1933-2487", "RI-1933-2611", "RI-1933-2630",
            "RI-1933-2687", "RI-1933-2701", "RI-1933-2712", "RI-1933-2852", "RI-1933-2868",
            "RI-1933-2872", "RI-1933-2908", "RI-1933-2914", "RI-1933-2921", "RI-1933-2935",
            "RI-1933-4379", "RI-1933-4455", "RI-1933-4489", "RI-1933-4545", "RI-1933-4546",
            "RI-1933-4547", "RI-1933-4711", "RI-1933-6136", "RI-1933-6636", "RI-1933-6766")


df <- read.csv("data/PCA/Trombocytopeni~rest.eigenvec", sep = "\t", stringsAsFactors = F)
df$case <- rep(0, nrow(df))
df[df$X.IID %in% cases, "case"] <- rep(1, length(cases))
g <- ggplot(df, aes(x=PC1, y=PC2, color=as.factor(case), label=as.factor(case) )) + geom_point()

ggplotly(g)
RI-1933-6121
RI-1933-2187

pPC2 <- 4.8*10^-3
pPC1 <- 3.5*10^-4

df_test <- read.csv("PCA_all.eigenvec", sep = "\t")
samples <- read.csv("samplefile.txt", sep = "\t")

df_test <- merge(df_test, samples, by.x = "X.IID", by.y = "samplename")

df_test <- df_test %>%  filter(PC1 > pPC1 | PC2 > pPC2)


assocs <- c(
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

phenofile <- read.csv("phenofile", sep = "\t", as.is=T, check.names = F)

for (assoc in assocs) {
  idx_unchange  <- sapply(df_test$adrtype, function(x) grepl(x, assoc))
  samp <- as.character(df_test[!idx_unchange, "X.IID"])
  phenofile[phenofile$IID %in% samp, assoc] = -9
}
