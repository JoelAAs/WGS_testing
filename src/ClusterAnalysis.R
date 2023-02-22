req_packages <- c("AnnotationDbi", "DESeq2", "devtools", "doParallel",
                  "dynamicTreeCut", "edgeR", "flashClust", "foreach",
                  "ggdendro", "ggrepel", "ggplot2", "igraph",
                  "limma", "MODA", "openxlsx", "org.Hs.eg.db",
                  "plyr", "preprocessCore", "Rcpp", "reticulate",
                  "RSQLite", "stackoverflow", "STRINGdb", "WGCNA")


#### INPUT LOAD
# Load required packages
for (i in 1:length(req_packages)) {
  library(req_packages[i], character.only = TRUE)
}

library(MODifieR)

cns_genes <- read.delim("../../Resultat/CNS-tox/indicated_genes",
                       stringsAsFactors = FALSE, header = F)

cns_genes <- cns_genes$V1

#BiocManager::install("biomaRt")
library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

entrezgene <- unique( getBM(attributes = "entrezgene",
                                values = "*",
                                mart = ensembl) )
ids <- getBM(attributes = c(
  "entrezgene_id",
  "ensembl_gene_id",
  "external_gene_name",
  "refseq_peptide"
),
filter = "entrezgene_accession",
values = cns_genes,
mart = mart)
ids$P <- 0
gene_ids <- unique(ids$entrezgene_id)
genes_ids <- data.frame(gene = gene_ids)
genes_ids$pvalue <- 0


ppi_network <- read.delim("data/PPI_network.txt",
                          stringsAsFactors = FALSE)
ppi_network[] <- lapply(ppi_network, as.character) # Format the data as "character"
class(ppi_network)
head(ppi_network)


library(DOSE)
disease_enrichment <- enrichDO(gene = gene_ids,
                               ont = "DO",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH")
print(disease_enrichment@result)

gene_input <- create_custom_microarray_input_object(diff_genes = genes_ids)

module_mcode <- mcode(MODifieR_input = gene_input,
                      ppi_network = ppi_network,
                      hierarchy = 1,
                      vwp = 0.1,
                      haircut = TRUE,
                      fluff = FALSE,
                      fdt = 0.2,
                      loops = TRUE)


module_mcode_genes <- module_mcode$module_genes

length(module_mcode_genes) # Number of genes in the module


#### Disease enrichment analysis of MCODE genes
# Disease enrichment analysis of the module genes
disease_enrichment2 <- enrichDO(gene = module_mcode_genes,
                                ont = "DO",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")

sort(disease_enrichment2@result$Description)
res2 <- disease_enrichment2@result

library(clusterProfiler)
pathEnrich <- enrichKEGG(gene = module_mcode_genes,
                         organism = "hsa", # Homo sapiens
                         keyType = "kegg", # KEGG database
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
View(pathEnrich@result)

# Visualize the top 10 significantly enriched

