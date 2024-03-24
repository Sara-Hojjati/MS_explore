library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# Read proteins in MS network
diagnosis_narmss <- read.table("data/ms_network.txt")
diagnosis_narmss <- unique(diagnosis_narmss$V1)
diagnosis_narmss_map <- mapIds(org.Hs.eg.db, diagnosis_narmss, 'ENTREZID', 'SYMBOL')
diagnosis_narmss_entrez <- unname(diagnosis_narmss_map[complete.cases(diagnosis_narmss_map)])

background <- as.character(read.table("data/string_proteins_background.txt")$V1)
  
### GO enrichment analysis
kk_go <- clusterProfiler::enrichGO(gene = diagnosis_narmss_entrez, ont="BP",
                                   OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                                   pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = background)
kk_go@result$significant <- (kk_go@result$p.adjust < 0.05)
write.table(kk_go@result[kk_go@result$significant,], "results/ms_network_GO_enrichment.tsv", sep='\t', quote = FALSE)
