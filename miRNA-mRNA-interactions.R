###miRNA-mRNA-interactions

# Load the miRTarBase miRNA-mRNA interaction data (adjust the file path accordingly)
miRNA_mRNA_interactions <- read.csv("hsa_MTI.csv")

# Inspect the data
head(miRNA_mRNA_interactions)

#####Filter miRNA-mRNA Interactions Using Differentially Expressed Genes

# Extract significant miRNAs and mRNAs from the DESeq2 results
sig_miRNA_symbols <- res_miRNA_annotated %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.5) %>%
  pull(hgnc_symbol)

sig_mRNA_symbols <- res_mRNA_annotated %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.5) %>%
  pull(hgnc_symbol)

# Filter miRTarBase for targets of these significant miRNAs
miRNA_target_genes <- miRNA_mRNA_interactions %>%
  filter(MIRNA %in% sig_miRNA_symbols)

# Save the miRNA-target interactions for reporting or further use
write.csv(miRNA_target_genes, "miRNA_target_genes_from_miRTarBase.csv")

# Filter miRTarBase for miRNAs that target these significant mRNAs
mRNA_targeted_by_miRNAs <- miRNA_mRNA_interactions %>%
  filter(`Target.Gene` %in% sig_mRNA_symbols)

# Save the mRNA-targeted interactions for reporting or further use
write.csv(mRNA_targeted_by_miRNAs, "miRNAs_targeting_sig_mRNAs_from_miRTarBase.csv")
