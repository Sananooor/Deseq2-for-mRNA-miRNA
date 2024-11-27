# Install BiocManager to manage Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "org.Hs.eg.db"))

# Install CRAN packages for machine learning and data manipulation
install.packages(c("dplyr", "ggplot2", "pheatmap"))

# Load required libraries
library(DESeq2)
library(biomaRt)
library(dplyr)
library(ggplot2)

# Load count data and metadata (adjust file paths accordingly)
counts <- read.csv("counts.csv", sep = ',', row.names = 1)
metadata <- read.csv("hcc-metadata.csv", row.names = 1)

# Ensure metadata contains a "Condition" column with "Normal" and "Tumor" labels
metadata$Condition <- factor(metadata$Condition, levels = c("Normal", "Tumor"))

# Check if column names in count data match the row names in metadata
all(colnames(counts) == rownames(metadata))  # Should return TRUE

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Condition)

# Pre-filtering: remove genes with very low counts
dds <- dds[rowSums(counts(dds)) > 10, ]


######Step 2: Annotate Genes as mRNA or miRNA Using biomaRt
# Load biomaRt and connect to Ensembl database for human genes
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve gene annotations for Ensembl gene IDs (e.g., mRNAs and miRNAs)
gene_annotations <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol', 'gene_biotype'),
                          filters = 'entrezgene_id',
                          values = rownames(dds),
                          mart = ensembl)

# Merge annotations with DESeq2 dataset
gene_annotations$entrezgene_id <- gene_annotations$entrezgene_id
dds_annotated <- dds[rownames(dds) %in% gene_annotations$entrezgene_id, ]

# Separate mRNAs and miRNAs based on the gene biotype
mrna_ids <- gene_annotations %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(entrezgene_id)

mirna_ids <- gene_annotations %>%
  filter(gene_biotype == "miRNA") %>%
  pull(entrezgene_id)

# Separate mRNA and miRNA counts
dds_mRNA <- dds[rownames(dds) %in% mrna_ids, ]
dds_miRNA <- dds[rownames(dds) %in% mirna_ids, ]


#####Step 3: Differential Expression Analysis for mRNA and miRNA
# Normalize and perform differential expression analysis for mRNA
dds_mRNA <- DESeq(dds_mRNA)
res_mRNA <- results(dds_mRNA, contrast = c("Condition", "Tumor", "Normal"))

# Filter significant mRNAs (adjusted p-value < 0.05 and log2FoldChange > 1)
sig_mRNA <- subset(res_mRNA, padj < 0.05 & abs(log2FoldChange) > 1)

# Normalize and perform differential expression analysis for miRNA
dds_miRNA <- DESeq(dds_miRNA)
res_miRNA <- results(dds_miRNA, contrast = c("Condition", "Tumor", "Normal"))

# Filter significant miRNAs (adjusted p-value < 0.05 and log2FoldChange > 1)
sig_miRNA <- subset(res_miRNA, padj < 0.05 & abs(log2FoldChange) > 1)


####Step 4: Annotate and Export the Results
# Annotate significant mRNAs using biomaRt
sig_mRNA_ensembl <- rownames(sig_mRNA)
mRNA_annotations <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
                          filters = 'entrezgene_id',
                          values = sig_mRNA_ensembl,
                          mart = ensembl)

# Annotate significant miRNAs using biomaRt
sig_miRNA_ensembl <- rownames(sig_miRNA)
miRNA_annotations <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
                           filters = 'entrezgene_id',
                           values = sig_miRNA_ensembl,
                           mart = ensembl)

# Merge annotations with DESeq2 results
res_mRNA$entrezgene_id <- rownames(res_mRNA)
res_mRNA_annotated <- merge(as.data.frame(res_mRNA), mRNA_annotations, by = "entrezgene_id")

res_miRNA$entrezgene_id <- rownames(res_miRNA)
res_miRNA_annotated <- merge(as.data.frame(res_miRNA), miRNA_annotations, by = "entrezgene_id")

# Save the annotated results
write.csv(res_mRNA_annotated, "DE_mRNA_annotated.csv")
write.csv(res_miRNA_annotated, "DE_miRNA_annotated.csv")


####Step 5: Visualization (Volcano Plot and Heatmap)
# Volcano Plot for mRNA
library(ggplot2)

# Create a volcano plot for mRNAs
res_mRNA_annotated$threshold <- as.factor(ifelse(res_mRNA_annotated$padj < 0.05 & res_mRNA_annotated$log2FoldChange > 0.5, "Upregulated",
                                                 ifelse(res_mRNA_annotated$padj < 0.05 & res_mRNA_annotated$log2FoldChange < -0.5, "Downregulated", "Not Significant")))

ggplot(res_mRNA_annotated, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = threshold), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "lightblue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot - mRNA", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme(legend.title = element_blank())

# Create a volcano plot for miRNAs
res_miRNA_annotated$threshold <- as.factor(ifelse(res_miRNA_annotated$padj < 0.05 & res_miRNA_annotated$log2FoldChange > 0.5, "Upregulated",
                                                  ifelse(res_miRNA_annotated$padj < 0.05 & res_miRNA_annotated$log2FoldChange < -0.5, "Downregulated", "Not Significant")))

ggplot(res_miRNA_annotated, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = threshold), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "darkred", "Downregulated" = "green", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot - miRNA", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme(legend.title = element_blank())


# Heatmap of top 10 differentially expressed mRNAs
library(pheatmap)

# Select the top 10 most differentially expressed mRNAs based on log2FoldChange
# Step 1: Ensure that the top 10 mRNA genes are correctly mapped using Ensembl IDs
# We will use the Ensembl IDs from res_mRNA_annotated

# Select the top 10 based on absolute log2 fold change
top10_mRNA <- res_mRNA_annotated %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10)

# Extract Ensembl IDs of the top 10 genes
top10_entrez_mRNA <- top10_mRNA$entrezgene_id

# Extract normalized counts for the top 10 mRNAs using Ensembl IDs
norm_mRNA <- assay(dds_mRNA)[rownames(dds_mRNA) %in% top10_entrez_mRNA, ]

# Match the rownames of norm_mRNA with top 10 HGNC symbols
rownames(norm_mRNA) <- top10_mRNA$hgnc_symbol

# Check if rownames are correctly assigned
head(rownames(norm_mRNA))

# Plot the heatmap for the top 10 mRNAs
pheatmap(norm_mRNA, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         main = "Top 10 Differentially Expressed mRNAs",
         annotation_col = metadata)


# Repeat for miRNAs
top10_miRNA <- res_miRNA_annotated %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10)

top10_entrez_miRNA <- top10_miRNA$entrezgene_id

norm_miRNA <- assay(dds_miRNA)[rownames(dds_miRNA) %in% top10_entrez_miRNA, ]

rownames(norm_miRNA) <- top10_miRNA$hgnc_symbol

# Plot the heatmap for the top 10 miRNAs
pheatmap(norm_miRNA, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         main = "Top 10 Differentially Expressed miRNAs",
         annotation_col = metadata)

