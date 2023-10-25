library(limma)
library(DESeq2)
library(readr)
library(dplyr)
library(readxl)
library(dplyr)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(gplots)
library(RColorBrewer)



setwd('/Users/alejandroadriaquelozano/Documents/Systems Biology/Programming/project/Data/')

#imporrt count data
lung_data <- read_delim('E-ENAD-46-raw-counts.tsv',delim = '\t')
lung_data_good <- lung_data[,3:40]
row.names(lung_data_good) <- lung_data$`Gene ID`


# colData
coldata <- read_tsv('E-ENAD-46-experiment-design.tsv')

# Create a data frame representing the design matrix
design <- data.frame(row.names = coldata$Run, Condition = coldata$`[disease]`,Tissue= coldata$`[organism part]`, Patient = coldata$`[individual]`)


#DESq2


# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = lung_data_good, colData = design, design = ~ Condition + Tissue + Patient)

# Perform differential expression analysis
dds <- DESeq(dds)

# pca
vsd <- vst(dds)

# Perform PCA on the transformed data
pca_result <- prcomp(t(assay(vsd)))

# Create a PCA plot
pca_data <- as.data.frame(pca_result$x)
rownames(pca_data) <- colnames(lung_data_good)
pca_data$Condition <- dds$Condition
pca_data$Tissue <- dds$Tissue

ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Tissue)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("circle", "triangle")) +
  theme_minimal() +
  geom_text(aes(label = rownames(pca_data)), size = 3, nudge_x = 0.1, nudge_y = -5)

# heatmap pero vamos no lo uso...
# Create a dataframe for annotation
annotation_df <- as.data.frame(colData(dds))
rownames(annotation_df) <- colnames(lung_data_good)

# Define color palettes for conditions and tissues
condition_colors <- brewer.pal(n = length(unique(annotation_df$Condition)), name = "Set1")
tissue_colors <- brewer.pal(n = length(unique(annotation_df$Tissue)), name = "Set2")

heatmap.2(as.matrix(lol), 
          scale = "row",  # Scale rows (genes)
          Rowv = FALSE,    # Do not perform hierarchical clustering on rows
          Colv = TRUE,    # Do not perform hierarchical clustering on columns
          col = colorRampPalette(c("blue", "white", "red"))(100),  # Define color palette
          dendrogram = "none",  # Do not display dendrograms
          margins = c(5, 10),   # Add extra space for row and column labels
          main = "Gene Expression Heatmap",
          xlab = "Samples",
          ylab = "Genes",
          ColSideColors = c(condition_colors[as.numeric(annotation_df$Condition)],
                  tissue_colors[as.numeric(annotation_df$Tissue)]))

# Volcano Plot
# Define custom legend labels

# Define custom legend labels
legend_labels <- c("Significant" = "DE Genes (|log2FC| > 1)",
                   "Not Significant" = "Non-DE Genes")

# Set log2 fold change cutoff
log2fc_cutoff <- 1

# Set significance cutoff (e.g., adjusted p-value < 0.05)
significance_cutoff <- 0.05

# Prepare the data, including applying the significance cutoff
lung_vs_normal$significant <- ifelse(lung_vs_normal$padj < significance_cutoff & abs(lung_vs_normal$log2FoldChange) > log2fc_cutoff, "Significant", "Not Significant")

# Customize your volcano plot
ggplot(lung_vs_normal, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), size = 2) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue"), 
                     labels = legend_labels, 
                     name = "Gene Expression Significance") +  # Specify custom legend title
  geom_hline(yintercept = -log10(significance_cutoff), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "gray") +  # Add fold change cutoff
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-log10(p-value)") +
  theme_minimal()+ 
  theme(legend.text = element_text(size = 14))


#contrast

contrast <- c("Condition", "COVID-19", "normal")

# Subset the DESeqDataSet to the specific tissue of interest
dds_lung <- dds[dds$Tissue == 'lung', ]
dds_colon <- dds[dds$Tissue == 'colon', ]

# Perform differential expression analysis for the specific tissue
results_lung <- results(dds_lung, contrast = contrast)
results_colon <- results(dds_colon, contrast = contrast)

lung_vs_normal <- data.frame(results_lung)
colon_vs_normal <- data.frame(results_colon)

# Contrast for Lung vs. Colon within Normal condition
contrast <- c("Tissue", "lung", "colon")
# Subset the DESeqDataSet 
dds_COVID <- dds[dds$Condition == 'COVID-19', ]
dds_normal <- dds[dds$Condition == 'normal', ]
# Perform differential expression analysis for the specific tissue
results_lungvscolon_covid <- data.frame(results(dds_COVID, contrast = contrast))
results_lungvscolon_normal <- data.frame(results(dds_normal, contrast = contrast))

significance_threshold <- 0.005

# Filter significant genes from the results objects

##### GO

#Filtering by p-value and getting up/down regulating genes
lung_filtered = lung_vs_normal %>% filter(padj< 0.05)
up_regulated_genes_lung = lung_filtered %>% filter(log2FoldChange>0)
down_regulated_genes_lung = lung_filtered %>% filter(log2FoldChange<0)

colon_filtered = colon_vs_normal %>% filter(padj< 0.05)
up_regulated_genes_colon = colon_filtered %>% filter(log2FoldChange>0)
down_regulated_genes_colon = colon_filtered %>%  filter(log2FoldChange<0)

lungcolon_filtered_covid = results_lungvscolon_covid %>% filter(padj< 0.05)

# I try to remove Differential expressed genes from lung vs colon (COVID) that are also present in lung vs colon (normal)
lungcolon_filtered_normal =results_lungvscolon_normal %>% filter(padj< 0.05)
lungcolon_filtered_covid$name <- rownames(lungcolon_filtered_covid)
lungcolon_filtered_normal$name <- rownames(lungcolon_filtered_normal)

lung_vs_colon <- anti_join(lungcolon_filtered_covid,lungcolon_filtered_normal, by='name')
#However no genes are present in both 'lungcolon_filtered_covid' and 'lungcolon_filtered_normal'
# We can continue the analysis using just 'lungcolon_filtered_covid'


up_regulated_genes_both = lungcolon_filtered_covid %>% filter(log2FoldChange>0)
down_regulated_genes_both = lungcolon_filtered_covid %>%  filter(log2FoldChange<0)



mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Get gene symbols
entrez_names <- getBM(c("entrezgene_id"), "ensembl_gene_id", test_input_ENSEMBL, mart)

test <- enrichGO(entrez_names$entrezgene_id, keyType = 'ENTREZID' ,ont = 'BP', OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1)
summary <- data.frame(summary(test))
dotplot(test, showCategory =10)