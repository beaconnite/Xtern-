# Load necessary libraries
library(DESeq2)
library(S4Vectors)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tximport)

# Sample names and conditions setup
sample_names <- paste("Sample", 1:18, sep="")
conditions <- rep(c("Treatment", "Control"), each = 9)  # Assuming an even split

# Create colData as a DataFrame with conditions
colData <- DataFrame(sample_name = sample_names, condition = factor(conditions))
rownames(colData) <- sample_names

# Define the directory containing quant.sf files
data_dir <- "C:/Users/wcp/Downloads/B646 Example Files/salmon_base/salmon_base"
quant_files <- list.files(path = data_dir, pattern = "quant.sf$", full.names = TRUE, recursive = TRUE)

# Define and load tx2gene mapping
tx2gene_path <- file.path(data_dir, "tx2gene.csv")
if (!file.exists(tx2gene_path)) {
  stop("tx2gene file not found. Check the path.")
}
tx2gene <- read.csv(tx2gene_path, header = TRUE)

# Load quantification data using tximport
txi <- tximport(files = quant_files, type = "salmon", txOut = FALSE, tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Make sure counts are integers
counts <- round(txi$counts)  # Rounding to ensure integer counts if necessary

# Set column names for counts based on sample_names

colnames(counts) <- sample_names

# Ensure colData is ordered the same as columns in counts
colData <- colData[colnames(counts), , drop = FALSE]

# Check for dimension consistency between counts and colData
if (ncol(counts) != nrow(colData)) {
  stop("The dimensions of colData do not match the counts matrix.")
}

# Initialize DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# Run the DESeq2 analysis
dds <- DESeq(dds)
results <- results(dds)
write.csv(as.data.frame(results), file = "DESeq2_Results.csv")

# Extract results for Treatment vs Control comparison
results_comparison <- results(dds, contrast = c("condition", "Treatment", "Control"))

# Visualizing results with MA Plot for the comparison
ma_plot <- plotMA(results_comparison, main = "MA Plot: Treatment vs Control", ylim = c(-2, 2))
ggsave("MA_Plot_Treatment_vs_Control.pdf", plot = ma_plot, width = 10, height = 8)

# EnhancedVolcano Plot
volcano_plot <- EnhancedVolcano(results_comparison,
                                lab = ifelse(results_comparison$padj < 0.05 & abs(results_comparison$log2FoldChange) > 0.5,
                                             rownames(results_comparison), ""),
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                title = 'Volcano plot: Treatment vs Control',
                                pCutoff = 0.05,
                                FCcutoff = 0.5)
ggsave("Volcano_Plot_Treatment_vs_Control.pdf", plot = volcano_plot, width = 10, height = 8)

# Gene enrichment analysis
sig_downregulated_genes <- subset(results, padj < 0.05 & log2FoldChange < -0.5)
geneList <- rownames(sig_downregulated_genes)
entrez_ids <- bitr(geneList, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- na.omit(entrez_ids$ENTREZID)

# Perform enrichment analysis
if (length(entrez_ids) > 0) {
  ego <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)
  
  # Dot plot
  dp <- dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment Analysis")
  ggsave("GO_Enrichment_Dot_Plot.pdf", plot = dp, width = 10, height = 8)
} else {
  warning("No valid ENTREZ IDs available after processing.")
}



# Load necessary libraries
library(DESeq2)

# Assuming dds is already initialized and prepped for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)

# Define the contrast of interest (adjust based on your column names)

results <- results(dds, contrast = c("condition", "Treatment", "Control"))
print(resultsNames(dds))
# Check for any errors or issues in contrasts
tryCatch({
  results <- results(dds, contrast = c("condition", "Treatment", "Control"))
  print(summary(results))  # Optionally print summary to check results
}, error = function(e) {
  print(e)
  print("Please check your factor levels and contrast definition.")
})


# Filtering results for significant downregulation
# Set thresholds
padj_threshold <- 0.05
log2FoldChangeThreshold <- -1  # Threshold for downregulation

# Filter results while handling NAs in the logical vector
downregulated_genes <- results[!is.na(results$padj) & results$padj < padj_threshold &
                                 !is.na(results$log2FoldChange) & results$log2FoldChange < log2FoldChangeThreshold, ]

# Convert the results to a data frame for easier handling
downregulated_genes_data <- as.data.frame(downregulated_genes)
downregulated_genes_data$gene <- rownames(downregulated_genes_data)

# Save to CSV
write.csv(downregulated_genes_data, "Downregulated_Genes_in_HPV.csv", row.names = FALSE)
print("Saved downregulated genes data to 'Downregulated_Genes_in_HPV.csv'.")
# Check for NA prevalence in padj and log2FoldChange
cat("NA in padj:", sum(is.na(results$padj)), "\n")
cat("NA in log2FoldChange:", sum(is.na(results$log2FoldChange)), "\n")

library(DESeq2)
library(ggplot2)

# Create a customized volcano plot for downregulated genes
volcano_plot <- ggplot(downregulated_genes_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < padj_threshold), alpha = 0.8) +
  scale_color_manual(values = c("#3498DB", "#E74C3C")) +  # Blue for non-significant, Red for significant
  labs(title = "Volcano Plot of Downregulated Genes",
       x = "Log2 Fold Change",
       y = "-log10(adjusted P-value)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  geom_vline(xintercept = log2FoldChangeThreshold, linetype = "dashed", color = "#95A5A6") +  # Light gray for threshold line
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "#95A5A6")

# Save the volcano plot
ggsave("volcano_plot_downregulated_genes.png", plot = volcano_plot, width = 10, height = 8)

# Bar plot of log2 fold changes for downregulated genes with red gradient colors
bar_plot <- ggplot(downregulated_genes_data, aes(x = reorder(gene, log2FoldChange), y = log2FoldChange)) +
  geom_bar(stat = "identity", aes(fill = log2FoldChange), show.legend = FALSE) +
  scale_fill_gradient(low = "#FADBD8", high = "#CB4335") +  # Light red to deep red
  coord_flip() +
  labs(title = "Log2 Fold Changes of Significantly Downregulated Genes",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Save the bar plot
ggsave("bar_plot_downregulated_genes.png", plot = bar_plot, width = 10, height = 8)
## Machine Learning Analysis ##
# Load necessary libraries
library(DESeq2)
library(caret)
library(readr)
library(pROC)
library(ggplot2)
library(tidyr)  # for handling missing values
library(dplyr)  # for data manipulation

# Read the dataset
downregulated_genes_path <- "C:/Users/wcp/Downloads/B646 Example Files/salmon_base/salmon_base/Enhanced_Downregulated_Genes_in_HPV.csv"
downregulated_genes <- read_csv(downregulated_genes_path)

# Remove rows with NA in critical columns
downregulated_genes <- downregulated_genes %>% 
  drop_na(c("HPV_status", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))

# Define the response and features
response <- "HPV_status"
features <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

# Ensure that HPV_status is a factor with correct levels
downregulated_genes[[response]] <- factor(downregulated_genes[[response]], levels = c("Negative", "Positive"))

# Split the data into training and testing sets
set.seed(123)  # for reproducibility
index <- createDataPartition(downregulated_genes[[response]], p = 0.8, list = FALSE)
train_data <- downregulated_genes[index, ]
test_data <- downregulated_genes[-index, ]

# Define training control
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE, savePredictions = "final", summaryFunction = twoClassSummary)

# Train the logistic regression model
model_formula <- as.formula(paste(response, "~", paste(features, collapse = " + ")))
model <- train(model_formula,
               data = train_data,
               method = "glm",
               family = "binomial",
               trControl = train_control)

# Print the model
print(model)

# Make predictions on the test set
predictions <- predict(model, newdata = test_data, type = "prob")
predicted_classes <- predict(model, newdata = test_data, type = "raw")

# Ensure predicted classes are a factor with the same levels as the test labels
predicted_classes <- factor(predicted_classes, levels = levels(test_data[[response]]))

# Evaluate model performance
conf_matrix <- confusionMatrix(predicted_classes, test_data[[response]])
print(conf_matrix)

# Calculate the AUC for the ROC curve
roc_result <- roc(response = test_data[[response]], predictor = as.numeric(predictions[, "Positive"]))
auc_value <- auc(roc_result)
print(auc_value)

# Plot the ROC Curve
roc_plot <- ggroc(roc_result)
# Correct the directory path
downregulated_genes_path <- "C:/Users/wcp/Downloads/B646 Example Files/salmon_base/salmon_base/"

# Assuming the rest of your setup is correct, here's the corrected ggsave call
# Ensure the directory exists or create it
if (!dir.exists(downregulated_genes_path)) {
  dir.create(downregulated_genes_path, recursive = TRUE)
}
# Save the ROC plot correctly
ggsave(
  file = paste0(downregulated_genes_path, "ROC_Curve.png"), 
  plot = roc_plot, 
  width = 7, 
  height = 7
)
# Save model results and predictions to CSV
write.csv(as.data.frame(conf_matrix$table), file.path(downregulated_genes_path, "confusion_matrix.csv"))
write.csv(as.data.frame(predictions), file.path(downregulated_genes_path, "prediction_probabilities.csv"))



# GO enrichment for significant top and bottom DEGs
sig_degs <- rownames(top_bottom_degs)[top_bottom_degs$pvalue < 0.05]
entrez_ids <- bitr(sig_degs, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
ego <- enrichGO(gene = entrez_ids$ENTREZID[!is.na(entrez_ids$ENTREZID)],
                OrgDb = org.Mm.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Check for GO term results and create a pie chart with labels in the legend
if (!is.null(ego) && nrow(ego) > 0) {
  ego_df <- as.data.frame(ego)
  
  # Select a reasonable number of terms for clarity in the pie chart
  top_ego_df <- head(ego_df, 10)
  
  # Calculate percentages
  top_ego_df$Percentage <- round(100 * top_ego_df$Count / sum(top_ego_df$Count), 1)
  
  # Create labels that combine the GO term description and the percentage for the legend
  top_ego_df$LegendLabels <- paste(top_ego_df$Description, "(", top_ego_df$Percentage, "%)", sep = "")
  
  # Create a pie chart with labeled slices
  p <- ggplot(top_ego_df, aes(x = "", y = Count, fill = LegendLabels)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +  # Transform bar chart into pie chart
    theme_void() +  # Clean chart look
    scale_fill_viridis_d(name = "GO Terms", guide = guide_legend(title.position = "top", title.hjust = 0.5)) +  # Color scale with readable legend
    ggtitle("GO Enrichment Analysis Pie Chart") +  # Add chart title
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center the title and adjust the font
          legend.position = "right",  # Adjust legend position
          legend.text = element_text(size = 10))  # Adjust legend text size
  
  # Save the pie chart with labels in the legend
  ggsave("GO_Enrichment_Pie_Chart_with_Legend.pdf", plot = p, width = 10, height = 8)
} else {
  print("No significant GO terms found.")
}
