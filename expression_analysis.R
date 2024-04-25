# Load necessary libraries
library(readr)
library(dplyr)
library(biomaRt)

# Read the dataset
dataset_path <- "C:/Users/wcp/Downloads/B646 Example Files/salmon_base/salmon_base/Downregulated_Genes_in_HPV.csv"
downregulated_genes <- read_csv(dataset_path)

# Initialize the ENSEMBL BioMart dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Assuming 'gene' column contains the ENSG identifiers
gene_ids <- downregulated_genes$gene

# Fetch additional gene information from BioMart
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)

# Merge the new information with your original dataset
enhanced_dataset <- left_join(downregulated_genes, gene_info, by = c("gene" = "ensembl_gene_id"))

# Define a placeholder function to determine HPV status - you will need actual criteria or a database for this
determine_HPV_status <- function(gene_id) {
  # Placeholder logic - replace with actual logic to determine HPV status
  # This could be a database query or some condition based on gene properties
  if (gene_id %in% c("specific_gene_ids_known_to_be_positive")) {
    return("Positive")
  } else {
    return("Negative")
  }
}

# Add HPV status column by applying the function to each gene ID
enhanced_dataset$HPV_status <- sapply(enhanced_dataset$gene, determine_HPV_status)

# Write the new dataset to a CSV file
output_path <- "C:/Users/wcp/Downloads/B646 Example Files/salmon_base/salmon_base/Enhanced_Downregulated_Genes_in_HPV.csv"
write_csv(enhanced_dataset, output_path)
