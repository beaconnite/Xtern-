# Load necessary libraries
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(readr)
library(randomForest)
library(e1071)  # for SVM

# Set the correct directory path for outputs
output_directory <- "C:/Users/wcp/Downloads/B646 Example Files/salmon_base/salmon_base/"
# Check if directory exists, if not create it
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Read the dataset
downregulated_genes_path <- file.path(output_directory, "Enhanced_Downregulated_Genes_in_HPV.csv")
downregulated_genes <- read_csv(downregulated_genes_path)

# Clean the data
downregulated_genes <- downregulated_genes %>%
  drop_na(c("HPV_status", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))

# Define the response and features
response <- "HPV_status"
features <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
downregulated_genes[[response]] <- factor(downregulated_genes[[response]], levels = c("Negative", "Positive"))

# Split the data
set.seed(123)
train_index <- createDataPartition(downregulated_genes[[response]], p = 0.8, list = FALSE)
train_data <- downregulated_genes[train_index, ]
test_data <- downregulated_genes[-train_index, ]

# Define training control
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE, savePredictions = "final", summaryFunction = twoClassSummary)

# Initialize an empty list to store models
models <- list()

# Train models
model_types <- c("glm", "rf", "svmRadial")
names(model_types) <- c("Logistic Regression", "Random Forest", "SVM")

for (model_type in model_types) {
  models[[model_type]] <- train(reformulate(features, response),
                                data = train_data,
                                method = model_type,
                                trControl = train_control,
                                preProcess = "scale",
                                tuneLength = 5)
}

# ... (rest of your code)

# Evaluate and plot models
roc_data <- data.frame()
for (model_type in names(models)) {
  # Predict probabilities for the positive class
  predictions <- predict(models[[model_type]], newdata = test_data, type = "prob")[, "Positive"]
  
  # Calculate ROC using the correct probability scores
  roc_result <- roc(test_data[[response]], predictions)
  
  # Store ROC data for plotting
  # The sensitivities are TPRs and (1 - specificities) are FPRs
  roc_data <- rbind(roc_data, data.frame(tpr = roc_result$sensitivities, fpr = 1 - roc_result$specificities, Model = model_type))
}

# Plot ROC curves with corrected FPR and TPR
roc_plot <- ggplot(roc_data, aes(x = fpr, y = tpr, color = Model)) +
  geom_line() +
  labs(title = "ROC Curve Comparison", x = "False Positive Rate", y = "True Positive Rate") +
  scale_color_manual(values = c("red", "blue", "green")) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))  # Ensure the axes are correct

# Save the ROC plot
ggsave(file.path(output_directory, "ROC_Curve_Comparison_Corrected.png"), plot = roc_plot, width = 10, height = 8)
