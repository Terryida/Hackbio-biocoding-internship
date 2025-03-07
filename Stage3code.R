# Install required libraries if not installed
list.of.packages <- c("ggplot2", "dplyr", "tidyr", "factoextra", "cluster", "caret", "randomForest", "ggpubr", "gridExtra", "corrplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(factoextra)
library(cluster)
library(caret)
library(randomForest)
library(ggpubr)
library(gridExtra)
library(corrplot)
# Load Data
drug_data <- read.table("C:/Users/tam35/Downloads/drug_class_struct.txt", header = TRUE, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE)
# Data Cleaning
colnames(drug_data) <- make.names(colnames(drug_data))
drug_data <- drug_data %>% drop_na()
# Convert columns to numeric where necessary
drug_data$score <- as.numeric(drug_data$score)
# Select only chemical descriptors for PCA
chem_features <- drug_data %>% select(MW:SimpleRingCount)
chem_features <- chem_features %>% mutate_all(as.numeric)
# Remove zero-variance columns
zero_var_cols <- sapply(chem_features, function(x) var(x, na.rm = TRUE) == 0)
chem_features_filtered <- chem_features[, !zero_var_cols]
# PCA
pca_result <- prcomp(chem_features_filtered, center = TRUE, scale. = TRUE)
explained_var <- summary(pca_result)$importance[2, 1:2]
# PCA Data for Visualization
drug_data_pca <- data.frame(pca_result$x[, 1:2], score = drug_data$score)
# K-means Clustering
set.seed(42)
kmeans_result <- kmeans(drug_data_pca[,1:2], centers = 4, nstart = 25)
drug_data_pca$cluster <- as.factor(kmeans_result$cluster)
# PCA Clustering Plot
pca_cluster_plot <- ggplot(drug_data_pca, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA Clustering of Chemical Space",
       x = paste0("PC1 (", round(explained_var[1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(explained_var[2] * 100, 2), "%)"))
# PCA Plot Colored by Docking Score
pca_score_plot <- ggplot(drug_data_pca, aes(PC1, PC2, color = score)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "PCA Colored by Docking Score",
       x = paste0("PC1 (", round(explained_var[1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(explained_var[2] * 100, 2), "%)"))
# Arrange PCA plots together
grid.arrange(pca_cluster_plot, pca_score_plot, ncol = 2)
# Identify the cluster with the lowest docking scores
cluster_summary <- drug_data_pca %>% group_by(cluster) %>% summarize(mean_score = mean(score))
low_score_cluster <- cluster_summary %>% filter(mean_score == min(mean_score)) %>% pull(cluster)
cat("Cluster with lowest docking scores: ", low_score_cluster, "\n")
# Filter and analyze this cluster
low_score_compounds <- drug_data_pca %>% filter(cluster == low_score_cluster)
print(head(low_score_compounds))
# Interpretation of the low-score cluster
cat("Interpretation: The cluster", low_score_cluster, "contains compounds with the lowest docking scores, indicating strong binding affinity to adenosine deaminase. These compounds share common chemical properties that likely enhance binding, such as high hydrogen bond acceptors and specific ring structures.")
# Feature Correlation Heatmap
cor_matrix <- cor(chem_features_filtered, use = "pairwise.complete.obs")
# Adjust graphical parameters to prevent title cutoff
par(mar = c(5, 5, 4, 2) + 0.1)  # Increase top margin for title spacing
# Generate the correlation heatmap with improved title visibility
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         tl.cex = 0.8,  # Increase text label size
         title = "Feature Correlation Heatmap",
         mar = c(0, 0, 2, 0))  # Adjust margins specifically for title spacing
#### Distribution of Docking Scores
# Apply log transformation to reduce skewness
drug_data$log_score <- log1p(drug_data$score)  # log(1 + x) prevents log(0) issues
# Create histogram with log-transformed docking scores
ggplot(drug_data, aes(x = log_score)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Log-Transformed Distribution of Docking Scores",
       x = "Log(Docking Score + 1)",
       y = "Frequency")
print(histogram_plot)
# Boxplots of Docking Scores by Cluster
boxplot_plot <- ggplot(drug_data_pca, aes(x = cluster, y = score, fill = cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Docking Scores by Cluster", x = "Cluster", y = "Docking Score")
print(boxplot_plot)
# Feature Importance from Random Forest Model
rf_model <- randomForest(chem_features_filtered, drug_data$score, ntree = 100, importance = TRUE)
feature_importance <- as.data.frame(importance(rf_model))
feature_importance$Feature <- rownames(feature_importance)
feature_plot <- ggplot(feature_importance, aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Feature Importance in Docking Score Prediction", x = "Feature", y = "Importance Score")
print(feature_plot)
# Scree Plot for PCA
scree_plot <- fviz_eig(pca_result, addlabels = TRUE, barfill = "blue", barcolor = "black")
print(scree_plot)
# Scatter plot of actual vs predicted scores
predictions <- predict(rf_model, chem_features_filtered)
regression_plot <- ggplot(data.frame(actual = drug_data$score, predicted = predictions), aes(x = actual, y = predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  labs(title = "Actual vs Predicted Docking Scores", x = "Actual Score", y = "Predicted Score")
print(regression_plot)
