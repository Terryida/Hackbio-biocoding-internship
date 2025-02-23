# Code for 2.3 Botany and plant science
#Some scientists are trying to engineer mutants for a crop to become resistant to a pesticide.They compared the metabolic response of the engineered mutants to the metabolic response of the wild type plants
#The took readings after 8 and 24 hours
# 1. Load and Orient the Dataset
# ------------------------------
file_path <- "C:/Users/tam35/Downloads/Pesticide_treatment_data.txt"
# Load without specifying row.names
data_raw <- read.table(file_path, header = TRUE, sep = "\t")
# If needed, transpose or set row names.
if (!("WT_DMSO_1" %in% rownames(data_raw))) {
  if ("WT_DMSO_1" %in% colnames(data_raw)) {
    data_fixed <- as.data.frame(t(data_raw))
  } else if ("WT_DMSO_1" %in% data_raw[[1]]) {
    rownames(data_raw) <- data_raw[[1]]
    data_fixed <- data_raw[ , -1, drop=FALSE]
  } else {
    data_fixed <- data_raw
  }
} else {
  data_fixed <- data_raw
}
# Double-check row names
if (!("WT_DMSO_1" %in% rownames(data_fixed))) {
  stop("Could not find 'WT_DMSO_1' in row names. Adjust your data or code.")
}
# 2. Calculate ΔM (24h pesticide - DMSO)
# ------------------------------
delta_WT <- as.numeric(data_fixed["WT_pesticide_24h_1", ]) -
  as.numeric(data_fixed["WT_DMSO_1", ])
delta_mutant <- as.numeric(data_fixed["mutant_pesticide_24h_1", ]) -
  as.numeric(data_fixed["mutant_DMSO_1", ])
cat("Length of delta_WT:", length(delta_WT), "\n")
cat("Length of delta_mutant:", length(delta_mutant), "\n")
# 3. Compute Residuals and Assign Colors
residuals <- delta_mutant - delta_WT
cutoff <- 0.3
point_colors <- ifelse(abs(residuals) <= cutoff, "grey", "salmon")
# Identify outlier metabolites
outlier_metabolites <- names(residuals[abs(residuals) > cutoff])
cat("Outlier metabolites (|residual| > ±", cutoff, "):\n", sep="")
print(outlier_metabolites)
# 4. Select 6 Metabolites for Time-Course
if (length(outlier_metabolites) < 6) {
  all_mets <- colnames(data_fixed)
  needed <- 6 - length(outlier_metabolites)
  add_mets <- setdiff(all_mets, outlier_metabolites)
  selected_mets <- unique(c(outlier_metabolites, head(add_mets, needed)))
} else {
  selected_mets <- head(outlier_metabolites, 6)
}
cat("\nSelected metabolites for time-course:\n")
print(selected_mets)
# 5. Scatter Plot: WT vs. Mutant ΔM
par(mfrow = c(1, 1))  # Ensure single-plot layout in RStudio
plot(delta_WT, delta_mutant, pch = 19, col = point_colors,
     xlab = "WT ΔM (24h - DMSO)",
     ylab = "Mutant ΔM (24h - DMSO)",
     main = "Scatter: WT vs. Mutant")
abline(0, 1, lwd = 2)
# 6. Time-Course Plots (Each Shows One at a Time in RStudio)
times <- c(0, 8, 24)
for (met in selected_mets) {
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 6))  # Increase right margin for legend
  wt_vals <- as.numeric(data_fixed[c("WT_DMSO_1",
                                     "WT_pesticide_8h_1",
                                     "WT_pesticide_24h_1"), met])
  mut_vals <- as.numeric(data_fixed[c("mutant_DMSO_1",
                                      "mutant_pesticide_8h_1",
                                      "mutant_pesticide_24h_1"), met])
  plot(times, wt_vals, type = "o", col = "blue", lwd = 2,
       ylim = range(c(wt_vals, mut_vals), na.rm = TRUE),
       xlab = "Time (h)", ylab = met,
       main = paste("Time course:", met))
  lines(times, mut_vals, type = "o", col = "red", lwd = 2)
  # Legend positioned outside the plot
  legend("topright", legend = c("WT", "Mutant"),
         col = c("blue", "red"), lty = 1, lwd = 2, cex = 0.8,
         xpd = TRUE, inset = c(-.8, 0))  # Moves it further right
  Sys.sleep(1.5)  # Optional: Pause for 1.5s before showing the next plot
}
# 1. **Calculate ΔM for WT and Mutants**
#    - We calculated the metabolic response difference (ΔM) by subtracting
#      the metabolite levels in the 24-hour pesticide treatment from the DMSO control.
# 2. **Scatter Plot of ΔM for WT vs Mutants**
#    - A scatter plot was generated where each point represents a metabolite,
#      with ΔM for WT on the x-axis and ΔM for mutants on the y-axis.
# 3. **Fit a line with y-intercept = 0 and slope = 1**
#    - The fitted line `y = x` (solid diagonal line) was added to the scatter plot.
# 4. **Calculate residuals from the fitted line**
#    - The residuals were computed as `residuals = delta_mutant - delta_WT`.
# 5. **Color metabolites based on residuals**
#    - **Grey:** Metabolites where residuals are within ±0.3.
#    - **Salmon:** Metabolites where residuals are outside ±0.3.
# 6. **What are these metabolites? Trends?**
#    - Outlier metabolites with large residuals were identified.
#    - A trend was observed where some metabolites had much higher ΔM in
#      mutants compared to WT (above the y=x line), while others showed
#      the opposite (below the y=x line).
#    - This suggests that the mutant may have **deregulated** metabolism,
#      affecting specific pathways more strongly.
# 7. **Plot time-course data for 6 outlier metabolites**
#    - Time-course plots were created for 6 metabolites that had large residuals.
#    - WT and Mutant data were plotted over **0h, 8h, and 24h**.
# 8. **What do the time-course plots show?**
#    - Some metabolites show a **drastic difference** between WT and Mutant.
#    - Certain metabolites exhibit a **strong increase or decrease** over time in mutants,
#      suggesting a mutation-driven impact on metabolism.
#    - These trends could indicate which metabolic pathways are most affected
#      by the mutation, potentially revealing **new biological insights**.
## Hackbio task 2.4 Biochemistry and oncology
# Proteins structures are known to be strongly connected to their functions. However, at the amino acid level, not all amino acids contribute to structure and function equally. Galardini and colleagues decided to investigate the impact of all possible individual, non synonymous nonsense mutations on the structure and function of protein.
# The functional impact was computed as SIFT scores and the structural impact was calculated as FoldX Score (in kCal/mol).
# Import both sift and foldx datasets;
sift_dataset = 'https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/sift.tsv'
sift_dataset1 = read.table(sift_dataset,header = TRUE)
head(sift_dataset1)
foldx_dataset = 'https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/foldX.tsv'
foldx_dataset1 = read.table(foldx_dataset,header = TRUE)
head(foldx_dataset1)
# In both datasets, create a column specific_Protein_aa which will be a cantenation of the Protein and Amino_acid columns such that If you have Protein A5A607 and Amino_acid E63D, you have specific_Protein_aa A5A607_E63D
sift_dataset1$specific_Protein_aa <- paste(sift_dataset1$Protein, sift_dataset1$Amino_Acid, sep = "_")
foldx_dataset1$specific_Protein_aa <- paste(foldx_dataset1$Protein, foldx_dataset1$Amino_Acid, sep = "_")
# Using the specific_Protein_aa column, merge sift and foldx dataset into one final dataframe.
new_mergeddata <- sift_dataset1 %>%
  dplyr::left_join(foldx_dataset1, by = join_by(specific_Protein_aa))
# According to the authors;A SIFT Score below 0.05 is deleterious and A FoldX score greater than 2 kCal/mol is deleterious. Using the criteria above, Find all mutations that have a SIFT score below 0.05 and FoldX Score above 2 (i.e: Mutations that affect both structure and function)
# Filter mutations that meet the deleterious criteria
deleterious_mutations <- subset(new_mergeddata, sift_Score< 0.05 & foldX_Score > 2)
# Print the result
print(deleterious_mutations)
# Study the amino acid substitution nomenclature
# Investigate for the amino acid that has the most functional and structural impact. Hint: Using the amino acid column, find a way to select the first amino acid. Solution here.Generate a frequency table for all the amino acids. Use table() to get the frequency table
# Extract the first amino acid from the Mutation column
deleterious_mutations$AA <- substr(deleterious_mutations$Amino_Acid.x, 1, 1)
# Generate a frequency table for the first amino acids
aa_frequency <- table(deleterious_mutations$AA)
# Print the frequency table
print(aa_frequency)
# Using the amino frequency table above, generate a barplot and pie chart to represent the frequency of the amino acids.
barplot(aa_frequency,col=aa_frequency)
pie(aa_frequency,col=aa_frequency)
# Briefly describe the amino acid with the highest impact on protein structure and function
# Glycine (Gly/G) is the amino acid with the shortest side chain, having an R-group consistent only of a single hydrogen. As a result, glycine is the only amino acid that is not chiral. Its small side chain allows it to readily fit into both hydrophobic and hydrophilic environments.This allows for greater flexibility in protein folding and enables it to fit into tight spaces within a protein structure that other amino acids cannot.
# What can you say about the structural property and functional property of amino acids with more than 100 occurences.
# Amino acids with more than 100 occurrences play crucial structural and functional roles: hydrophobic residues (A, F, I, L, V, W) stabilize protein cores, while flexible residues (G, P) influence loops and turns, and polar/charged residues (D, Q, R, S, T, Y) contribute to enzymatic activity and interactions. Functionally, these amino acids are essential for protein folding, stability, binding interactions, and post-translational modifications, making them critical for both structure and function.
## Hackbio task 2.6
##Transcriptomics
##This is a processed RNAseq dataset involving reading in quantitated gene expression data from an RNA-seq experiment, exploring the data using base R functions and then interpretation. The dataset contains an experiment between a diseased cell line and diseased cell lines treated with compound X. The difference in expression change between the two health status is computed as Fold change to log 2 (Log2FC) and the significance of each is computed in p-value.
##Load the datasets
Assignment = 'https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt'
Assignment1 = read.table(Assignment,header = TRUE)
head(Assignment1)
# Load the libraries
library(ggplot2)
library(dplyr)
# # Generate a volcano plot.
ggplot(Assignment1, aes(x = log2FoldChange, y = -log10(pvalue)
                       )) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  labs(title = "Volcano Plot", x = "Expression", y = "P-Value") +
  theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
Assignment1$diffexpressed <- "NO"
# Determine the upregulated genes (Genes with Log2FC > 1 and pvalue < 0.01)
Assignment1$diffexpressed[Assignment1$log2FoldChange > 1 & Assignment1$pvalue < 0.01] <- "UP"
# Determine the downregulated genes (Genes with Log2FC < -1 and pvalue < 0.01)
Assignment1$diffexpressed[Assignment1$log2FoldChange < -1 & Assignment1$pvalue < 0.01] <- "DOWN"
# What are the functions of the top 5 upregulated genes (Use genecards)
head(Assignment1[order(Assignment1$pvalue) & Assignment1$diffexpressed == 'UP', ])
# 1. EMILIN2 (Elastin Microfibril Interfacer 2): is an extracellular matrix protein that plays a crucial role in anchoring smooth muscle cells to elastic fibers, contributing to the formation and assembly of blood vessels.
# 2. POU3F4 (POU Class 3 Homeobox 4):is a transcription factor involved in the development of the inner ear and certain neural tissues. Mutations in this gene are associated with X-linked nonsyndromic hearing loss.
# 3. LOC285954:refers to a long non-coding RNA (lncRNA) whose specific function is not well-characterized. However, some studies suggest it may play a role in tumor suppression and regulation of cell proliferation.
# 4. VEPH1 (Ventricular Zone Expressed PH Domain-Containing 1):is believed to be involved in neurodevelopment, particularly in the formation of the ventricular zone during brain development.
# 5. DTHD1 (Death Domain-Containing 1):contains a death domain, suggesting a role in apoptotic processes.
# What are the functions of the top 5 downregulated genes.
head(Assignment1[order(Assignment1$pvalue) & Assignment1$diffexpressed == 'DOWN', ])
# 1. TBX5 (T-box transcription factor 5):is a transcription factor crucial for the development of the heart and upper limbs during embryonic growth.
# 2. IFITM1 (Interferon-Induced Transmembrane Protein 1):is part of the interferon-induced transmembrane protein family and plays a role in the immune response against viral infections.
# 3. LAMA2 (Laminin Subunit Alpha 2):encodes the alpha-2 chain of laminin, a protein essential for the structural scaffolding of basement membranes in various tissues.
# 4. CAV2 (Caveolin 2): encodes caveolin-2, a protein integral to the formation of caveolae—small invaginations on the cell membrane involved in various cellular processes, including signal transduction, lipid regulation, and endocytosis.
# 5. TNN (Tenascin N):encodes tenascin-N, an extracellular matrix glycoprotein involved in tissue development and repair. It plays a role in cell adhesion, migration, and proliferation, particularly within the nervous system and during wound healing.
## Hackbio task 2.7
##Public Health
##NHANES is a program run by the CDC to assess the health and nutritional status of adults and children in the US. It combines survey questions and physical examinations, including medical and physiological measurements and laboratory tests, and examines a representative sample of about 5,000 people each year. The data is used to determine the prevalence of diseases and risk factors, establish national standards, and support epidemiology studies and health sciences research. This information helps to develop public health policy, design health programs and services, and expand the nation's health knowledge.
##Load the datasets
task = 'https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/nhanes.csv'
task4 = read.csv(task)
head(task4)
summary(factor(is.na(task4)))
#Process all NA (either by deleting or by converting to zero)
# Replace NAs in data frame with 0
task4[is.na(task4)] <- 0
# View data frame
task4
summary(factor(is.na(task4)))
# Visualize the distribution of BMI, Weight, Weight in pounds (weight *2.2) and Age with an histogram.
ggplot(task4, aes(BMI)) +
  geom_histogram(bins = 20,fill = "black") +
  labs(x = "BMI(Kgs/m2)", y = "Frequency") +
  theme_classic() +
  theme(text = element_text(size = 12, face = "bold"))
ggplot(task4, aes(Weight)) +
  geom_histogram(bins = 20,fill = "blue") +
   labs(x = "Weight(Kgs)", y = "Frequency") +
  theme_classic() +
  theme(text = element_text(size = 12, face = "bold"))
ggplot(task4, aes(Weight*2.2)) +
  geom_histogram(bins = 20,fill = "green") +
  labs(x = "Weight(pounds)", y = "Frequency") +
  theme_classic() +
  theme(text = element_text(size = 12, face = "bold"))
ggplot(task4, aes(Age)) +
  geom_histogram(bins = 20,fill="yellow") +
   labs(x = "Age(Years)", y = "Frequency") +
  theme_classic() +
  theme(text = element_text(size = 12, face = "bold"))
#What’s the mean 60-second pulse rate for all participants in the data?
mean(task4$Pulse)
# What’s the range of values for diastolic blood pressure in all participants? (Hint: see help for min(), max()).
range(task4$BPDia)
# What’s the variance and standard deviation for income among all participants?
var(task4$Income)
sd(task4$Income)
# Visualize the relationship between weight and height ?
# Color the points by
# gender
# diabetes
# smoking status
library(ggpubr)
ggscatter(task4, x = "Weight", y ="Height",
   color = "Gender", #shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.y=250, label.x =120,
                         label.sep = "\n")) +
  labs(x = "Weight(in Kgs)", y = "Height(in cm)")+
 theme(text = element_text(size = 10))
ggscatter(task4, x = "Weight", y ="Height",
   color = "Diabetes", #shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.y=250, label.x =120,
                         label.sep = "\n")) +
  labs(x = "Weight(in Kgs)", y = "Height(in cm)")+
 theme(text = element_text(size = 10))
ggscatter(task4, x = "Weight", y ="Height",
   color = "SmokingStatus", #shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.y=250, label.x =120,
                         label.sep = "\n")) +
  labs(x = "Weight(in Kgs)", y = "Height(in cm)")+
 theme(text = element_text(size = 10))
# Conduct t-test between the following variables and make conclusions on the relationship between them based on P-Value
# Age and Gender
# BMI and Diabetes
# Alcohol Year and Relationship Status
# Independent t-test
t.test(Age ~ Gender, data = task4)
# Since the p-value (p=0.08022) is greater than 0.05, we fail to reject the null hypothesis and conclude that there is no statistically significant difference in age between females and males.
t.testBMI <- task4 %>%
  filter(!Diabetes %in% "0")
t.test(BMI ~ Diabetes, data = t.testBMI)
# Since the p-value is less than 0.05 (p<2.2×10 −16), we reject the null hypothesis and conclude that there is a statistically significant difference in BMI between individuals with and without diabetes.
t.testAlcohol <- task4 %>%
  filter(!RelationshipStatus %in% "0")
t.test(AlcoholYear ~ RelationshipStatus, data = t.testAlcohol)
# Since the p-value (p=3.388×10−8) is less than 0.05, we reject the null hypothesis and conclude that there is a statistically significant difference in yearly alcohol consumption between individuals who are in a committed relationship and those who are single.
