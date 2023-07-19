#----Transcriptomics----

# Load the DESeq2 library for differential expression analysis of transcriptomics
library(DESeq2)

countData <- read.csv("./datasets/gene_count.csv")
getDERResults <- function(count_data, i) {
  
  # Extract the sample names from the count data
  sample <- names(count_data[c(2:7)])
  
  # Create a vector representing the experimental conditions
  dex <- c("control", "control", "control", "treatment", "treatment", "treatment")
  
  # Create a metadata data frame 
  metadata <- data.frame(sample, dex)
  
  # Create a DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~dex, tidy = TRUE)
  
  # Run the DESeq 
  dds <- DESeq(dds)
  
  # Perform differential expression analysis and obtain results
  res <- results(dds, lfcThreshold = 1, alpha = 0.05)
  
  # Filter the results based on adjusted p-value and fold change thresholds
  filter_de <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
  
  return(filter_de@rownames)
}

# T10
DEG_list_T10 <- list()
run_names <- c("C001", "C01", "C1", "C10", "C100", "C1000")

for (i in 1:6) {

  # Subset the count data for the current run
  count_data <- countData[c(1:4, (i * 3 + 2):(i * 3 + 4))]

  # Perform differential expression analysis and obtain filtered gene names
  filtered_row_names <- getDERResults(count_data, i)

  # Store the filtered gene names in the DEG_list_T10
  DEG_list_T10[[run_names[i]]] <- filtered_row_names

  # Print the number of differentially expressed genes for the current run
  print(length(filtered_row_names))
}

# T24
DEG_list_T24 <- list()
for (i in 8:13) {

  # Subset the count data for the current run
  count_data <- countData[c(1, 23:25, (i * 3 + 2):(i * 3 + 4))]
  
  # Perform differential expression analysis and obtain filtered gene names
  filtered_row_names <- getDERResults(count_data, i)

  # Store the filtered gene names in the DEG_list_T24
  DEG_list_T24[[run_names[i-7]]] <- filtered_row_names

  # Print the number of differentially expressed genes for the current run
  print(length(filtered_row_names))
}

#------Proteomics------
# Load the limma library for differential expression analysis for proteomics and phosphoproteomics
library(limma)

proteomics <- read.csv("./datasets/proteins.csv")
getLimmaResults <- function(quant_data, group_list) {

  # Filter out rows with too many missing values
  filtered_df <- quant_data[rowSums(quant_data == 0) < 3, ]
  
  # Apply log2 transformation to the expression values
  filtered_df[, 2:7] <- log2(filtered_df[, 2:7] + 1)
  
  # Create the design matrix for the linear model
  design <- model.matrix(~group_list)

  # Set column and row names of the design matrix
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(filtered_df)[-1]

  # Fit the linear model
  fit <- lmFit(filtered_df, design)

  # Perform empirical Bayes moderation of the t-statistics
  fit <- eBayes(fit, trend = TRUE)

  # Extract the top differentially expressed proteins based on adjusted p-value
  result_limma <- topTable(fit, coef = 2, n = Inf)

  # Return the filtered differentially expressed proteins
  return(result_limma[result_limma$adj.P.Val < 0.05 & abs(result_limma$logFC) > 1, ])
}

# 10 min time point
# Create an empty list to store the differentially expressed proteins at T10
DEP_list_T10 <- list()

for (i in 1:6) {
  # Subset the proteomics data for the current run at T10
  quant_data <- proteomics[c(1, 3:5, (i * 3 + 3):(i * 3 + 5))]

  # Create a factor representing the experimental conditions
  group_list <- factor(c(rep("control", 3), rep("treat", 3)))

  # Perform differential expression analysis and obtain filtered protein names
  filtered_limma_results <- getLimmaResults(quant_data, group_list)

  # Store the filtered protein names in DEP_list_T10
  DEP_list_T10[[run_names[i]]] <- filtered_limma_results

  # Print the number of differentially expressed proteins for the current run
  print(sum(filtered_limma_results$adj.P.Val < 0.05 & abs(filtered_limma_results$logFC) > 1))
}

# 24 hr time point
# Create an empty list to store the differentially expressed proteins at T24
DEP_list_T24 <- list()

for (i in 8:13) {
  # Subset the proteomics data for the current run at T24
  quant_data <- proteomics[c(1, 24:26, (i * 3 + 3):(i * 3 + 5))]

  # Create a factor representing the experimental conditions
  group_list <- factor(c(rep("control", 3), rep("treat", 3)))

  # Perform differential expression analysis and obtain filtered protein names
  filtered_limma_results <- getLimmaResults(quant_data, group_list)

  # Store the filtered protein names in DEP_list_T24
  DEP_list_T24[[run_names[i-7]]] <- filtered_limma_results
  
  # Print the number of differentially expressed proteins for the current run
  print(sum(filtered_limma_results$adj.P.Val < 0.05 & abs(filtered_limma_results$logFC) > 1))
}

#------Phosphoproteomics------
phos <- read.csv("./datasets/PhosphoSites.csv")
# Create an empty list to store the differentially expressed phosphosites at T10
DEPhos_list_T10 <- list()

for (i in 1:6) {
  # Subset the phosphoproteomics data for the current run at T10
  quant_data <- phos[c(1:4, (i * 3 + 2):(i * 3 + 4))]

  # Create a factor representing the experimental conditions
  group_list <- factor(c(rep("control", 3), rep("treat", 3)))

  # Perform differential expression analysis and obtain filtered phosphosites
  filtered_limma_results <- getLimmaResults(quant_data, group_list)

  # Store the filtered phosphoprotein names inD EPhos_list_T10
  DEPhos_list_T10[[run_names[i]]] <- filtered_limma_results
  
  # Print the number of differentially expressed phosphosites for the current run
  print(sum(filtered_limma_results$adj.P.Val < 0.05 & abs(filtered_limma_results$logFC) > 1))
}

# Runs 8 to 13
# Create an empty list to store the differentially expressed phosphosites at T24
DEPhos_list_T24 <- list()
for (i in 8:13) {
  # Subset the phosphoproteomics data for the current run at T24
  quant_data <- phos[c(1, 23:25, (i * 3 + 2):(i * 3 + 4))]

  # Create a factor representing the experimental conditions
  group_list <- factor(c(rep("control", 3), rep("treat", 3)))

  # Perform differential expression analysis and obtain filtered phosphosites
  filtered_limma_results <- getLimmaResults(quant_data, group_list)

  # Store the filtered phosphosites names in DEPhos_list_T24
  DEPhos_list_T24[[run_names[i-7]]] <- filtered_limma_results

  # Print the number of differentially expressed phosphoproteins for the current run
  print(sum(filtered_limma_results$adj.P.Val < 0.05 & abs(filtered_limma_results$logFC) > 1))
}
