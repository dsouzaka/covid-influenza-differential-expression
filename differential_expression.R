# ============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS: COVID-19 vs Influenza
# Using limma/voom pipeline
# Author: Katelyn
# Date: December 2024
# Course: BINF 6310 - Introduction to Bioinformatics
# ============================================================================

# Load required libraries
library(limma) # Linear Models for Microarray and RNA-seq Data for statistical testing
library(edgeR) # Used for normalization and handling count data
library(ggplot2) # Makes plots
library(pheatmap) # Makes heatmaps

# ============================================================================
# STEP 1: Load Data
# ============================================================================

# Load raw counts and metadata
# Reads the TSV files into R
# header = True means that first row is column names
# columns separated by tabs, and row.names = 1 means that first column is gene IDs
counts <- read.table("selected_raw_counts.tsv", header=TRUE, sep="\t", row.names=1) # creates a matrix with genes in rows, samples in columns, and numbers are read counts
metadata <- read.table("metadata_6samples.tsv", header=TRUE, sep="\t") # creates a table telling you which sample is COVID vs Influenza


# ============================================================================
# STEP 2: Fix sample name matching
# ============================================================================

# Fix metadata sample names to match counts column names
# R automatically changes "-" to "." and adds X to numbers; must make sure names match exactly so analysis knows which sample is which
# match() function reorders metadata rows to match the column order in counts
metadata$sample <- gsub("-", ".", metadata$sample)  # Replace dashes with periods
metadata$sample[metadata$sample == "94189"] <- "X94189"  # Add X to number

# Reorder metadata to match counts column order
metadata <- metadata[match(colnames(counts), metadata$sample), ]

# Verify they match
print("Sample name check:")
print(data.frame(counts_cols = colnames(counts), 
                 metadata_samples = metadata$sample,
                 match = colnames(counts) == metadata$sample))

# ============================================================================
# STEP 3: Set up design matrix
# ============================================================================

# Create group factor
group <- factor(metadata$cohort, levels = c("COVID", "Influenza"))

# Create design matrix (no intercept model for easier contrasts)
# This tells limma which samples belong to which group
design <- model.matrix(~0 + group) # 0 + group means no intercept model and makes it easier to set up contrasts later
colnames(design) <- levels(group)



print("Design matrix:")
print(design)

# ============================================================================
# STEP 4: Filter low-count genes
# ============================================================================

# Keep genes with at least 10 counts in at least 3 samples
# Removes genes with very low counts to improve statistical power
keep <- rowSums(counts >= 10) >= 3 # must have at least 10 reads and gene must have more than 10 reads in at least 3 samples
counts_filtered <- counts[keep, ]

cat("\nGenes before filtering:", nrow(counts))
cat("\nGenes after filtering:", nrow(counts_filtered))

# ============================================================================
# STEP 5: Create DGEList and normalize
# ============================================================================

# Create DGEList object that edgeR uses; combines counts and sample info
dge <- DGEList(counts = counts_filtered, group = group)

# Calculate normalization factors (Trimmed Mean of M-values (TMM)) 
# TMM normalization adjusts for library size differences; bigger library doesn't mean more expressed
dge <- calcNormFactors(dge)

cat("\nNormalization factors:")
print(dge$samples)

# ============================================================================
# STEP 6: Voom transformation
# ============================================================================

# Apply voom transformation
# This converts counts to log2-CPM with precision weights; makes data more normally distributed which is better for statistics
# Models mean-variance relationship; low-count = more variable, high-count = less variable
# Assigns precision weights to each gene; less variable genes are stable with a high weight, low count are less table with a lower weight
pdf("voom_mean_variance_plot.pdf", width=8, height=6)
v <- voom(dge, design, plot=TRUE) # makes RNA-seq data work with limma
dev.off()
# Produces a plot that shows mean expression vs standard deviation; red line is voom's fitted model
# This is necessary because limma needs normally distributed microarray data and RNAseq count data is nt normal
# Voom transforms counts so limma can work

cat("\nVoom transformation complete. Plot saved to 'voom_mean_variance_plot.pdf'\n")

# ============================================================================
# STEP 7: Fit linear model with limma
# ============================================================================

# Fit the linear model
fit <- lmFit(v, design)
# How much of gene expression variance is explained by disease group
# For each gene, calculate the mean in COVID and Flu, use the design matrix to differentiate the samples, and accounts for precision weights voom calculated
# Produces an object containing estimated expression level for each gene in each group, variance estimates, and residuals (unexplained variation)

# Set up contrast (COVID vs Influenza)
# Defines specific comparison between COVID and Influenza
# Average expression in COVID - Average expression in Flue = log fold change
contrast.matrix <- makeContrasts(
  COVID_vs_Influenza = COVID - Influenza,
  levels = design
)

print("Contrast matrix:")
print(contrast.matrix)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrast.matrix) # applies contrast to fitted model

# Empirical Bayes moderation
# Makes statistics more reliable even with smaller sample sizes
# Looks at variance across all genes and uses information to improve variance estimates for each individual gene
# Shrinks extreme variance estimates toward the average
fit2 <- eBayes(fit2)

# ============================================================================
# STEP 8: Extract results
# ============================================================================

# Get all results
# coef = "COVID_vs_Influenza" which contrast to report
# number = Inf = return all genes
# output columns: logFC = log2fold change in COVID vs Flu, AveExpr = average expression across all samples
# cont. t = t-statistic, P.Value = raw p-value, adj.P.val = adjusted p-value (corrected for 15,556 genes)
# cont. b = log-odds that gene is differentially expressed
# sorted by most significan genes first (lowst adj.P.Val)
results <- topTable(fit2, coef="COVID_vs_Influenza", number=Inf, adjust.method="BH")

# adjust.method = "BH"
# if p < 0.05 is used for each gene, we'd get 778 significant genes expected by chance, false positives
# Benjamini-Hochberg correction controls the false discovery rate and adjusts p-value so that among significan genes onlt 5% are false positives


# Add gene names as a column
results$gene <- rownames(results)

# Add gene symbols
library(org.Hs.eg.db)
results$symbol <- mapIds(org.Hs.eg.db, 
                         keys=results$gene,
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")

# Look at top 20 genes
cat("\nTop 20 differentially expressed genes:\n")
print(head(results[, c("gene", "symbol", "logFC", "AveExpr", "P.Value", "adj.P.Val")], 20))

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("Total genes tested:", nrow(results), "\n")
cat("Significant genes (adj.P.Val < 0.05):", sum(results$adj.P.Val < 0.05), "\n")
cat("Significant genes (adj.P.Val < 0.05 & |logFC| > 1):", 
    sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 1), "\n")
cat("Upregulated in COVID:", sum(results$adj.P.Val < 0.05 & results$logFC > 1), "\n")
cat("Downregulated in COVID:", sum(results$adj.P.Val < 0.05 & results$logFC < -1), "\n")

# ============================================================================
# STEP 9: Volcano Plot
# ============================================================================

# Classify genes for coloring
# Categories for coloring points
# Grey = not significant (adj.o >= 0.05), Orange = significant (adj.p < 0.05 but |logFC| < 1)
# Red = upregulated in COVID (adj.p < 0.05 and logFC > 1), Blue = downregulated in COVID (adj.p < 0.05 and logFC < -1)
# X-axis = log FC (magnitude of change), Y-axis = -log10(p-value) = significance
results$significance <- "Not Significant"
results$significance[results$adj.P.Val < 0.05] <- "Significant"
results$significance[results$adj.P.Val < 0.05 & results$logFC > 1] <- "Upregulated in COVID"
results$significance[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Downregulated in COVID"

volcano_plot <- ggplot(results, aes(x=logFC, y=-log10(P.Value), color=significance)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("Upregulated in COVID"="red", 
                              "Downregulated in COVID"="blue",
                              "Significant"="orange",
                              "Not Significant"="gray")) +
  theme_minimal() +
  theme(legend.position="bottom") +
  labs(title="Volcano Plot: COVID-19 vs Influenza",
       x="Log2 Fold Change",
       y="-Log10 P-value",
       color="") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black")

ggsave("volcano_plot.pdf", volcano_plot, width=10, height=7)
print(volcano_plot)
cat("\nVolcano plot saved to 'volcano_plot.pdf'\n")

# ============================================================================
# STEP 10: MA Plot
# ============================================================================

pdf("MA_plot.pdf", width=8, height=6)
plotMA(fit2, coef=1, main="MA Plot: COVID-19 vs Influenza")
abline(h=c(-1, 1), col="blue", lty=2)
dev.off()
cat("MA plot saved to 'MA_plot.pdf'\n")

# ============================================================================
# STEP 11: Heatmap of top DEGs
# ============================================================================
# Visualize expression of top 50 genes across all samples
# Groups similar samples and genes together
# Get top 50 significant genes
top_genes <- rownames(results)[1:min(50, nrow(results))]

# Load logCPM for visualization
logcpm <- read.table("logCPM_normalized.tsv", header=TRUE, sep="\t", row.names=1)

# Make sure column names match
colnames(logcpm) <- gsub("-", ".", colnames(logcpm))
if("94189" %in% colnames(logcpm)) {
  colnames(logcpm)[colnames(logcpm) == "94189"] <- "X94189"
}

# Extract data for top genes
heatmap_data <- logcpm[top_genes, ]

# Create annotation for samples
annotation_col <- data.frame(
  Disease = metadata$cohort,
  row.names = metadata$sample
)

# Create heatmap
pdf("heatmap_top50_DEGs.pdf", width=10, height=12)
pheatmap(heatmap_data,
         scale="row",
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         show_rownames=TRUE,
         show_colnames=TRUE,
         annotation_col=annotation_col,
         main="Top 50 Differentially Expressed Genes\nCOVID-19 vs Influenza",
         fontsize_row=6)
dev.off()
cat("Heatmap saved to 'heatmap_top50_DEGs.pdf'\n")

# ============================================================================
# STEP 12: Save results
# ============================================================================

# Save full results table
write.table(results, "DE_results_COVID_vs_Influenza_full.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

# Save significant genes only (for pathway enrichment teammate)
sig_genes <- results[results$adj.P.Val < 0.05, ]
write.table(sig_genes, "significant_DEGs_for_pathway.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

# Save top 100 genes
top100 <- results[1:min(100, nrow(results)), ]
write.table(top100, "top100_DEGs.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

cat("\n=== FILES SAVED ===\n")
cat("1. DE_results_COVID_vs_Influenza_full.tsv - All results\n")
cat("2. significant_DEGs_for_pathway.tsv - Significant genes for pathway enrichment\n")
cat("3. top100_DEGs.tsv - Top 100 DEGs\n")
cat("4. volcano_plot.pdf\n")
cat("5. MA_plot.pdf\n")
cat("6. heatmap_top50_DEGs.pdf\n")
cat("7. voom_mean_variance_plot.pdf\n")

cat("\n=== ANALYSIS COMPLETE! ===\n")