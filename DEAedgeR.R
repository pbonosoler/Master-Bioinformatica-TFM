setwd("/home/pbonosoler/TFM/glycerol/counts/")

library(edgeR)
library(sva)
library(UpSetR)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)
library(pathview)
library(scales)
library(pheatmap)
library(DT)

# Load merged data
merged_data <- read.table("counts_table.tsv", header=TRUE, row.names=1, sep="\t")
# Load phenotype data
pdata <- read.table("pdata.txt", header=TRUE, row.names=1, sep=",")
# Load gene origin list
gene_origin <- read.table("Scerevisiae.geneOrigin.MAR18.txt", header=TRUE, row.names=1, sep="\t")
gene_origin$combined_origin <- ifelse(gene_origin$origin %in% c("WGD", "SSD"), "duplications", gene_origin$origin)
gene_origin$names <- rownames(gene_origin)
# Interaction counts (genes)
interaction_counts <- read.table("gene_interactions.tsv", header=FALSE, row.names=1, sep="\t")
colnames(interaction_counts) <- "nInteractions"
interaction_counts$Gene <- rownames(interaction_counts)
rownames(interaction_counts) <- NULL
# Interaction counts (protein)
p_interaction_counts <- read.table("protein_interactions.tsv", header=FALSE, row.names=1, sep="\t")
colnames(p_interaction_counts) <- "p_nInteractions"
p_interaction_counts$Gene <- rownames(p_interaction_counts)
rownames(p_interaction_counts) <- NULL

# Create DGEList object
y <- DGEList(counts=merged_data, genes=row.names(merged_data))

# Get the sample names (row names) from the DGEList object
sample_names_dge <- rownames(y$samples)

# Reorder the pdata dataframe to match the sample order in the DGEList
pdata_reordered <- pdata[match(sample_names_dge, rownames(pdata)), ]

# Add phenotypic data to DGEList
y$samples <- cbind(y$samples, pdata_reordered)

# Eliminate rows corresponding to no genes
eliminate <- grepl("_", rownames (y$counts), )
y <- y[!eliminate,]

# Get the names of all genes and their interactions
all_genes <- rownames(y$counts)
all_genes <- as.data.frame(all_genes)
all_genes_interactions <- merge(all_genes, interaction_counts, by.x="all_genes", by.y="Gene", all.x = TRUE)
all_genes_interactions$nInteractions[is.na(all_genes_interactions$nInteractions)] <- 0
all_genes_interactions <- merge(all_genes_interactions, gene_origin, by.x="all_genes", by.y = "names", all.x = TRUE)
all_genes_interactions <- merge(all_genes_interactions, p_interaction_counts, by.x="all_genes", by.y = "Gene", all.x = TRUE)
all_genes_interactions$p_nInteractions[is.na(all_genes_interactions$p_nInteractions)] <- 0
rownames(all_genes_interactions) <- all_genes_interactions$all_genes

#Group the samples
interaction_group <- interaction(pdata_reordered$line, pdata_reordered$challenge)
y$samples$group <- interaction_group

# Boxplot of raw counts
boxplot(log2(y$counts + 1), las = 2, main = "Boxplot of Raw Counts", ylab = "Log2 Counts")

# Density plot of raw counts
plot(density(log2(y$counts[, 1] + 1)), col = 1, main = "Density Plot of Raw Counts", xlab = "Log2 Counts")
for (i in 2:ncol(y$counts)) {
  lines(density(log2(y$counts[, i] + 1)), col = i)
}

# Normalization
y <- calcNormFactors(y, method = "TMM")

#Filter by expression
keep <- filterByExpr(y, min.total.count=15, min.prop=0.7, min.count=10)
y <- y[keep,]

design_all <- model.matrix(~y$samples$group)
y <- estimateDisp(y, design = design_all)

# Perform batch correction with ComBat_seq
combat_data <- ComBat_seq(y$counts, batch = y$samples$sequencing)

# Replace counts with batch corrected counts
y$counts <- combat_data

# Re-normalization
y <- calcNormFactors(y, method = "TMM")

# Boxplot of processed counts
boxplot(log2(y$counts + 1), las = 2, main = "Boxplot of Processed Counts", ylab = "Log2 Counts")

# Density plot of processed counts
plot(density(log2(y$counts[, 1] + 1)), col = 1, main = "Density Plot of Processed Counts", xlab = "Log2 Counts")
for (i in 2:ncol(y$counts)) {
  lines(density(log2(y$counts[, i] + 1)), col = i)
}

# Perform exact test
Da1 <- exactTest(y, pair = c("Da1.YPD","Da1.Glycerol"))
Da1 <- Da1$table
p.adjustBH <- p.adjust(Da1$PValue, method = "BH")
Da1 <- cbind(Da1, p.adjustBH)
keep <- (Da1$p.adjustBH < 0.05 & abs(Da1$logFC) > 1)
Da1_significant <- Da1[keep,]

Da2 <- exactTest(y, pair = c("Da2.YPD","Da2.Glycerol"))
Da2 <- Da2$table
p.adjustBH <- p.adjust(Da2$PValue, method = "BH")
Da2 <- cbind(Da2, p.adjustBH)
keep <- (Da2$p.adjustBH < 0.05 & abs(Da2$logFC) > 1)
Da2_significant <- Da2[keep,]

Da3 <- exactTest(y, pair = c("Da3.YPD","Da3.Glycerol"))
Da3 <- Da3$table
p.adjustBH <- p.adjust(Da3$PValue, method = "BH")
Da3 <- cbind(Da3, p.adjustBH)
keep <- (Da3$p.adjustBH < 0.05 & abs(Da3$logFC) > 1)
Da3_significant <- Da3[keep,]

Ga1 <- exactTest(y, pair = c("Ga1.YPD","Ga1.Glycerol"))
Ga1 <- Ga1$table
p.adjustBH <- p.adjust(Ga1$PValue, method = "BH")
Ga1 <- cbind(Ga1, p.adjustBH)
keep <- (Ga1$p.adjustBH < 0.05 & abs(Ga1$logFC) > 1)
Ga1_significant <- Ga1[keep,]

Ga2 <- exactTest(y, pair = c("Ga2.YPD","Ga2.Glycerol"))
Ga2 <- Ga2$table
p.adjustBH <- p.adjust(Ga2$PValue, method = "BH")
Ga2 <- cbind(Ga2, p.adjustBH)
keep <- (Ga2$p.adjustBH < 0.05 & abs(Ga2$logFC) > 1)
Ga2_significant <- Ga2[keep,]

Ga3 <- exactTest(y, pair = c("Ga3.YPD","Ga3.Glycerol"))
Ga3 <- Ga3$table
p.adjustBH <- p.adjust(Ga3$PValue, method = "BH")
Ga3 <- cbind(Ga3, p.adjustBH)
keep <- (Ga3$p.adjustBH < 0.05 & abs(Ga3$logFC) > 1)
Ga3_significant <- Ga3[keep,]

# Count DEGs for Da lines
n_Da1 <- nrow(Da1_significant)
n_Da2 <- nrow(Da2_significant)
n_Da3 <- nrow(Da3_significant)

# Count DEGs for Ga lines
n_Ga1 <- nrow(Ga1_significant)
n_Ga2 <- nrow(Ga2_significant)
n_Ga3 <- nrow(Ga3_significant)

# Combine the counts into groups
Da_DEGs <- c(n_Da1, n_Da2, n_Da3)
Ga_DEGs <- c(n_Ga1, n_Ga2, n_Ga3)

# Perform a Mann-Whitney U test (Wilcoxon rank-sum test) to compare the distributions
test_result_DEGs <- wilcox.test(Da_DEGs, Ga_DEGs)

#Volcano plots DEGs

plot_volcano <- function(data, title) {
  p <-ggplot(data, aes(x = logFC, y = -log10(PValue))) +
    geom_point(alpha = 0.4) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(p-value)") +
    geom_vline(xintercept = c(-1, 1), col = "red") +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    scale_y_continuous(limits = c(0, 400))
    
    print(p)
}

plot_volcano(Da1, "Da1 DEGs Volcano Plot")
plot_volcano(Da2, "Da2 DEGs Volcano Plot")
plot_volcano(Da3, "Da3 DEGs Volcano Plot")
plot_volcano(Ga1, "Ga1 DEGs Volcano Plot")
plot_volcano(Ga2, "Ga2 DEGs Volcano Plot")
plot_volcano(Ga3, "Ga3 DEGs Volcano Plot")

# Add a column for the gene identifiers (row names) in each data frame
Da1_significant$gene_id <- rownames(Da1_significant)
Da2_significant$gene_id <- rownames(Da2_significant)
Da3_significant$gene_id <- rownames(Da3_significant)
Ga1_significant$gene_id <- rownames(Ga1_significant)
Ga2_significant$gene_id <- rownames(Ga2_significant)
Ga3_significant$gene_id <- rownames(Ga3_significant)

# Add columns for the number of interactions (genes and proteins)
Da1_significant <- merge(Da1_significant, all_genes_interactions, by.x = "gene_id", by.y = "all_genes", all.x = TRUE)
Da2_significant <- merge(Da2_significant, all_genes_interactions, by.x = "gene_id", by.y = "all_genes", all.x = TRUE)
Da3_significant <- merge(Da3_significant, all_genes_interactions, by.x = "gene_id", by.y = "all_genes", all.x = TRUE)
Ga1_significant <- merge(Ga1_significant, all_genes_interactions, by.x = "gene_id", by.y = "all_genes", all.x = TRUE)
Ga2_significant <- merge(Ga2_significant, all_genes_interactions, by.x = "gene_id", by.y = "all_genes", all.x = TRUE)
Ga3_significant <- merge(Ga3_significant, all_genes_interactions, by.x = "gene_id", by.y = "all_genes", all.x = TRUE)

# Create a list of gene sets
listInput <- list(
  Da1 = Da1_significant$gene_id,
  Da2 = Da2_significant$gene_id,
  Da3 = Da3_significant$gene_id,
  Ga1 = Ga1_significant$gene_id,
  Ga2 = Ga2_significant$gene_id,
  Ga3 = Ga3_significant$gene_id
)

# Create the UpSet plot
upset(fromList(listInput),
      sets = c("Da1", "Da2", "Da3", "Ga1", "Ga2", "Ga3"),
      order.by = "freq", 
      main.bar.color = "blue",
      sets.bar.color = "gray")

# Divide the DEG into up and down-regulated
Da1_up <- Da1_significant$logFC > 1
Da1_up <- Da1_significant[Da1_up,]
Da1_down <- Da1_significant$logFC < -1
Da1_down <- Da1_significant[Da1_down,]

Da2_up <- Da2_significant$logFC > 1
Da2_up <- Da2_significant[Da2_up,]
Da2_down <- Da2_significant$logFC < -1
Da2_down <- Da2_significant[Da2_down,]

Da3_up <- Da3_significant$logFC > 1
Da3_up <- Da3_significant[Da3_up,]
Da3_down <- Da3_significant$logFC < -1
Da3_down <- Da3_significant[Da3_down,]

Ga1_up <- Ga1_significant$logFC > 1
Ga1_up <- Ga1_significant[Ga1_up,]
Ga1_down <- Ga1_significant$logFC < -1
Ga1_down <- Ga1_significant[Ga1_down,]

Ga2_up <- Ga2_significant$logFC > 1
Ga2_up <- Ga2_significant[Ga2_up,]
Ga2_down <- Ga2_significant$logFC < -1
Ga2_down <- Ga2_significant[Ga2_down,]

Ga3_up <- Ga3_significant$logFC > 1
Ga3_up <- Ga3_significant[Ga3_up,]
Ga3_down <- Ga3_significant$logFC < -1
Ga3_down <- Ga3_significant[Ga3_down,]


# Create a list of upregulated gene sets
listInput_up <- list(
  Da1_up = Da1_up$gene_id,
  Da2_up = Da2_up$gene_id,
  Da3_up = Da3_up$gene_id,
  Ga1_up = Ga1_up$gene_id,
  Ga2_up = Ga2_up$gene_id,
  Ga3_up = Ga3_up$gene_id
)

# Create the UpSet plot for upregulated genes
upset(fromList(listInput_up),
      sets = c("Da1_up", "Da2_up", "Da3_up", "Ga1_up", "Ga2_up", "Ga3_up"),
      order.by = "freq", 
      main.bar.color = "blue",
      sets.bar.color = "gray")


# Create a list of down-regulated gene sets
listInput_down <- list(
  Da1_down = Da1_down$gene_id,
  Da2_down = Da2_down$gene_id,
  Da3_down = Da3_down$gene_id,
  Ga1_down = Ga1_down$gene_id,
  Ga2_down = Ga2_down$gene_id,
  Ga3_down = Ga3_down$gene_id
)

# Create the UpSet plot for downregulated genes
upset(fromList(listInput_down),
      sets = c("Da1_down", "Da2_down", "Da3_down", "Ga1_down", "Ga2_down", "Ga3_down"),
      order.by = "freq", 
      main.bar.color = "blue",
      sets.bar.color = "gray")


# Define lines to analyze for functional annotation
all_lines <- list(Da1_significant$gene_id, Da2_significant$gene_id, Da3_significant$gene_id, 
                  Ga1_significant$gene_id, Ga2_significant$gene_id, Ga3_significant$gene_id)
names(all_lines) <- c("Da1_significant", "Da2_significant", "Da3_significant",
                  "Ga1_significant", "Ga2_significant", "Ga3_significant")


all_lines_go <- compareCluster(all_lines, fun = "enrichGO", universe = rownames(y$counts), 
                               OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = "BP", 
                               pAdjustMethod = "BH", pvalueCutoff = 0.05)

all_lines_go_s <- simplify(all_lines_go, cutoff = 0.7,  by = "p.adjust",
                           select_fun = min, measure = "Wang")

# Generate the dotplot
p <- dotplot(all_lines_go_s, showCategory = 10, label_format = 80)
# Save the plot
ggsave("dotplot_all_go.png", plot = p, width = 15, height = 10)

lines_up <- list(Da1_up$gene_id, Da2_up$gene_id, Da3_up$gene_id, Ga1_up$gene_id, Ga2_up$gene_id, Ga3_up$gene_id)
names(lines_up) <- c("Da1_up", "Da2_up", "Da3_up", "Ga1_up", "Ga2_up", "Ga3_up")

up_lines_go <- compareCluster(lines_up, fun = "enrichGO", universe = rownames(y$counts), 
                              OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = "BP", 
                              pAdjustMethod = "BH", pvalueCutoff = 0.05)

up_lines_go_s <- simplify(up_lines_go, cutoff = 0.7,  by = "p.adjust",
                          select_fun = min, measure = "Wang")

p <- dotplot(up_lines_go_s, showCategory = 10, label_format = 80)
ggsave("dotplot_up_go.png", plot = p, width = 15, height = 10)

lines_down <- list(Da1_down$gene_id, Da2_down$gene_id, Da3_down$gene_id, Ga1_down$gene_id, Ga2_down$gene_id, Ga3_down$gene_id)
names(lines_down) <- c("Da1_down", "Da2_down", "Da3_down", "Ga1_down", "Ga2_down", "Ga3_down")

down_lines_go <- compareCluster(lines_down, fun = "enrichGO", universe = rownames(y$counts), 
                              OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = "BP", 
                              pAdjustMethod = "BH", pvalueCutoff = 0.05)

down_lines_go_s <- simplify(down_lines_go, cutoff = 0.7,  by = "p.adjust",
                          select_fun = min, measure = "Wang")

p <- dotplot(down_lines_go_s, showCategory = 10, label_format = 80)
ggsave("dotplot_down_go.png", plot = p, width = 15, height = 10)

all_lines_kegg <- compareCluster(all_lines, fun = "enrichKEGG", organism = "sce", 
                              pAdjustMethod = "BH", pvalueCutoff = 0.05)
all_lines_kegg@compareClusterResult$Description <- sapply(strsplit(all_lines_kegg@compareClusterResult$Description, " - S"), `[`, 1)

p <- dotplot(all_lines_kegg, showCategory = 10, label_format = 80)
ggsave("dotplot_all_KEGG.png", plot = p, width = 15, height = 10)

up_lines_kegg <- compareCluster(lines_up, fun = "enrichKEGG", organism = "sce",
                                pAdjustMethod = "BH", pvalueCutoff = 0.05)
up_lines_kegg@compareClusterResult$Description <- sapply(strsplit(up_lines_kegg@compareClusterResult$Description, " - S"), `[`, 1)

p <- dotplot(up_lines_kegg, showCategory = 10, label_format = 80)
ggsave("dotplot_up_KEGG.png", plot = p, width = 15, height = 10)

down_lines_kegg <- compareCluster(lines_down, fun = "enrichKEGG", organism = "sce",
                                pAdjustMethod = "BH", pvalueCutoff = 0.05)
down_lines_kegg@compareClusterResult$Description <- sapply(strsplit(down_lines_kegg@compareClusterResult$Description, " - S"), `[`, 1)


p <- dotplot(down_lines_kegg, showCategory = 10, label_format = 80)
ggsave("dotplot_down_KEGG.png", plot = p, width = 15, height = 10)

# Gene origin analysis
# Calculate the proportions of duplicates in the entire yeast genome
total_counts <- table(gene_origin$combined_origin)

# Perform Fisher's exact test for significant genes
Da1_counts <- table(Da1_significant$combined_origin)
Da1_contingency_table <- rbind(Da1_counts, total_counts)
Da1_fisher_test <- fisher.test(Da1_contingency_table)

Da2_counts <- table(Da2_significant$combined_origin)
Da2_contingency_table <- rbind(Da2_counts, total_counts)
Da2_fisher_test <- fisher.test(Da2_contingency_table)

Da3_counts <- table(Da3_significant$combined_origin)
Da3_contingency_table <- rbind(Da3_counts, total_counts)
Da3_fisher_test <- fisher.test(Da3_contingency_table)

Ga1_counts <- table(Ga1_significant$combined_origin)
Ga1_contingency_table <- rbind(Ga1_counts, total_counts)
Ga1_fisher_test <- fisher.test(Ga1_contingency_table)

Ga2_counts <- table(Ga2_significant$combined_origin)
Ga2_contingency_table <- rbind(Ga2_counts, total_counts)
Ga2_fisher_test <- fisher.test(Ga2_contingency_table)

Ga3_counts <- table(Ga3_significant$combined_origin)
Ga3_contingency_table <- rbind(Ga3_counts, total_counts)
Ga3_fisher_test <- fisher.test(Ga3_contingency_table)

# Extract odds ratios (estimates) and p-values
odds_ratios <- c(1,
                 Da1_fisher_test$estimate, 
                 Da2_fisher_test$estimate, 
                 Da3_fisher_test$estimate, 
                 Ga1_fisher_test$estimate, 
                 Ga2_fisher_test$estimate, 
                 Ga3_fisher_test$estimate)

p_values <- c(1, 
              Da1_fisher_test$p.value, 
              Da2_fisher_test$p.value, 
              Da3_fisher_test$p.value, 
              Ga1_fisher_test$p.value, 
              Ga2_fisher_test$p.value, 
              Ga3_fisher_test$p.value)

# Create a data frame for ggplot
plot_data <- data.frame(
  Set = c("Total", "Da1", "Da2", "Da3", "Ga1", "Ga2", "Ga3"),
  OddsRatio = unlist(odds_ratios),
  p_value = p_values
)

# Plot the data using ggplot2
ggplot(plot_data, aes(x = Set, y = OddsRatio)) +
  geom_bar(stat = "identity", fill = "grey", width = 0.8, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = ifelse(!is.na(p_value), paste("p =", scales::pvalue(p_value)), ""), 
                group = Set), 
            vjust = -0.5,  
            size = 5) +  
  labs(title = "Odds Ratios of Significant Genes by Origin",
       x = "Sample Set",
       y = "Odds Ratio Duplicates vs Singletons") +
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Perform Fisher's exact test for up-regulated genes
Da1_up_counts <- table(Da1_up$combined_origin)
Da1_up_contingency_table <- rbind(Da1_up_counts, total_counts)
Da1_up_fisher_test <- fisher.test(Da1_up_contingency_table)

Da2_up_counts <- table(Da2_up$combined_origin)
Da2_up_contingency_table <- rbind(Da2_up_counts, total_counts)
Da2_up_fisher_test <- fisher.test(Da2_up_contingency_table)

Da3_up_counts <- table(Da3_up$combined_origin)
Da3_up_contingency_table <- rbind(Da3_up_counts, total_counts)
Da3_up_fisher_test <- fisher.test(Da3_up_contingency_table)

Ga1_up_counts <- table(Ga1_up$combined_origin)
Ga1_up_contingency_table <- rbind(Ga1_up_counts, total_counts)
Ga1_up_fisher_test <- fisher.test(Ga1_up_contingency_table)

Ga2_up_counts <- table(Ga2_up$combined_origin)
Ga2_up_contingency_table <- rbind(Ga2_up_counts, total_counts)
Ga2_up_fisher_test <- fisher.test(Ga2_up_contingency_table)

Ga3_up_counts <- table(Ga3_up$combined_origin)
Ga3_up_contingency_table <- rbind(Ga3_up_counts, total_counts)
Ga3_up_fisher_test <- fisher.test(Ga3_up_contingency_table)

# Extract odds ratios (estimates) and p-values
odds_ratios <- c(1,
                 Da1_up_fisher_test$estimate, 
                 Da2_up_fisher_test$estimate, 
                 Da3_up_fisher_test$estimate, 
                 Ga1_up_fisher_test$estimate, 
                 Ga2_up_fisher_test$estimate, 
                 Ga3_up_fisher_test$estimate)

p_values <- c(1, 
              Da1_up_fisher_test$p.value, 
              Da2_up_fisher_test$p.value, 
              Da3_up_fisher_test$p.value, 
              Ga1_up_fisher_test$p.value, 
              Ga2_up_fisher_test$p.value, 
              Ga3_up_fisher_test$p.value)

# Create a data frame for ggplot
plot_data <- data.frame(
  Set = c("Total", "Da1_up", "Da2_up", "Da3_up", "Ga1_up", "Ga2_up", "Ga3_up"),
  OddsRatio = unlist(odds_ratios),
  p_value = p_values
)

# Plot the data using ggplot2
ggplot(plot_data, aes(x = Set, y = OddsRatio)) +
  geom_bar(stat = "identity", fill = "grey", width = 0.8, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = ifelse(!is.na(p_value), paste("p =", scales::pvalue(p_value)), ""), 
                group = Set), 
            vjust = -0.5,  
            size = 5) +  
  labs(title = "Odds Ratios of Upregulated Genes by Origin",
       x = "Sample Set",
       y = "Odds Ratio Duplicates vs Singletons") +
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Perform Fisher's exact test for down-regulated genes
Da1_down_counts <- table(Da1_down$combined_origin)
Da1_down_contingency_table <- rbind(Da1_down_counts, total_counts)
Da1_down_fisher_test <- fisher.test(Da1_down_contingency_table)

Da2_down_counts <- table(Da2_down$combined_origin)
Da2_down_contingency_table <- rbind(Da2_down_counts, total_counts)
Da2_down_fisher_test <- fisher.test(Da2_down_contingency_table)

Da3_down_counts <- table(Da3_down$combined_origin)
Da3_down_contingency_table <- rbind(Da3_down_counts, total_counts)
Da3_down_fisher_test <- fisher.test(Da3_down_contingency_table)

Ga1_down_counts <- table(Ga1_down$combined_origin)
Ga1_down_contingency_table <- rbind(Ga1_down_counts, total_counts)
Ga1_down_fisher_test <- fisher.test(Ga1_down_contingency_table)

Ga2_down_counts <- table(Ga2_down$combined_origin)
Ga2_down_contingency_table <- rbind(Ga2_down_counts, total_counts)
Ga2_down_fisher_test <- fisher.test(Ga2_down_contingency_table)

Ga3_down_counts <- table(Ga3_down$combined_origin)
Ga3_down_contingency_table <- rbind(Ga3_down_counts, total_counts)
Ga3_down_fisher_test <- fisher.test(Ga3_down_contingency_table)

# Extract odds ratios (estimates) and p-values
odds_ratios <- c(1,
                 Da1_down_fisher_test$estimate, 
                 Da2_down_fisher_test$estimate, 
                 Da3_down_fisher_test$estimate, 
                 Ga1_down_fisher_test$estimate, 
                 Ga2_down_fisher_test$estimate, 
                 Ga3_down_fisher_test$estimate)

p_values <- c(1, 
              Da1_down_fisher_test$p.value, 
              Da2_down_fisher_test$p.value, 
              Da3_down_fisher_test$p.value, 
              Ga1_down_fisher_test$p.value, 
              Ga2_down_fisher_test$p.value, 
              Ga3_down_fisher_test$p.value)

# Create a data frame for ggplot
plot_data <- data.frame(
  Set = c("Total", "Da1_down", "Da2_down", "Da3_down", "Ga1_down", "Ga2_down", "Ga3_down"),
  OddsRatio = unlist(odds_ratios),
  p_value = p_values
)

# Plot the data using ggplot2
ggplot(plot_data, aes(x = Set, y = OddsRatio)) +
  geom_bar(stat = "identity", fill = "grey", width = 0.8, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = ifelse(!is.na(p_value), paste("p =", scales::pvalue(p_value)), ""), 
                group = Set), 
            vjust = -0.5,  
            size = 5) +  
  labs(title = "Odds Ratios of Downregulated Genes by Origin",
       x = "Sample Set",
       y = "Odds Ratio Duplicates vs Singletons") +
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Create contingency tables and perform Fisher's exact test for SSD vs WGD
total_counts_2 <- table(gene_origin$origin)

Da1_wgd_ssd_counts <- table(Da1_significant$origin)
Da1_wgd_ssd_counts <- Da1_wgd_ssd_counts[c("WGD", "SSD")]
Da1_wgd_ssd_table <- rbind(Da1_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da1_wgd_ssd_fisher <- fisher.test(Da1_wgd_ssd_table)

Da2_wgd_ssd_counts <- table(Da2_significant$origin)
Da2_wgd_ssd_counts <- Da2_wgd_ssd_counts[c("WGD", "SSD")]
Da2_wgd_ssd_table <- rbind(Da2_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da2_wgd_ssd_fisher <- fisher.test(Da2_wgd_ssd_table)

Da3_wgd_ssd_counts <- table(Da3_significant$origin)
Da3_wgd_ssd_counts <- Da3_wgd_ssd_counts[c("WGD", "SSD")]
Da3_wgd_ssd_table <- rbind(Da3_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da3_wgd_ssd_fisher <- fisher.test(Da3_wgd_ssd_table)

Ga1_wgd_ssd_counts <- table(Ga1_significant$origin)
Ga1_wgd_ssd_counts <- Ga1_wgd_ssd_counts[c("WGD", "SSD")]
Ga1_wgd_ssd_table <- rbind(Ga1_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga1_wgd_ssd_fisher <- fisher.test(Ga1_wgd_ssd_table)

Ga2_wgd_ssd_counts <- table(Ga2_significant$origin)
Ga2_wgd_ssd_counts <- Ga2_wgd_ssd_counts[c("WGD", "SSD")]
Ga2_wgd_ssd_table <- rbind(Ga2_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga2_wgd_ssd_fisher <- fisher.test(Ga2_wgd_ssd_table)

Ga3_wgd_ssd_counts <- table(Ga3_significant$origin)
Ga3_wgd_ssd_counts <- Ga3_wgd_ssd_counts[c("WGD", "SSD")]
Ga3_wgd_ssd_table <- rbind(Ga3_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga3_wgd_ssd_fisher <- fisher.test(Ga3_wgd_ssd_table)

# Extract odds ratios (estimates) and p-values
odds_ratios <- c(1,
                 Da1_wgd_ssd_fisher$estimate, 
                 Da2_wgd_ssd_fisher$estimate, 
                 Da3_wgd_ssd_fisher$estimate, 
                 Ga1_wgd_ssd_fisher$estimate, 
                 Ga2_wgd_ssd_fisher$estimate, 
                 Ga3_wgd_ssd_fisher$estimate)

p_values <- c(1, 
              Da1_wgd_ssd_fisher$p.value, 
              Da2_wgd_ssd_fisher$p.value, 
              Da3_wgd_ssd_fisher$p.value, 
              Ga1_wgd_ssd_fisher$p.value, 
              Ga2_wgd_ssd_fisher$p.value, 
              Ga3_wgd_ssd_fisher$p.value)

# Create a data frame for ggplot
plot_data <- data.frame(
  Set = c("Total", "Da1", "Da2", "Da3", "Ga1", "Ga2", "Ga3"),
  OddsRatio = unlist(odds_ratios),
  p_value = p_values
)

# Plot the data using ggplot2
ggplot(plot_data, aes(x = Set, y = OddsRatio)) +
  geom_bar(stat = "identity", fill = "grey", width = 0.8, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = ifelse(!is.na(p_value), paste("p =", scales::pvalue(p_value)), ""), 
                group = Set), 
            vjust = -0.5,  
            size = 5) +  
  labs(title = "Odds Ratios of Significant Duplicated Genes",
       x = "Sample Set",
       y = "Odds Ratio WGDs vs SSDs") +
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Create contingency tables and perform Fisher's exact test for up-regulated genes
Da1_up_wgd_ssd_counts <- table(Da1_up$origin)
Da1_up_wgd_ssd_counts <- Da1_up_wgd_ssd_counts[c("WGD", "SSD")]
Da1_up_wgd_ssd_table <- rbind(Da1_up_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da1_up_wgd_ssd_fisher <- fisher.test(Da1_up_wgd_ssd_table)

Da2_up_wgd_ssd_counts <- table(Da2_up$origin)
Da2_up_wgd_ssd_counts <- Da2_up_wgd_ssd_counts[c("WGD", "SSD")]
Da2_up_wgd_ssd_table <- rbind(Da2_up_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da2_up_wgd_ssd_fisher <- fisher.test(Da2_up_wgd_ssd_table)

Da3_up_wgd_ssd_counts <- table(Da3_up$origin)
Da3_up_wgd_ssd_counts <- Da3_up_wgd_ssd_counts[c("WGD", "SSD")]
Da3_up_wgd_ssd_table <- rbind(Da3_up_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da3_up_wgd_ssd_fisher <- fisher.test(Da3_up_wgd_ssd_table)

Ga1_up_wgd_ssd_counts <- table(Ga1_up$origin)
Ga1_up_wgd_ssd_counts <- Ga1_up_wgd_ssd_counts[c("WGD", "SSD")]
Ga1_up_wgd_ssd_table <- rbind(Ga1_up_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga1_up_wgd_ssd_fisher <- fisher.test(Ga1_up_wgd_ssd_table)

Ga2_up_wgd_ssd_counts <- table(Ga2_up$origin)
Ga2_up_wgd_ssd_counts <- Ga2_up_wgd_ssd_counts[c("WGD", "SSD")]
Ga2_up_wgd_ssd_table <- rbind(Ga2_up_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga2_up_wgd_ssd_fisher <- fisher.test(Ga2_up_wgd_ssd_table)

Ga3_up_wgd_ssd_counts <- table(Ga3_up$origin)
Ga3_up_wgd_ssd_counts <- Ga3_up_wgd_ssd_counts[c("WGD", "SSD")]
Ga3_up_wgd_ssd_table <- rbind(Ga3_up_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga3_up_wgd_ssd_fisher <- fisher.test(Ga3_up_wgd_ssd_table)

# Extract odds ratios (estimates) and p-values
odds_ratios <- c(1,
                 Da1_up_wgd_ssd_fisher$estimate, 
                 Da2_up_wgd_ssd_fisher$estimate, 
                 Da3_up_wgd_ssd_fisher$estimate, 
                 Ga1_up_wgd_ssd_fisher$estimate, 
                 Ga2_up_wgd_ssd_fisher$estimate, 
                 Ga3_up_wgd_ssd_fisher$estimate)

p_values <- c(1, 
              Da1_up_wgd_ssd_fisher$p.value, 
              Da2_up_wgd_ssd_fisher$p.value, 
              Da3_up_wgd_ssd_fisher$p.value, 
              Ga1_up_wgd_ssd_fisher$p.value, 
              Ga2_up_wgd_ssd_fisher$p.value, 
              Ga3_up_wgd_ssd_fisher$p.value)

# Create a data frame for ggplot
plot_data <- data.frame(
  Set = c("Total", "Da1_up", "Da2_up", "Da3_up", "Ga1_up", "Ga2_up", "Ga3_up"),
  OddsRatio = unlist(odds_ratios),
  p_value = p_values
)

# Plot the data using ggplot2
ggplot(plot_data, aes(x = Set, y = OddsRatio)) +
  geom_bar(stat = "identity", fill = "grey", width = 0.8, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = ifelse(!is.na(p_value), paste("p =", scales::pvalue(p_value)), ""), 
                group = Set), 
            vjust = -0.5,  
            size = 5) +  
  labs(title = "Odds Ratios of Upregulated Duplicated Genes",
       x = "Sample Set",
       y = "Odds Ratio WGDs vs SSDs") +
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Create contingency tables and perform Fisher's exact test for down-regulated genes
Da1_down_wgd_ssd_counts <- table(Da1_down$origin)
Da1_down_wgd_ssd_counts <- Da1_down_wgd_ssd_counts[c("WGD", "SSD")]
Da1_down_wgd_ssd_table <- rbind(Da1_down_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da1_down_wgd_ssd_fisher <- fisher.test(Da1_down_wgd_ssd_table)

Da2_down_wgd_ssd_counts <- table(Da2_down$origin)
Da2_down_wgd_ssd_counts <- Da2_down_wgd_ssd_counts[c("WGD", "SSD")]
Da2_down_wgd_ssd_table <- rbind(Da2_down_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da2_down_wgd_ssd_fisher <- fisher.test(Da2_down_wgd_ssd_table)

Da3_down_wgd_ssd_counts <- table(Da3_down$origin)
Da3_down_wgd_ssd_counts <- Da3_down_wgd_ssd_counts[c("WGD", "SSD")]
Da3_down_wgd_ssd_table <- rbind(Da3_down_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Da3_down_wgd_ssd_fisher <- fisher.test(Da3_down_wgd_ssd_table)

Ga1_down_wgd_ssd_counts <- table(Ga1_down$origin)
Ga1_down_wgd_ssd_counts <- Ga1_down_wgd_ssd_counts[c("WGD", "SSD")]
Ga1_down_wgd_ssd_table <- rbind(Ga1_down_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga1_down_wgd_ssd_fisher <- fisher.test(Ga1_down_wgd_ssd_table)

Ga2_down_wgd_ssd_counts <- table(Ga2_down$origin)
Ga2_down_wgd_ssd_counts <- Ga2_down_wgd_ssd_counts[c("WGD", "SSD")]
Ga2_down_wgd_ssd_table <- rbind(Ga2_down_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga2_down_wgd_ssd_fisher <- fisher.test(Ga2_down_wgd_ssd_table)

Ga3_down_wgd_ssd_counts <- table(Ga3_down$origin)
Ga3_down_wgd_ssd_counts <- Ga3_down_wgd_ssd_counts[c("WGD", "SSD")]
Ga3_down_wgd_ssd_table <- rbind(Ga3_down_wgd_ssd_counts, total_counts_2[c("WGD", "SSD")])
Ga3_down_wgd_ssd_fisher <- fisher.test(Ga3_down_wgd_ssd_table)

# Extract odds ratios (estimates) and p-values
odds_ratios <- c(1,
                 Da1_down_wgd_ssd_fisher$estimate, 
                 Da2_down_wgd_ssd_fisher$estimate, 
                 Da3_down_wgd_ssd_fisher$estimate, 
                 Ga1_down_wgd_ssd_fisher$estimate, 
                 Ga2_down_wgd_ssd_fisher$estimate, 
                 Ga3_down_wgd_ssd_fisher$estimate)

p_values <- c(1, 
              Da1_down_wgd_ssd_fisher$p.value, 
              Da2_down_wgd_ssd_fisher$p.value, 
              Da3_down_wgd_ssd_fisher$p.value, 
              Ga1_down_wgd_ssd_fisher$p.value, 
              Ga2_down_wgd_ssd_fisher$p.value, 
              Ga3_down_wgd_ssd_fisher$p.value)

# Create a data frame for ggplot
plot_data <- data.frame(
  Set = c("Total", "Da1_down", "Da2_down", "Da3_down", "Ga1_down", "Ga2_down", "Ga3_down"),
  OddsRatio = unlist(odds_ratios),
  p_value = p_values
)

# Plot the data using ggplot2
ggplot(plot_data, aes(x = Set, y = OddsRatio)) +
  geom_bar(stat = "identity", fill = "grey", width = 0.8, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = ifelse(!is.na(p_value), paste("p =", scales::pvalue(p_value)), ""), 
                group = Set), 
            vjust = -0.5,  
            size = 5) +  
  labs(title = "Odds Ratios of Downregulated Duplicated Genes",
       x = "Sample Set",
       y = "Odds Ratio WGDs vs SSDs") +
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Perform Wilcoxon rank-sum test to compare log fold changes (logFC) between duplicates and singletons
Da1_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da1_significant, subset = combined_origin %in% c("singleton", "duplications"))
Da2_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da2_significant, subset = combined_origin %in% c("singleton", "duplications"))
Da3_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da3_significant, subset = combined_origin %in% c("singleton", "duplications"))
Ga1_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga1_significant, subset = combined_origin %in% c("singleton", "duplications"))
Ga2_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga2_significant, subset = combined_origin %in% c("singleton", "duplications"))
Ga3_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga3_significant, subset = combined_origin %in% c("singleton", "duplications"))
Da1_up_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da1_up, subset = combined_origin %in% c("singleton", "duplications"))
Da1_down_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da1_down, subset = combined_origin %in% c("singleton", "duplications"))
Da2_up_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da2_up, subset = combined_origin %in% c("singleton", "duplications"))
Da2_down_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da2_down, subset = combined_origin %in% c("singleton", "duplications"))
Da3_up_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da3_up, subset = combined_origin %in% c("singleton", "duplications"))
Da3_down_wilcox <- wilcox.test(logFC ~ combined_origin, data = Da3_down, subset = combined_origin %in% c("singleton", "duplications"))
Ga1_up_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga1_up, subset = combined_origin %in% c("singleton", "duplications"))
Ga1_down_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga1_down, subset = combined_origin %in% c("singleton", "duplications"))
Ga2_up_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga2_up, subset = combined_origin %in% c("singleton", "duplications"))
Ga2_down_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga2_down, subset = combined_origin %in% c("singleton", "duplications"))
Ga3_up_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga3_up, subset = combined_origin %in% c("singleton", "duplications"))
Ga3_down_wilcox <- wilcox.test(logFC ~ combined_origin, data = Ga3_down, subset = combined_origin %in% c("singleton", "duplications"))

# Perform Wilcoxon rank-sum test to compare log fold changes (logFC) between WGDs and SSDs
Da1_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da1_significant, subset = origin %in% c("WGD", "SSD"))
Da2_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da2_significant, subset = origin %in% c("WGD", "SSD"))
Da3_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da3_significant, subset = origin %in% c("WGD", "SSD"))
Ga1_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga1_significant, subset = origin %in% c("WGD", "SSD"))
Ga2_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga2_significant, subset = origin %in% c("WGD", "SSD"))
Ga3_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga3_significant, subset = origin %in% c("WGD", "SSD"))

Da1_up_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da1_up, subset = origin %in% c("WGD", "SSD"))
Da1_down_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da1_down, subset = origin %in% c("WGD", "SSD"))
Da2_up_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da2_up, subset = origin %in% c("WGD", "SSD"))
Da2_down_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da2_down, subset = origin %in% c("WGD", "SSD"))
Da3_up_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da3_up, subset = origin %in% c("WGD", "SSD"))
Da3_down_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Da3_down, subset = origin %in% c("WGD", "SSD"))

Ga1_up_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga1_up, subset = origin %in% c("WGD", "SSD"))
Ga1_down_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga1_down, subset = origin %in% c("WGD", "SSD"))
Ga2_up_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga2_up, subset = origin %in% c("WGD", "SSD"))
Ga2_down_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga2_down, subset = origin %in% c("WGD", "SSD"))
Ga3_up_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga3_up, subset = origin %in% c("WGD", "SSD"))
Ga3_down_wgd_ssd_wilcox <- wilcox.test(logFC ~ origin, data = Ga3_down, subset = origin %in% c("WGD", "SSD"))

# Combine the results into summary tables
# Create a summary data frame for singletons vs duplicates
results_duplicates_singletons <- data.frame(
  Line = c("Da1", "Da2", "Da3", "Ga1", "Ga2", "Ga3",
           "Da1_up", "Da2_up", "Da3_up", "Ga1_up", "Ga2_up", "Ga3_up",
           "Da1_down", "Da2_down", "Da3_down", "Ga1_down", "Ga2_down", "Ga3_down"),
  Fisher_p.value = c(Da1_fisher_test$p.value, Da2_fisher_test$p.value, Da3_fisher_test$p.value,
                     Ga1_fisher_test$p.value, Ga2_fisher_test$p.value, Ga3_fisher_test$p.value,
                     Da1_up_fisher_test$p.value, Da2_up_fisher_test$p.value, Da3_up_fisher_test$p.value,
                     Ga1_up_fisher_test$p.value, Ga2_up_fisher_test$p.value, Ga3_up_fisher_test$p.value,
                     Da1_down_fisher_test$p.value, Da2_down_fisher_test$p.value, Da3_down_fisher_test$p.value,
                     Ga1_down_fisher_test$p.value, Ga2_down_fisher_test$p.value, Ga3_down_fisher_test$p.value),
  
  Odds.Ratio = c(Da1_fisher_test$estimate, Da2_fisher_test$estimate, Da3_fisher_test$estimate,
                 Ga1_fisher_test$estimate, Ga2_fisher_test$estimate, Ga3_fisher_test$estimate,
                 Da1_up_fisher_test$estimate, Da2_up_fisher_test$estimate, Da3_up_fisher_test$estimate,
                 Ga1_up_fisher_test$estimate, Ga2_up_fisher_test$estimate, Ga3_up_fisher_test$estimate,
                 Da1_down_fisher_test$estimate, Da2_down_fisher_test$estimate, Da3_down_fisher_test$estimate,
                 Ga1_down_fisher_test$estimate, Ga2_down_fisher_test$estimate, Ga3_down_fisher_test$estimate),
  
  Wilcox_p.value = c(Da1_wilcox$p.value, Da2_wilcox$p.value, Da3_wilcox$p.value,
                     Ga1_wilcox$p.value, Ga2_wilcox$p.value, Ga3_wilcox$p.value,
                     Da1_up_wilcox$p.value, Da2_up_wilcox$p.value, Da3_up_wilcox$p.value,
                     Ga1_up_wilcox$p.value, Ga2_up_wilcox$p.value, Ga3_up_wilcox$p.value,
                     Da1_down_wilcox$p.value, Da2_down_wilcox$p.value, Da3_down_wilcox$p.value,
                     Ga1_down_wilcox$p.value, Ga2_down_wilcox$p.value, Ga3_down_wilcox$p.value)
)

# Create a summary data frame for WGDs vs. SSDs
results_wgd_vs_ssd <- data.frame(
  Line = c("Da1", "Da2", "Da3", "Ga1", "Ga2", "Ga3",
           "Da1_up", "Da2_up", "Da3_up", "Ga1_up", "Ga2_up", "Ga3_up",
           "Da1_down", "Da2_down", "Da3_down", "Ga1_down", "Ga2_down", "Ga3_down"),
  
  Fisher_p.value = c(Da1_wgd_ssd_fisher$p.value, Da2_wgd_ssd_fisher$p.value, Da3_wgd_ssd_fisher$p.value,
                     Ga1_wgd_ssd_fisher$p.value, Ga2_wgd_ssd_fisher$p.value, Ga3_wgd_ssd_fisher$p.value,
                     Da1_up_wgd_ssd_fisher$p.value, Da2_up_wgd_ssd_fisher$p.value, Da3_up_wgd_ssd_fisher$p.value,
                     Ga1_up_wgd_ssd_fisher$p.value, Ga2_up_wgd_ssd_fisher$p.value, Ga3_up_wgd_ssd_fisher$p.value,
                     Da1_down_wgd_ssd_fisher$p.value, Da2_down_wgd_ssd_fisher$p.value, Da3_down_wgd_ssd_fisher$p.value,
                     Ga1_down_wgd_ssd_fisher$p.value, Ga2_down_wgd_ssd_fisher$p.value, Ga3_down_wgd_ssd_fisher$p.value),
  
  Odds.Ratio = c(Da1_wgd_ssd_fisher$estimate, Da2_wgd_ssd_fisher$estimate, Da3_wgd_ssd_fisher$estimate,
                 Ga1_wgd_ssd_fisher$estimate, Ga2_wgd_ssd_fisher$estimate, Ga3_wgd_ssd_fisher$estimate,
                 Da1_up_wgd_ssd_fisher$estimate, Da2_up_wgd_ssd_fisher$estimate, Da3_up_wgd_ssd_fisher$estimate,
                 Ga1_up_wgd_ssd_fisher$estimate, Ga2_up_wgd_ssd_fisher$estimate, Ga3_up_wgd_ssd_fisher$estimate,
                 Da1_down_wgd_ssd_fisher$estimate, Da2_down_wgd_ssd_fisher$estimate, Da3_down_wgd_ssd_fisher$estimate,
                 Ga1_down_wgd_ssd_fisher$estimate, Ga2_down_wgd_ssd_fisher$estimate, Ga3_down_wgd_ssd_fisher$estimate),
  
  Wilcox_p.value = c(Da1_wgd_ssd_wilcox$p.value, Da2_wgd_ssd_wilcox$p.value, Da3_wgd_ssd_wilcox$p.value,
                     Ga1_wgd_ssd_wilcox$p.value, Ga2_wgd_ssd_wilcox$p.value, Ga3_wgd_ssd_wilcox$p.value,
                     Da1_up_wgd_ssd_wilcox$p.value, Da2_up_wgd_ssd_wilcox$p.value, Da3_up_wgd_ssd_wilcox$p.value,
                     Ga1_up_wgd_ssd_wilcox$p.value, Ga2_up_wgd_ssd_wilcox$p.value, Ga3_up_wgd_ssd_wilcox$p.value,
                     Da1_down_wgd_ssd_wilcox$p.value, Da2_down_wgd_ssd_wilcox$p.value, Da3_down_wgd_ssd_wilcox$p.value,
                     Ga1_down_wgd_ssd_wilcox$p.value, Ga2_down_wgd_ssd_wilcox$p.value, Ga3_down_wgd_ssd_wilcox$p.value)
)


# Analyze the number of interactions
# Genetic Interactions
# Calculate means and standard deviations for each group
# Define group names
names <- c("Da1_up",  "Da2_up",  "Da3_up",  "Ga1_up",  "Ga2_up",  "Ga3_up",
  "Da1_down",  "Da2_down",  "Da3_down",  "Ga1_down",  "Ga2_down",  "Ga3_down",
  "Da1_up_duplicates",  "Da1_up_singletons",  "Da2_up_duplicates",
  "Da2_up_singletons",  "Da3_up_duplicates",  "Da3_up_singletons",
  "Ga1_up_duplicates",  "Ga1_up_singletons",  "Ga2_up_duplicates",
  "Ga2_up_singletons",  "Ga3_up_duplicates",  "Ga3_up_singletons",
  "Da1_down_duplicates",  "Da1_down_singletons",  "Da2_down_duplicates",
  "Da2_down_singletons",  "Da3_down_duplicates",  "Da3_down_singletons",
  "Ga1_down_duplicates",  "Ga1_down_singletons",  "Ga2_down_duplicates",
  "Ga2_down_singletons",  "Ga3_down_duplicates",  "Ga3_down_singletons"
)
means <- c(
  Da1_up_mean = mean(Da1_up$nInteractions, na.rm = TRUE),
  Da2_up_mean = mean(Da2_up$nInteractions, na.rm = TRUE),
  Da3_up_mean = mean(Da3_up$nInteractions, na.rm = TRUE),
  Ga1_up_mean = mean(Ga1_up$nInteractions, na.rm = TRUE),
  Ga2_up_mean = mean(Ga2_up$nInteractions, na.rm = TRUE),
  Ga3_up_mean = mean(Ga3_up$nInteractions, na.rm = TRUE),
  Da1_down_mean = mean(Da1_down$nInteractions, na.rm = TRUE),
  Da2_down_mean = mean(Da2_down$nInteractions, na.rm = TRUE),
  Da3_down_mean = mean(Da3_down$nInteractions, na.rm = TRUE),
  Ga1_down_mean = mean(Ga1_down$nInteractions, na.rm = TRUE),
  Ga2_down_mean = mean(Ga2_down$nInteractions, na.rm = TRUE),
  Ga3_down_mean = mean(Ga3_down$nInteractions, na.rm = TRUE),
  Da1_up_duplicates_mean = mean(Da1_up$nInteractions[Da1_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da1_up_singletons_mean = mean(Da1_up$nInteractions[Da1_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da2_up_duplicates_mean = mean(Da2_up$nInteractions[Da2_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da2_up_singletons_mean = mean(Da2_up$nInteractions[Da2_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da3_up_duplicates_mean = mean(Da3_up$nInteractions[Da3_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da3_up_singletons_mean = mean(Da3_up$nInteractions[Da3_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga1_up_duplicates_mean = mean(Ga1_up$nInteractions[Ga1_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga1_up_singletons_mean = mean(Ga1_up$nInteractions[Ga1_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga2_up_duplicates_mean = mean(Ga2_up$nInteractions[Ga2_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga2_up_singletons_mean = mean(Ga2_up$nInteractions[Ga2_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga3_up_duplicates_mean = mean(Ga3_up$nInteractions[Ga3_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga3_up_singletons_mean = mean(Ga3_up$nInteractions[Ga3_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da1_down_duplicates_mean = mean(Da1_down$nInteractions[Da1_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da1_down_singletons_mean = mean(Da1_down$nInteractions[Da1_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da2_down_duplicates_mean = mean(Da2_down$nInteractions[Da2_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da2_down_singletons_mean = mean(Da2_down$nInteractions[Da2_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da3_down_duplicates_mean = mean(Da3_down$nInteractions[Da3_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da3_down_singletons_mean = mean(Da3_down$nInteractions[Da3_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga1_down_duplicates_mean = mean(Ga1_down$nInteractions[Ga1_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga1_down_singletons_mean = mean(Ga1_down$nInteractions[Ga1_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga2_down_duplicates_mean = mean(Ga2_down$nInteractions[Ga2_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga2_down_singletons_mean = mean(Ga2_down$nInteractions[Ga2_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga3_down_duplicates_mean = mean(Ga3_down$nInteractions[Ga3_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga3_down_singletons_mean = mean(Ga3_down$nInteractions[Ga3_down$combined_origin %in% c("singleton")], na.rm = TRUE)
)

sds <- c(
  Da1_up_sd = sd(Da1_up$nInteractions, na.rm = TRUE),
  Da2_up_sd = sd(Da2_up$nInteractions, na.rm = TRUE),
  Da3_up_sd = sd(Da3_up$nInteractions, na.rm = TRUE),
  Ga1_up_sd = sd(Ga1_up$nInteractions, na.rm = TRUE),
  Ga2_up_sd = sd(Ga2_up$nInteractions, na.rm = TRUE),
  Ga3_up_sd = sd(Ga3_up$nInteractions, na.rm = TRUE),
  Da1_down_sd = sd(Da1_down$nInteractions, na.rm = TRUE),
  Da2_down_sd = sd(Da2_down$nInteractions, na.rm = TRUE),
  Da3_down_sd = sd(Da3_down$nInteractions, na.rm = TRUE),
  Ga1_down_sd = sd(Ga1_down$nInteractions, na.rm = TRUE),
  Ga2_down_sd = sd(Ga2_down$nInteractions, na.rm = TRUE),
  Ga3_down_sd = sd(Ga3_down$nInteractions, na.rm = TRUE),
  Da1_up_duplicates_sd = sd(Da1_up$nInteractions[Da1_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da1_up_singletons_sd = sd(Da1_up$nInteractions[Da1_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da2_up_duplicates_sd = sd(Da2_up$nInteractions[Da2_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da2_up_singletons_sd = sd(Da2_up$nInteractions[Da2_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da3_up_duplicates_sd = sd(Da3_up$nInteractions[Da3_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da3_up_singletons_sd = sd(Da3_up$nInteractions[Da3_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga1_up_duplicates_sd = sd(Ga1_up$nInteractions[Ga1_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga1_up_singletons_sd = sd(Ga1_up$nInteractions[Ga1_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga2_up_duplicates_sd = sd(Ga2_up$nInteractions[Ga2_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga2_up_singletons_sd = sd(Ga2_up$nInteractions[Ga2_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga3_up_duplicates_sd = sd(Ga3_up$nInteractions[Ga3_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga3_up_singletons_sd = sd(Ga3_up$nInteractions[Ga3_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da1_down_duplicates_sd = sd(Da1_down$nInteractions[Da1_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da1_down_singletons_sd = sd(Da1_down$nInteractions[Da1_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da2_down_duplicates_sd = sd(Da2_down$nInteractions[Da2_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da2_down_singletons_sd = sd(Da2_down$nInteractions[Da2_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da3_down_duplicates_sd = sd(Da3_down$nInteractions[Da3_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da3_down_singletons_sd = sd(Da3_down$nInteractions[Da3_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga1_down_duplicates_sd = sd(Ga1_down$nInteractions[Ga1_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga1_down_singletons_sd = sd(Ga1_down$nInteractions[Ga1_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga2_down_duplicates_sd = sd(Ga2_down$nInteractions[Ga2_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga2_down_singletons_sd = sd(Ga2_down$nInteractions[Ga2_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga3_down_duplicates_sd = sd(Ga3_down$nInteractions[Ga3_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga3_down_singletons_sd = sd(Ga3_down$nInteractions[Ga3_down$combined_origin %in% c("singleton")], na.rm = TRUE)
)

size <- c(
  Da1_up_length = length(Da1_up$nInteractions),
  Da2_up_length = length(Da2_up$nInteractions),
  Da3_up_length = length(Da3_up$nInteractions),
  Ga1_up_length = length(Ga1_up$nInteractions),
  Ga2_up_length = length(Ga2_up$nInteractions),
  Ga3_up_length = length(Ga3_up$nInteractions),
  Da1_down_length = length(Da1_down$nInteractions),
  Da2_down_length = length(Da2_down$nInteractions),
  Da3_down_length = length(Da3_down$nInteractions),
  Ga1_down_length = length(Ga1_down$nInteractions),
  Ga2_down_length = length(Ga2_down$nInteractions),
  Ga3_down_length = length(Ga3_down$nInteractions),
  Da1_up_duplicates_length = length(Da1_up$nInteractions[Da1_up$combined_origin %in% c("duplications")]),
  Da1_up_singletons_length = length(Da1_up$nInteractions[Da1_up$combined_origin %in% c("singleton")]),
  Da2_up_duplicates_length = length(Da2_up$nInteractions[Da2_up$combined_origin %in% c("duplications")]),
  Da2_up_singletons_length = length(Da2_up$nInteractions[Da2_up$combined_origin %in% c("singleton")]),
  Da3_up_duplicates_length = length(Da3_up$nInteractions[Da3_up$combined_origin %in% c("duplications")]),
  Da3_up_singletons_length = length(Da3_up$nInteractions[Da3_up$combined_origin %in% c("singleton")]),
  Ga1_up_duplicates_length = length(Ga1_up$nInteractions[Ga1_up$combined_origin %in% c("duplications")]),
  Ga1_up_singletons_length = length(Ga1_up$nInteractions[Ga1_up$combined_origin %in% c("singleton")]),
  Ga2_up_duplicates_length = length(Ga2_up$nInteractions[Ga2_up$combined_origin %in% c("duplications")]),
  Ga2_up_singletons_length = length(Ga2_up$nInteractions[Ga2_up$combined_origin %in% c("singleton")]),
  Ga3_up_duplicates_length = length(Ga3_up$nInteractions[Ga3_up$combined_origin %in% c("duplications")]),
  Ga3_up_singletons_length = length(Ga3_up$nInteractions[Ga3_up$combined_origin %in% c("singleton")]),
  Da1_down_duplicates_length = length(Da1_down$nInteractions[Da1_down$combined_origin %in% c("duplications")]),
  Da1_down_singletons_length = length(Da1_down$nInteractions[Da1_down$combined_origin %in% c("singleton")]),
  Da2_down_duplicates_length = length(Da2_down$nInteractions[Da2_down$combined_origin %in% c("duplications")]),
  Da2_down_singletons_length = length(Da2_down$nInteractions[Da2_down$combined_origin %in% c("singleton")]),
  Da3_down_duplicates_length = length(Da3_down$nInteractions[Da3_down$combined_origin %in% c("duplications")]),
  Da3_down_singletons_length = length(Da3_down$nInteractions[Da3_down$combined_origin %in% c("singleton")]),
  Ga1_down_duplicates_length = length(Ga1_down$nInteractions[Ga1_down$combined_origin %in% c("duplications")]),
  Ga1_down_singletons_length = length(Ga1_down$nInteractions[Ga1_down$combined_origin %in% c("singleton")]),
  Ga2_down_duplicates_length = length(Ga2_down$nInteractions[Ga2_down$combined_origin %in% c("duplications")]),
  Ga2_down_singletons_length = length(Ga2_down$nInteractions[Ga2_down$combined_origin %in% c("singleton")]),
  Ga3_down_duplicates_length = length(Ga3_down$nInteractions[Ga3_down$combined_origin %in% c("duplications")]),
  Ga3_down_singletons_length = length(Ga3_down$nInteractions[Ga3_down$combined_origin %in% c("singleton")])
  
)

# Combine means and standard deviations into a data frame
stats_interactions <- data.frame(
  name = names,
  mean = means,
  sd = sds,
  size = size
)

rownames(stats_interactions) <- names

# Perform bootstrapping for each subset and compare means
n_iterations <- 10000
ci_results <- list()

for (i in 1:nrow(stats_interactions)) {
  subset_data_name <- stats_interactions$name[i]
  subset_mean <- stats_interactions$mean[i]
  subset_size <- stats_interactions$size[i]
  
  # Determine appropriate subset
  if (grepl("duplicates", subset_data_name)) {
    subset_data <- all_genes_interactions$nInteractions[all_genes_interactions$combined_origin %in% c("duplications")]
  } else if (grepl("singletons", subset_data_name)) {
    subset_data <- all_genes_interactions$nInteractions[all_genes_interactions$combined_origin %in% c("singleton")]
  } else {
    subset_data <- all_genes_interactions$nInteractions
  }
  
  # Bootstrapping
  bootstrapped_means <- numeric(n_iterations)
  for (j in 1:n_iterations) {
    resample <- sample(subset_data, size = subset_size, replace = FALSE)
    bootstrapped_means[j] <- mean(resample, na.rm = TRUE)
  }
  
  ci_lower <- quantile(bootstrapped_means, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrapped_means, 0.975, na.rm = TRUE)
  
  ci_results[[i]] <- list(
    name = subset_data_name,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    significant_low = subset_mean < ci_lower,
    significant_up = subset_mean > ci_upper
  )
}

# Convert list of CI results to data frame
ci_results_df <- do.call(rbind, lapply(ci_results, function(x) as.data.frame(t(unlist(x)))))

# Merge the results with the original statistics
stats_interactions <- merge(stats_interactions, ci_results_df, by = "name", all.x = TRUE)

datatable(stats_interactions, options = list(pageLength = 5)) %>%
  formatStyle(
    'significant_low.2.5%',
    backgroundColor = styleEqual(c("TRUE", "FALSE"), c('green', 'white')),
    color = styleEqual(c("TRUE", "FALSE"), c('white', 'black'))
  ) %>%
  formatStyle(
    'significant_up.97.5%',
    backgroundColor = styleEqual(c("TRUE", "FALSE"), c('green', 'white')),
    color = styleEqual(c("TRUE", "FALSE"), c('white', 'black'))
  )



# Protein-protein
# Calculate means and standard deviations for each group
# Define group names
names <- c("Da1_up",  "Da2_up",  "Da3_up",  "Ga1_up",  "Ga2_up",  "Ga3_up",
           "Da1_down",  "Da2_down",  "Da3_down",  "Ga1_down",  "Ga2_down",  "Ga3_down",
           "Da1_up_duplicates",  "Da1_up_singletons",  "Da2_up_duplicates",
           "Da2_up_singletons",  "Da3_up_duplicates",  "Da3_up_singletons",
           "Ga1_up_duplicates",  "Ga1_up_singletons",  "Ga2_up_duplicates",
           "Ga2_up_singletons",  "Ga3_up_duplicates",  "Ga3_up_singletons",
           "Da1_down_duplicates",  "Da1_down_singletons",  "Da2_down_duplicates",
           "Da2_down_singletons",  "Da3_down_duplicates",  "Da3_down_singletons",
           "Ga1_down_duplicates",  "Ga1_down_singletons",  "Ga2_down_duplicates",
           "Ga2_down_singletons",  "Ga3_down_duplicates",  "Ga3_down_singletons"
)
means <- c(
  Da1_up_mean = mean(Da1_up$p_nInteractions, na.rm = TRUE),
  Da2_up_mean = mean(Da2_up$p_nInteractions, na.rm = TRUE),
  Da3_up_mean = mean(Da3_up$p_nInteractions, na.rm = TRUE),
  Ga1_up_mean = mean(Ga1_up$p_nInteractions, na.rm = TRUE),
  Ga2_up_mean = mean(Ga2_up$p_nInteractions, na.rm = TRUE),
  Ga3_up_mean = mean(Ga3_up$p_nInteractions, na.rm = TRUE),
  Da1_down_mean = mean(Da1_down$p_nInteractions, na.rm = TRUE),
  Da2_down_mean = mean(Da2_down$p_nInteractions, na.rm = TRUE),
  Da3_down_mean = mean(Da3_down$p_nInteractions, na.rm = TRUE),
  Ga1_down_mean = mean(Ga1_down$p_nInteractions, na.rm = TRUE),
  Ga2_down_mean = mean(Ga2_down$p_nInteractions, na.rm = TRUE),
  Ga3_down_mean = mean(Ga3_down$p_nInteractions, na.rm = TRUE),
  Da1_up_duplicates_mean = mean(Da1_up$p_nInteractions[Da1_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da1_up_singletons_mean = mean(Da1_up$p_nInteractions[Da1_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da2_up_duplicates_mean = mean(Da2_up$p_nInteractions[Da2_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da2_up_singletons_mean = mean(Da2_up$p_nInteractions[Da2_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da3_up_duplicates_mean = mean(Da3_up$p_nInteractions[Da3_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da3_up_singletons_mean = mean(Da3_up$p_nInteractions[Da3_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga1_up_duplicates_mean = mean(Ga1_up$p_nInteractions[Ga1_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga1_up_singletons_mean = mean(Ga1_up$p_nInteractions[Ga1_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga2_up_duplicates_mean = mean(Ga2_up$p_nInteractions[Ga2_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga2_up_singletons_mean = mean(Ga2_up$p_nInteractions[Ga2_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga3_up_duplicates_mean = mean(Ga3_up$p_nInteractions[Ga3_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga3_up_singletons_mean = mean(Ga3_up$p_nInteractions[Ga3_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da1_down_duplicates_mean = mean(Da1_down$p_nInteractions[Da1_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da1_down_singletons_mean = mean(Da1_down$p_nInteractions[Da1_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da2_down_duplicates_mean = mean(Da2_down$p_nInteractions[Da2_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da2_down_singletons_mean = mean(Da2_down$p_nInteractions[Da2_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da3_down_duplicates_mean = mean(Da3_down$p_nInteractions[Da3_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da3_down_singletons_mean = mean(Da3_down$p_nInteractions[Da3_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga1_down_duplicates_mean = mean(Ga1_down$p_nInteractions[Ga1_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga1_down_singletons_mean = mean(Ga1_down$p_nInteractions[Ga1_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga2_down_duplicates_mean = mean(Ga2_down$p_nInteractions[Ga2_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga2_down_singletons_mean = mean(Ga2_down$p_nInteractions[Ga2_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga3_down_duplicates_mean = mean(Ga3_down$p_nInteractions[Ga3_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga3_down_singletons_mean = mean(Ga3_down$p_nInteractions[Ga3_down$combined_origin %in% c("singleton")], na.rm = TRUE)
)

sds <- c(
  Da1_up_sd = sd(Da1_up$p_nInteractions, na.rm = TRUE),
  Da2_up_sd = sd(Da2_up$p_nInteractions, na.rm = TRUE),
  Da3_up_sd = sd(Da3_up$p_nInteractions, na.rm = TRUE),
  Ga1_up_sd = sd(Ga1_up$p_nInteractions, na.rm = TRUE),
  Ga2_up_sd = sd(Ga2_up$p_nInteractions, na.rm = TRUE),
  Ga3_up_sd = sd(Ga3_up$p_nInteractions, na.rm = TRUE),
  Da1_down_sd = sd(Da1_down$p_nInteractions, na.rm = TRUE),
  Da2_down_sd = sd(Da2_down$p_nInteractions, na.rm = TRUE),
  Da3_down_sd = sd(Da3_down$p_nInteractions, na.rm = TRUE),
  Ga1_down_sd = sd(Ga1_down$p_nInteractions, na.rm = TRUE),
  Ga2_down_sd = sd(Ga2_down$p_nInteractions, na.rm = TRUE),
  Ga3_down_sd = sd(Ga3_down$p_nInteractions, na.rm = TRUE),
  Da1_up_duplicates_sd = sd(Da1_up$p_nInteractions[Da1_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da1_up_singletons_sd = sd(Da1_up$p_nInteractions[Da1_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da2_up_duplicates_sd = sd(Da2_up$p_nInteractions[Da2_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da2_up_singletons_sd = sd(Da2_up$p_nInteractions[Da2_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da3_up_duplicates_sd = sd(Da3_up$p_nInteractions[Da3_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da3_up_singletons_sd = sd(Da3_up$p_nInteractions[Da3_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga1_up_duplicates_sd = sd(Ga1_up$p_nInteractions[Ga1_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga1_up_singletons_sd = sd(Ga1_up$p_nInteractions[Ga1_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga2_up_duplicates_sd = sd(Ga2_up$p_nInteractions[Ga2_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga2_up_singletons_sd = sd(Ga2_up$p_nInteractions[Ga2_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga3_up_duplicates_sd = sd(Ga3_up$p_nInteractions[Ga3_up$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga3_up_singletons_sd = sd(Ga3_up$p_nInteractions[Ga3_up$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da1_down_duplicates_sd = sd(Da1_down$p_nInteractions[Da1_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da1_down_singletons_sd = sd(Da1_down$p_nInteractions[Da1_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da2_down_duplicates_sd = sd(Da2_down$p_nInteractions[Da2_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da2_down_singletons_sd = sd(Da2_down$p_nInteractions[Da2_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Da3_down_duplicates_sd = sd(Da3_down$p_nInteractions[Da3_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Da3_down_singletons_sd = sd(Da3_down$p_nInteractions[Da3_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga1_down_duplicates_sd = sd(Ga1_down$p_nInteractions[Ga1_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga1_down_singletons_sd = sd(Ga1_down$p_nInteractions[Ga1_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga2_down_duplicates_sd = sd(Ga2_down$p_nInteractions[Ga2_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga2_down_singletons_sd = sd(Ga2_down$p_nInteractions[Ga2_down$combined_origin %in% c("singleton")], na.rm = TRUE),
  Ga3_down_duplicates_sd = sd(Ga3_down$p_nInteractions[Ga3_down$combined_origin %in% c("duplications")], na.rm = TRUE),
  Ga3_down_singletons_sd = sd(Ga3_down$p_nInteractions[Ga3_down$combined_origin %in% c("singleton")], na.rm = TRUE)
)

size <- c(
  Da1_up_length = length(Da1_up$p_nInteractions),
  Da2_up_length = length(Da2_up$p_nInteractions),
  Da3_up_length = length(Da3_up$p_nInteractions),
  Ga1_up_length = length(Ga1_up$p_nInteractions),
  Ga2_up_length = length(Ga2_up$p_nInteractions),
  Ga3_up_length = length(Ga3_up$p_nInteractions),
  Da1_down_length = length(Da1_down$p_nInteractions),
  Da2_down_length = length(Da2_down$p_nInteractions),
  Da3_down_length = length(Da3_down$p_nInteractions),
  Ga1_down_length = length(Ga1_down$p_nInteractions),
  Ga2_down_length = length(Ga2_down$p_nInteractions),
  Ga3_down_length = length(Ga3_down$p_nInteractions),
  Da1_up_duplicates_length = length(Da1_up$p_nInteractions[Da1_up$combined_origin %in% c("duplications")]),
  Da1_up_singletons_length = length(Da1_up$p_nInteractions[Da1_up$combined_origin %in% c("singleton")]),
  Da2_up_duplicates_length = length(Da2_up$p_nInteractions[Da2_up$combined_origin %in% c("duplications")]),
  Da2_up_singletons_length = length(Da2_up$p_nInteractions[Da2_up$combined_origin %in% c("singleton")]),
  Da3_up_duplicates_length = length(Da3_up$p_nInteractions[Da3_up$combined_origin %in% c("duplications")]),
  Da3_up_singletons_length = length(Da3_up$p_nInteractions[Da3_up$combined_origin %in% c("singleton")]),
  Ga1_up_duplicates_length = length(Ga1_up$p_nInteractions[Ga1_up$combined_origin %in% c("duplications")]),
  Ga1_up_singletons_length = length(Ga1_up$p_nInteractions[Ga1_up$combined_origin %in% c("singleton")]),
  Ga2_up_duplicates_length = length(Ga2_up$p_nInteractions[Ga2_up$combined_origin %in% c("duplications")]),
  Ga2_up_singletons_length = length(Ga2_up$p_nInteractions[Ga2_up$combined_origin %in% c("singleton")]),
  Ga3_up_duplicates_length = length(Ga3_up$p_nInteractions[Ga3_up$combined_origin %in% c("duplications")]),
  Ga3_up_singletons_length = length(Ga3_up$p_nInteractions[Ga3_up$combined_origin %in% c("singleton")]),
  Da1_down_duplicates_length = length(Da1_down$p_nInteractions[Da1_down$combined_origin %in% c("duplications")]),
  Da1_down_singletons_length = length(Da1_down$p_nInteractions[Da1_down$combined_origin %in% c("singleton")]),
  Da2_down_duplicates_length = length(Da2_down$p_nInteractions[Da2_down$combined_origin %in% c("duplications")]),
  Da2_down_singletons_length = length(Da2_down$p_nInteractions[Da2_down$combined_origin %in% c("singleton")]),
  Da3_down_duplicates_length = length(Da3_down$p_nInteractions[Da3_down$combined_origin %in% c("duplications")]),
  Da3_down_singletons_length = length(Da3_down$p_nInteractions[Da3_down$combined_origin %in% c("singleton")]),
  Ga1_down_duplicates_length = length(Ga1_down$p_nInteractions[Ga1_down$combined_origin %in% c("duplications")]),
  Ga1_down_singletons_length = length(Ga1_down$p_nInteractions[Ga1_down$combined_origin %in% c("singleton")]),
  Ga2_down_duplicates_length = length(Ga2_down$p_nInteractions[Ga2_down$combined_origin %in% c("duplications")]),
  Ga2_down_singletons_length = length(Ga2_down$p_nInteractions[Ga2_down$combined_origin %in% c("singleton")]),
  Ga3_down_duplicates_length = length(Ga3_down$p_nInteractions[Ga3_down$combined_origin %in% c("duplications")]),
  Ga3_down_singletons_length = length(Ga3_down$p_nInteractions[Ga3_down$combined_origin %in% c("singleton")])
  
)

# Combine means and standard deviations into a data frame
stats_interactions_p <- data.frame(
  name = names,
  mean = means,
  sd = sds,
  size = size
)

rownames(stats_interactions_p) <- names

# Perform bootstrapping for each subset and compare means
n_iterations <- 10000
ci_results <- list()

for (i in 1:nrow(stats_interactions_p)) {
  subset_data_name <- stats_interactions_p$name[i]
  subset_mean <- stats_interactions_p$mean[i]
  subset_size <- stats_interactions_p$size[i]
  
  # Determine appropriate subset
  if (grepl("duplicates", subset_data_name)) {
    subset_data <- all_genes_interactions$p_nInteractions[all_genes_interactions$combined_origin %in% c("duplications")]
  } else if (grepl("singletons", subset_data_name)) {
    subset_data <- all_genes_interactions$p_nInteractions[all_genes_interactions$combined_origin %in% c("singleton")]
  } else {
    subset_data <- all_genes_interactions$p_nInteractions
  }
  
  # Bootstrapping
  bootstrapped_means <- numeric(n_iterations)
  for (j in 1:n_iterations) {
    resample <- sample(subset_data, size = subset_size, replace = FALSE)
    bootstrapped_means[j] <- mean(resample, na.rm = TRUE)
  }
  
  ci_lower <- quantile(bootstrapped_means, 0.025, na.rm = TRUE)
  ci_upper <- quantile(bootstrapped_means, 0.975, na.rm = TRUE)
  
  ci_results[[i]] <- list(
    name = subset_data_name,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    significant_low = subset_mean < ci_lower,
    significant_up = subset_mean > ci_upper
  )
}

# Convert list of CI results to data frame
ci_results_df <- do.call(rbind, lapply(ci_results, function(x) as.data.frame(t(unlist(x)))))

# Merge the results with the original statistics
stats_interactions_p <- merge(stats_interactions_p, ci_results_df, by = "name", all.x = TRUE)

datatable(stats_interactions_p, options = list(pageLength = 5)) %>%
  formatStyle(
    'significant_low.2.5%',
    backgroundColor = styleEqual(c("TRUE", "FALSE"), c('green', 'white')),
    color = styleEqual(c("TRUE", "FALSE"), c('white', 'black'))
  ) %>%
  formatStyle(
    'significant_up.97.5%',
    backgroundColor = styleEqual(c("TRUE", "FALSE"), c('green', 'white')),
    color = styleEqual(c("TRUE", "FALSE"), c('white', 'black'))
  )


# Study of the glycerol pathway

# Prepare gene data vectors for each line
Da1_gene_data <- setNames(Da1_significant$logFC, Da1_significant$gene_id)
Da2_gene_data <- setNames(Da2_significant$logFC, Da2_significant$gene_id)
Da3_gene_data <- setNames(Da3_significant$logFC, Da3_significant$gene_id)
Ga1_gene_data <- setNames(Ga1_significant$logFC, Ga1_significant$gene_id)
Ga2_gene_data <- setNames(Ga2_significant$logFC, Ga2_significant$gene_id)
Ga3_gene_data <- setNames(Ga3_significant$logFC, Ga3_significant$gene_id)

# Convert named vectors to data frames
Da1_df <- data.frame(gene_id = names(Da1_gene_data), Da1 = Da1_gene_data, stringsAsFactors = FALSE)
Da2_df <- data.frame(gene_id = names(Da2_gene_data), Da2 = Da2_gene_data, stringsAsFactors = FALSE)
Da3_df <- data.frame(gene_id = names(Da3_gene_data), Da3 = Da3_gene_data, stringsAsFactors = FALSE)
Ga1_df <- data.frame(gene_id = names(Ga1_gene_data), Ga1 = Ga1_gene_data, stringsAsFactors = FALSE)
Ga2_df <- data.frame(gene_id = names(Ga2_gene_data), Ga2 = Ga2_gene_data, stringsAsFactors = FALSE)
Ga3_df <- data.frame(gene_id = names(Ga3_gene_data), Ga3 = Ga3_gene_data, stringsAsFactors = FALSE)

# Merge all data frames by 'gene_id' using Reduce and merge
merged_data_pathview <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
                      list(Da1_df, Da2_df, Da3_df, Ga1_df, Ga2_df, Ga3_df))

# Set the 'gene_id' column as row names
rownames(merged_data_pathview) <- merged_data_pathview$gene_id

# Remove the 'gene_id' column from the data frame
merged_data_pathview <- merged_data_pathview[ ,-1]

pathway_id <- "sce00561"
# Generate pathview for the glycerolipid metabolism
pathview(gene.data = merged_data_pathview, pathway.id = pathway_id, gene.idtype = "orf", species = "sce", limit = c(-2, 2), out.suffix = "Combined")

pathway_id <- "sce00564"
# Generate pathview for the glycerophospholipid metabolism
pathview(gene.data = merged_data_pathview, pathway.id = pathway_id, gene.idtype = "orf", species = "sce", limit = c(-2, 2), out.suffix = "Combined")

#pathway sce00561
all_lines_kegg2 <- compareCluster(all_genes, fun = "enrichKEGG", organism = "sce", pvalueCutoff = 1)
all_lines_kegg2@compareClusterResult$Description <- sapply(strsplit(all_lines_kegg2@compareClusterResult$Description, " - S"), `[`, 1)

kegg_results_df <- as.data.frame(all_lines_kegg2)
sce00561_results <- subset(kegg_results_df, ID == "sce00561")

glycerol_gene_ids <- sce00561_results$geneID
glycerol_gene_ids <- unlist(strsplit("YBL011W/YCR068W/YCR105W/YDL052C/YDR058C/YDR284C/YDR503C/YER037W/YER062C/YER073W/YFL053W/YGR205W/YHL032C/YIL053W/YKL094W/YKR067W/YKR089C/YML070W/YMR110C/YMR165C/YMR169C/YMR170C/YMR313C/YMR318C/YNR008W/YOR081C/YOR120W/YOR175C/YOR245C/YOR374W/YPL061W/YPR139C", "/"))

# Initialize a list of DEG sets with logFC values
deg_sets <- list(
  Da1 = Da1_significant[, c("gene_id", "logFC")],
  Da2 = Da2_significant[, c("gene_id", "logFC")],
  Da3 = Da3_significant[, c("gene_id", "logFC")],
  Ga1 = Ga1_significant[, c("gene_id", "logFC")],
  Ga2 = Ga2_significant[, c("gene_id", "logFC")],
  Ga3 = Ga3_significant[, c("gene_id", "logFC")]
)

# Create a matrix to store logFC values
logFC_matrix <- sapply(deg_sets, function(set) {
  # Match logFC values for glycerol genes
  logFC_values <- set$logFC[match(glycerol_gene_ids, set$gene_id)]
  return(logFC_values)
})
logFC_matrix[is.na(logFC_matrix)] <- 0

# Set row and column names for the matrix
rownames(logFC_matrix) <- glycerol_gene_ids
colnames(logFC_matrix) <- names(deg_sets)

# Ensure rownames are gene IDs for gene_origin
gene_origin_subset <- gene_origin[glycerol_gene_ids, , drop = FALSE]  # Subset to only include glycerolipid genes

# Create the annotation data frame
origin_info <- data.frame(combined_origin = gene_origin_subset$combined_origin)

# Set row names to match the glycerol gene IDs (they should already match if subset was successful)
rownames(origin_info) <- rownames(gene_origin_subset)
colnames(origin_info)[1] <- "Origin"


# Get the enzyme activity codes for each gene

EC_glycerol <- c("2.3.1.15", "3.1.1.3", "1.1.1.2", "2.3.1.51", "3.1.1.3", "3.1.3.4", "3.1.3.4", "3.1.3.106",
                 "3.1.3.21", "1.2.1.3", "2.7.1.29", "2.7.1.31", "2.7.1.30", "3.1.3.21", "3.1.1.23", "2.3.1.15",
                 "3.1.1.3", "2.7.1.29", "1.2.1.3", "3.1.3.4", "1.2.1.3", "1.2.1.3", "3.1.1.3", "1.1.1.2", 
                 "2.3.1.158", "3.1.1.3", "1.1.1.156", "2.3.1.51", "2.3.1.22", "1.2.1.3", "1.2.1.3", 
                 "2.3.1.51")

EC_glycerol_df <- data.frame(Enzyme = EC_glycerol, row.names = glycerol_gene_ids)

# Combine the two annotations into a single data frame
combined_annotation <- cbind(origin_info, EC_glycerol_df)
combined_annotation_sorted <- combined_annotation[order(combined_annotation$Enzyme, na.last = TRUE),, drop =FALSE]

sorted_rownames <- rownames(combined_annotation_sorted)
logFC_matrix_sorted <- logFC_matrix[sorted_rownames, , drop = FALSE]

# Create the heatmap
pheatmap(logFC_matrix_sorted,
         cluster_rows = FALSE,
         cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color gradient: blue = downregulated, red = upregulated
         annotation_row = combined_annotation_sorted,
         na_col = "grey",  # Color for NA values
         main = "logFC of Glycerolipid Metabolism Genes in DEGs Sets")

#glycerol fisher test

contingency_glycerol <- matrix(c(15, 3, 8, 6), nrow = 2, byrow =TRUE)
rownames(contingency_glycerol) <- c("Duplicated", "Singleton")
colnames(contingency_glycerol) <- c("DEG", "No-DEG")

test_glycerol <- fisher.test(contingency_glycerol)

