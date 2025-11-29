# DEPM Project by Domenico Azzarito, Federico Lattanzio and Michele Pezza


# Project Steps

# 1. Selection of the disease: Kidney renal papillary cell carcinoma (KIRP)

# 2. Differentially expressed genes (DEGs)

# 3. CoExpression Networks

# 4. Differential CoExpression Network

# 5. Patient Similarity Network


##########################################################################



# TCGA-KIRP

rm(list = ls())      # clear workspace
gc()                 # garbage collection


# The first step is to loading the required packages for the analysis


library(BiocGenerics) # Contains generic Bioconductor functions. # Used internally by SummarizedExperiment and DESeq2

library(DESeq2) # Store RNA-Seq counts and normalize counts

library(psych) # for correlation analysis

library(NetworkToolbox) # Additional functions for network analysis

library(ggplot2) # volcano plot, degree distribution plot and any additional plot

# library(ggnet)    

library(GGally) # Package containing ggnet2(), essential function for plotting

library(sna) # functions for clustering coefficient, components..

library(network) # defines the network() object type

library(TCGAbiolinks) # queries the GDC portal and downloads TCGA data

library(GenomicRanges)  # Gene annotation is stored as GRanges

library(SummarizedExperiment)  # TCGA RNA-Seq is stored in a SummarizedExperiment object

library(DT)  # Interactive visualization

sessionInfo()


# We proceed to download the data

############### Download tumor RNA-Seq data (Primary Tumor)

proj <- "TCGA-KIRP"
dir.create(file.path(proj))

# RNA-Seq data: primary tumor samples
rna.query.C <- TCGAbiolinks::GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# Download the data from GDC to local cache folder "GDCdata" (run only once)
#GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")

# Convert the downloaded files into a SummarizedExperiment object
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata")

# Extract the raw count matrix (genes x tumor samples)
rna.expr.data.C <- assay(rna.data.C)

# Extract gene annotation (Ensembl IDs, gene symbols, chromosomes, etc.)
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))

# Checks

head(rna.query.C$results[[1]]$cases)

rna.data.C

dim(rna.expr.data.C)

head(genes.info)

############### Download normal samples

# RNA-Seq data: solid tissue normal samples
rna.query.N <- TCGAbiolinks::GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)

# Download normal data
#GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")

# Prepare SummarizedExperiment for normal samples
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata")

# Extract raw counts for normal samples
rna.expr.data.N <- assay(rna.data.N)

# Extract gene annotation for normal dataset
genes.info2 <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))

# Sanity check: tumor and normal gene annotations should match
all(na.omit(genes.info2) == na.omit(genes.info))

# Checks

dim(rna.expr.data.N)
rna.data.N

# Results:
# Tumor: 60,660 genes × 290 samples
# Normal: 60,660 genes × 32 samples
# Gene annotation: identical between tumor and normal

#########################################################################

# Clinical Data (To be done, not needed now)


#######################################################################


#########################################################################
# Data cleaning 
#########################################################################

## Inspect basic structure ------------------------------------------

# Dimensions: genes x samples
dim(rna.expr.data.C)  # tumor
dim(rna.expr.data.N)  # normal

# Quick look at column names (TCGA barcodes)
head(colnames(rna.expr.data.C))
head(colnames(rna.expr.data.N))

# Extract patient IDs (first 12 characters of barcode)
patients.C <- substr(colnames(rna.expr.data.C), 1, 12)
patients.N <- substr(colnames(rna.expr.data.N), 1, 12)

# How many unique patients in each group?
length(unique(patients.C))
length(unique(patients.N))

# Distribution of tumor samples per patient (to detect duplicates)
sort(table(patients.C))


## Keep only patients with exactly ONE tumor sample ------------------

# Table: how many tumor samples per patient
tab.C <- table(patients.C)

# Select patient IDs with exactly 1 primary tumor sample
single.C <- names(tab.C[tab.C == 1])

# Indices of columns in tumor matrix corresponding to single-sample patients
idx.single.C <- which(patients.C %in% single.C)

# Tumor expression matrix: only one tumor sample per patient
expr.C <- as.data.frame(rna.expr.data.C[, idx.single.C])

# Normal expression matrix: start with all normals
expr.N <- as.data.frame(rna.expr.data.N)


## Rename columns to pure patient IDs --------------------------------

# Use only the 12-character patient ID as column names
colnames(expr.C) <- substr(colnames(expr.C), 1, 12)
colnames(expr.N) <- substr(colnames(expr.N), 1, 12)

# quick check
head(colnames(expr.C))
head(colnames(expr.N))

# Make sure there are no duplicated column names now
sum(duplicated(colnames(expr.C)))
sum(duplicated(colnames(expr.N)))


## Keep only patients with BOTH tumor and normal samples -------------

# Patient IDs present in BOTH expr.C (tumor) and expr.N (normal)
common_ids <- intersect(colnames(expr.C), colnames(expr.N))
length(common_ids)   # this is the number of paired patients

# Subset both matrices to the same set of patients
expr.C <- expr.C[, common_ids, drop = FALSE]
expr.N <- expr.N[, common_ids, drop = FALSE]

# Sanity checks
ncol(expr.C)                      # number of paired tumor samples
ncol(expr.N)                      # number of paired normal samples
all(colnames(expr.C) == colnames(expr.N))   # should be TRUE


## Final technical checks on counts ----------------------------------

# Check that genes (rows) are aligned between tumor and normal
all(rownames(expr.C) == rownames(expr.N))   # must be TRUE

# Check data types and missingness
typeof(expr.C[1, 1])                # should be "integer" or "double"
any(is.na(expr.C))                  # should be FALSE
any(is.nan(as.matrix(expr.C)))      # should be FALSE

typeof(expr.N[1, 1])
any(is.na(expr.N))
any(is.nan(as.matrix(expr.N)))

#########################################################################
# At this point:
# - expr.C = raw tumor counts, genes x N_paired patients
# - expr.N = raw normal counts, genes x N_paired patients
# - Columns in expr.C and expr.N are the SAME patients, SAME order
# These are the matrices we will feed into DESeq2 for normalization.
#########################################################################


#########################################################################
# Normalizing data with DESeq2 
#########################################################################

## Check alignment of tumor and normal matrices ----------------------

# Must be TRUE (same gene order)
all(rownames(expr.C) == rownames(expr.N))

# Combine into full dataset: normals first, then tumors
full.data <- cbind(expr.N, expr.C)
full.data <- data.frame(full.data)

dim(full.data)   # should be 60660 x 64 (32 normal + 32 tumor)


## Build metadata -----------------------------------------------------

metad <- data.frame(
  condition = factor(c(
    rep("normal", ncol(expr.N)),
    rep("cancer", ncol(expr.C))
  ))
)

# Row names must match full.data column names
rownames(metad) <- colnames(full.data)

# DESeq2 expects gene_id as a column when tidy=TRUE
full.data <- cbind(gene_id = rownames(full.data), full.data)


## Build DESeq2 dataset ----------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = full.data,
  colData   = metad,
  design    = ~ condition,
  tidy      = TRUE
)

# Sanity check
dim(counts(dds))   # should be 60660 x 64


## Filter low-expression genes ---------------------------------------

# At least 10 counts in 90% of patients for each group
# We have 32 normal and 32 cancer → threshold = 0.90 * 32 ≈ 28.8 → floor = 28

threshold <- floor(0.9 * ncol(expr.N))   # = 28

keep <- rowSums(counts(dds) >= 10) >= threshold
dds <- dds[keep, ]

dim(counts(dds))  # reduced number of genes 


##  Normalize using DESeq2 size factors -------------------------------

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)

# No gene should be zero in all samples
sum(rowSums(normalized_counts == 0) == ncol(normalized_counts))


## Split normalized matrix back into normal/tumor ---------------------

filtr.expr.n <- as.data.frame(
  normalized_counts[, 1:ncol(expr.N)]
)

filtr.expr.c <- as.data.frame(
  normalized_counts[, (ncol(expr.N) + 1):ncol(normalized_counts)]
)

# Keep original patient IDs (DESeq2 sometimes adds .1, .2)
colnames(filtr.expr.n) <- colnames(expr.N)
colnames(filtr.expr.c) <- colnames(expr.C)


## Final checks ------------------------------------------------------

# Same genes
all(rownames(filtr.expr.c) == rownames(filtr.expr.n))

# Same patients
all(colnames(filtr.expr.c) == colnames(filtr.expr.n))

nrow(filtr.expr.c)   # number of filtered/normalized genes
ncol(filtr.expr.c)   # should be 32

sizeFactors(dds)
summary(sizeFactors(dds))


#########################################################################
# Output:
#   filtr.expr.n = normalized expression for 32 normal samples
#   filtr.expr.c = normalized expression for 32 tumor samples
#
# These matrices are what we will use for:
#   • DEGs
#   • Volcano Plot
#   • Co-expression networks
#   • Differential networks
#   • PSN and SNF
#########################################################################

##############################################################################

##############################################################################
# Differential Expression Analysis (DEGs) — UPDATED THRESHOLDS
##############################################################################

## Compute log2 Fold Change --------------------------------------------

log2fc <- log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n))
names(log2fc) <- rownames(filtr.expr.c)

## Compute paired t-test p-values --------------------------------------

pval <- sapply(
  1:nrow(filtr.expr.c),
  function(i) t.test(
    as.numeric(filtr.expr.c[i, ]),
    as.numeric(filtr.expr.n[i, ]),
    paired = TRUE
  )$p.value
)

## FDR correction -------------------------------------------------------

pval_fdr <- p.adjust(pval, method = "fdr")

## Build DEG table ------------------------------------------------------

deg.table <- data.frame(
  gene        = rownames(filtr.expr.c),
  log2FC      = log2fc,
  pvalue      = pval,
  pvalue_fdr  = pval_fdr,
  row.names   = rownames(filtr.expr.c)
)

## Apply NEW thresholds -------------------------------------------------

FDR_threshold <- 0.001
FC_threshold  <- 2.0

deg.genes <- subset(
  deg.table,
  abs(log2FC) >= FC_threshold & pvalue_fdr <= FDR_threshold
)

cat("Number of DEGs:", nrow(deg.genes), "\n")

## Volcano plot ---------------------------------------------------------

deg.table$color <- "NO"
deg.table$color[deg.table$log2FC >= FC_threshold & deg.table$pvalue_fdr <= FDR_threshold] <- "UP"
deg.table$color[deg.table$log2FC <= -FC_threshold & deg.table$pvalue_fdr <= FDR_threshold] <- "DOWN"
deg.table$color <- factor(deg.table$color, levels = c("DOWN", "NO", "UP"))

ggplot(deg.table, aes(x = log2FC, y = -log10(pvalue_fdr), color = color)) +
  geom_point(alpha = 0.6) +
  xlab("log2 Fold Change") +
  ylab("-log10 FDR") +
  geom_vline(xintercept = c(-FC_threshold, FC_threshold), col = "red") +
  geom_hline(yintercept = -log10(FDR_threshold), col = "red") +
  theme_minimal()


DEG_list <- rownames(deg.genes)
length(DEG_list)


#########################################################################

# CoExpression Networks

########################################################################

## Extract the DEG expression matrices -----------------------------------

# Keep only DEGs in the normalized matrices
expr_deg_cancer <- filtr.expr.c[DEG_list, ]
expr_deg_normal <- filtr.expr.n[DEG_list, ]

# Log transform (as asked by the professor)
expr_deg_cancer <- log2(expr_deg_cancer + 1)
expr_deg_normal <- log2(expr_deg_normal + 1)

# Check dimensions
dim(expr_deg_cancer)
dim(expr_deg_normal)

# The number of columns must be 32 in both
# The number of rows must match number of DEGs
all(rownames(expr_deg_cancer) == rownames(expr_deg_normal))


## Compute adjacency matrices for cancer and normal -----------------------

# Set correlation threshold 
corr_threshold <- 0.7

## Compute adjacency matrices for cancer and normal -----------------------

# Cancer network
# cor is faster then test.corr if you have only to compute the rho as here
cor_mat_c <- cor(t(expr_deg_cancer), method = "pearson")  # Pearson correlation
diag(cor_mat_c) <- 0  # Remove self-correlations
adj_mat_c <- (abs(cor_mat_c) >= corr_threshold) * 1 # Binary adjacency

# Normal network
cor_mat_n <- cor(t(expr_deg_normal), method = "pearson")
diag(cor_mat_n) <- 0
adj_mat_n <- (abs(cor_mat_n) >= corr_threshold) * 1

## Compute degree for each network ---------------------------------------

degree_c <- colSums(adj_mat_c != 0)
names(degree_c) <- rownames(adj_mat_c)
degree_c <- sort(degree_c, decreasing = TRUE)

degree_n <- colSums(adj_mat_n != 0)
names(degree_n) <- rownames(adj_mat_n)
degree_n <- sort(degree_n, decreasing = TRUE)

## Check if scale-free network (plot degree distribution) ----------------

# Cancer degree distribution
hist(degree_c, breaks = 50, main = "Cancer Network Degree Distribution",
     xlab = "Degree", col = "lightblue", freq = FALSE)

# To check scale-free: plot log-log (degree vs frequency)
deg_freq_c <- table(degree_c)
plot(log(as.numeric(names(deg_freq_c))), log(deg_freq_c), 
     main = "Cancer: Log-Log Degree Distribution", xlab = "log(Degree)", ylab = "log(Frequency)", 
     pch = 16, col = "blue")

# Normal degree distribution
hist(degree_n, breaks = 50, main = "Normal Network Degree Distribution",
     xlab = "Degree", col = "lightgreen", freq = FALSE)

deg_freq_n <- table(degree_n)
plot(log(as.numeric(names(deg_freq_n))), log(deg_freq_n), 
     main = "Normal: Log-Log Degree Distribution", xlab = "log(Degree)", ylab = "log(Frequency)", 
     pch = 16, col = "green")

## Find hubs (top 5% highest degree) -------------------------------------

y <- quantile(degree_n[degree_n>0], 0.95)
x <- quantile(degree_c[degree_c>0], 0.95)

hubs_c <- degree_c[degree_c>=x]

hubs_n <- degree_n[degree_n>=y]


## Compare hubs between conditions ---------------------------------------

# Common hubs
common_hubs <- intersect(names(hubs_c), names(hubs_n))

# Hubs unique to cancer
selective_hubs_c <- setdiff(names(hubs_c), names(hubs_n))

# Hubs unique to normal
selective_hubs_n <- setdiff(names(hubs_n), names(hubs_c))

# Output
cat("Number of cancer hubs:", length(hubs_c), "\n")
cat("Number of normal hubs:", length(hubs_n), "\n")
cat("Number of common hubs:", length(common_hubs), "\n")
cat("Number of selective hubs in cancer:", length(selective_hubs_c), "\n")
cat("Number of selective hubs in normal:", length(selective_hubs_n), "\n")

#########################################################################
# Output:
#   adj_mat_c = adjacency matrix for cancer co-expression network
#   adj_mat_n = adjacency matrix for normal co-expression network
#   hubs_c = hubs in cancer network
#   hubs_n = hubs in normal network
#########################################################################

#########################################################################
#
# Differential Co-expressed Network
#
#########################################################################

# Set Z-score threshold 
z_threshold <- 4 # SCELTO PIU' RESTRITTIVO PER FAR VENIRE RETE MENO DENSA -> POWER-LAW +-

# Sample sizes
n_c <- ncol(expr_deg_cancer)  # Number of cancer samples
n_n <- ncol(expr_deg_normal)  # Number of normal samples

## Compute differential adjacency matrix ----------------------------------

# Function to compute Fisher's Z-transform
fisher_z <- function(r) {
  0.5 * log((1 + r) / (1 - r))
}

# Compute Z-transforms for cancer and normal correlations
z_c <- fisher_z(cor_mat_c)
z_n <- fisher_z(cor_mat_n)

# Compute Z-score for difference: (z_c - z_n) / sqrt(1/(n_c-3) + 1/(n_n-3))
se <- sqrt(1/(n_c - 3) + 1/(n_n - 3))
z_diff <- (z_c - z_n) / se

# Binary adjacency
adj_mat_diff <- z_diff
adj_mat_diff <- adj_mat_diff_bin <- (abs(z_diff) >= z_threshold) * 1

# Remove diagonal
diag(adj_mat_diff) <- 0

## Compute degree for differential network -------------------------------

degree_diff <- colSums(adj_mat_diff != 0)
names(degree_diff) <- rownames(adj_mat_diff)
degree_diff <- sort(degree_diff, decreasing = TRUE)

## Check if scale-free network (plot degree distribution) ----------------

# Differential network degree distribution
hist(degree_diff, breaks = 50, main = "Differential Network Degree Distribution", 
     xlab = "Degree", col = "lightcoral", freq = FALSE)

# Log-log plot
deg_freq_diff <- table(degree_diff)
if (length(deg_freq_diff) > 1) {
  plot(log(as.numeric(names(deg_freq_diff))), log(deg_freq_diff), 
       main = "Differential: Log-Log Degree Distribution", xlab = "log(Degree)", ylab = "log(Frequency)", 
       pch = 16, col = "red")
}

## Find hubs (top 5% highest degree) -------------------------------------

z <- quantile(degree_diff[degree_diff>0], 0.95)
hubs_diff <- degree_diff[degree_diff>=z]


## Compare hubs with task 3 (cancer and normal hubs) ---------------------

# Common with cancer hubs
common_hubs_diff_c <- intersect(names(hubs_diff), names(hubs_c))

# Common with normal hubs
common_hubs_diff_n <- intersect(names(hubs_diff), names(hubs_n))

# Unique to differential
selective_hubs_diff <- setdiff(names(hubs_diff), union(names(hubs_c), names(hubs_n)))

# Output
cat("Number of differential hubs:", length(hubs_diff), "\n")
cat("Number of hubs common with cancer:", length(common_hubs_diff_c), "\n")
cat("Number of hubs common with normal:", length(common_hubs_diff_n), "\n")
cat("Number of selective hubs in differential:", length(selective_hubs_diff), "\n")

## Subnetwork Plot of the Most Relevant Gene (Highest Degree in Differential Network) --------------------

# Find the gene with the highest degree in the differential network
top_gene <- names(degree_diff)[1]  # Highest degree
top_gene_degree <- degree_diff[1]

cat("Most relevant gene:", top_gene, "with degree:", top_gene_degree, "\n")

# Find neighbors: genes connected to top_gene in adj_mat_diff
neighbors <- which(adj_mat_diff[top_gene, ] != 0)
neighbor_names <- colnames(adj_mat_diff)[neighbors]

# Subset: top_gene and its neighbors
subnet_nodes <- unique(c(top_gene, neighbor_names))
subnet_adj <- adj_mat_diff[subnet_nodes, subnet_nodes]

# Create network object
net_subnet <- network(subnet_adj, matrix.type = "adjacency", ignore.eval = FALSE, names.eval = "weights", directed = FALSE)

# Set attributes for plotting
net_subnet %v% "type" <- ifelse(network.vertex.names(net_subnet) == top_gene, "top_gene", "neighbor")
net_subnet %v% "color" <- ifelse(net_subnet %v% "type" == "top_gene", "red", "blue")
network::set.edge.attribute(net_subnet, "edgecolor", "green") 

# Plot subnetwork
ggnet2(net_subnet, color = "color", alpha = 0.7, size = 3,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.5,
       node.label = NULL, label.color = "black", label.size = 4) +
  ggtitle(paste("Subnetwork of", top_gene, "(Highest Degree in Differential Network)")) +
  guides(size = "none")


#########################################################################
# Output:
#   adj_mat_diff = adjacency matrix for differential co-expression network
#   hubs_diff = hubs in differential network
#########################################################################




