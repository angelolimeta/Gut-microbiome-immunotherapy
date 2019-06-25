# ==== IMPORT DATA ====
library("phyloseq")
# Import Phyloseq object for all studies
physeq = readRDS(file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/physeq.rds")

# ==== SETUP DESIGN MATRIX ====

# Let's take a look at the sample data
head(sample_data(physeq), 25)
# Check factors
lapply(sample_data(physeq), class)

# Load DESeq2
#BiocManager::install("DESeq2")
library("DESeq2")

# ==== RUN DESeq2 FOR ALL SAMPLES ====

# Convert phyloseq object to DEseq object
ds_full = phyloseq_to_deseq2(physeq, ~ Study + Treatment + Response)
# OTU matrix has many zeros. We need to re-estimate the size factors.
# Iterative normalisation does not converge, lets use poscounts instead.
ds_full = estimateSizeFactors(ds_full, type = "poscounts")
# Run DESeq2
ds_full = DESeq(ds_full, test = "Wald", fitType = "parametric")

# Extract results
res_full = results(ds_full, cooksCutoff = FALSE)
# Set p-value threshold and filter out significant OTUs
alpha = 0.055
sigTab_full = res_full[which(res_full$padj < alpha),]
sigTab_full = cbind(as(sigTab_full, "data.frame"), as(tax_table(physeq)[rownames(sigTab_full), ], "matrix"))
head(sigTab_full)

# Lets make a heatmap of the top differentially abundant microbes
hm.data = log10(otu_table(physeq)[rownames(sigTab_full),]+1e-10)
# Define metrics for clustering
drows1 <- "euclidean"
dcols1 <- "euclidean"
# Create annotations for heatmap
hm.rows = data.frame(Taxonomy = sigTab_full$Species)
hm.cols = data.frame(Response = gsub(".*_R","R",gsub(".*_NR","NR",colnames(otu_table(physeq)))))
hm.cols$Response = factor(hm.cols$Response,levels(hm.cols$Response)[c(2,1)])
rownames(hm.cols) = colnames(otu_table(physeq))
# Create heatmap
# type "?pheatmap()" for more help
# CairoWin()
library("pheatmap")
hm_full = pheatmap(hm.data,
                kmeans_k = NA,
                show_rownames = TRUE, show_colnames = TRUE,
                main = "Top differentially abundant OTUs, log10(relative abundance)",
                clustering_method = "average",
                cluster_rows = TRUE, cluster_cols = TRUE,
                clustering_distance_rows = drows1, 
                clustering_distance_cols = dcols1,
                labels_row = hm.rows$Taxonomy,
                annotation_col = hm.cols,
                fontsize_row = 10,
                cellwidth = 8,
                cellheight = 18,
                legend = F,
                annotation_legend = F
)

# Lets plot some violin plots comparing the relaticve abundance for 
# each of the significant differntially abundant OTUs.
library("reshape2")
# Store mOTUs IDs for all the significant OTUs
sigOTUs = rownames(sigTab_full)
# Convert OTU table to relative abundance
physeqR = transform_sample_counts(physeq, function(x) x / sum(x) )
# Subset OTU table in phyloseq object
physeq_subset = subset(otu_table(physeqR), rownames(otu_table(physeqR)) %in%  sigOTUs)
physeq_sig = merge_phyloseq(physeq_subset, tax_table(physeqR), sample_data(physeqR))
# Format data for ggplot
plotData = as.matrix(otu_table(physeq_sig))
plotData = melt(plotData)
plotData = cbind(plotData, Response = gsub(".*_R","R",gsub(".*_NR","NR",plotData$Var2)))
colnames(plotData) = c("mOTU","Patient","Abundance","Response")

library("ggplot2")
#plotData = plotData[plotData$Abundance != 0,]
plotData$Abundance = log10(as.numeric(plotData$Abundance) + 1e-8)

p_full = ggplot(data = plotData,aes(x=Response,y=Abundance,color=Response)) +
  facet_wrap(facets = "mOTU",ncol = 4, scales = "free") +
  geom_violin() +
  geom_jitter()
  
p_full

library("RColorBrewer")
library("ggplot2")