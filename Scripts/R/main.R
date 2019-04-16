library("pheatmap")
library("RColorBrewer")
library("ggplot2")
## if not installed, quickly add it as follows:
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("RColorBrewer", "pheatmap"))

# ==== IMPORT DATA ====

# Abundance table for each OTU
OTU = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/OTU.rds")
# Clinical metadata
clin = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/clin.rds")
# Pyholgenetic metadata
phylo_meta = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/phylo_meta.rds")

# ==== DATA PRE-PROCESSING ====

# Filter out microbes that are not identified across at least 3 samples
OTU_filtered = OTU[apply(OTU, 1, function(x) sum(x > 0)) > 3,]
# Sort according to most abundant OTUs
OTU_filtered = OTU_filtered[order(rowSums(OTU_filtered),decreasing = TRUE),]



# ==== OTU Level analysis ====

# Initialise empty data frame for storing differentially abundant OTUs
diffOTUs = data.frame(pvalue = 1:nrow(OTU_filtered), fold_change = 1:nrow(OTU_filtered), taxonomy = phylo_meta[rownames(OTU_filtered),1])
rownames(diffOTUs) = rownames(OTU_filtered)

# Loop over all OTUs
for (i in 1:nrow(OTU_filtered)) {
  # Extract abundances of each OTU for both R and NR
  R_abundance = as.numeric(OTU_filtered[i,grepl("_R",colnames(OTU_filtered))])
  NR_abundance = as.numeric(OTU_filtered[i,grepl("_NR",colnames(OTU_filtered))])
  
  diffOTUs$pvalue[i] = wilcox.test(R_abundance,NR_abundance, alternative = "two.sided")$p.value
  diffOTUs$fold_change[i] = log2((median(R_abundance) + 0.001)/(median(NR_abundance) + 0.001))
}

## Adjust for multiple testing
#diffOTUs$pvalue = p.adjust(diffOTUs$pvalue,method = "fdr")

# Order according to most significant
diffOTUs = diffOTUs[order(diffOTUs$pvalue),]

# Keep only significant microbes (P < 0.05)
diffOTUs = diffOTUs[diffOTUs$pvalue < 0.05,]

# Lets make a heatmap of the top differentially abundant microbes
hm.data = log10(OTU_filtered[rownames(diffOTUs),]+1e-10)

# define metrics for clustering
drows1 <- "correlation"
dcols1 <- "correlation"

hm.rows = data.frame(Taxonomy = as.character(gsub(".*s__","",diffOTUs$taxonomy)))

hm.cols = data.frame(Response = gsub(".*_R","R",gsub(".*_NR","NR",colnames(OTU))))
rownames(hm.cols) = colnames(OTU_filtered)

# Create heatmap
# type "?pheatmap()" for more help
# CairoWin()
test = pheatmap(hm.data,
         kmeans_k = NA,
         show_rownames = TRUE, show_colnames = TRUE,
         main = "Top differentially abundant OTUs, -log10(relative abundance)",
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
         annotation_legend = T
)

# ==== PLOT INDIVIDUAL OTUs ====
OTUforPlot = "ref_mOTU_v2_3198"
plotData = data.frame(Abundance = log10(as.numeric(OTU[OTUforPlot,])+1e-8)+8, Response = gsub(".*_R","R",gsub(".*_NR","NR",colnames(OTU))))

p1 = ggplot(data = plotData, aes(x=Response,y=Abundance, color=Response)) + geom_boxplot() +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="red") +
  labs(y = "Relative abundance", x = "") +
  labs(title = as.character(gsub(".*s__","",phylo_meta[OTUforPlot,1])))
quartz()
p1


p2 = ggplot(data = plotData, aes(x=Response,y=Abundance, color=Response)) + geom_violin() + 
  labs(title=as.character(gsub(".*s__","",phylo_meta[OTUforPlot,1])), 
       x="",
       y="Relative abundance")
quartz()
p2
# ==== Phylum Level analysis ====

# ==== Diversity indices ====

# Extract phyla
a = OTU_filtered[1,2]
a = gsub(".*p__", "",a)
a = gsub("\\|c__.*","",a)

row.names(OTU_filtered) = OTU_filtered[,2]
OTU_filtered = OTU_filtered[,-c(1,2,3)]


