library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("Rtsne")
library("dplyr")
## if not installed, quickly add it as follows:
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("RColorBrewer", "pheatmap", "ggplot2", "Rtsne"))

# ==== IMPORT DATA ====

# Abundance table for each OTU
OTU = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/OTU.rds")
# Clinical metadata
clin = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/clin.rds")
# Pyholgenetic metadata
phylo_meta = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/phylo_meta.rds")

# ==== DATA PRE-PROCESSING ====
# Remove Routy et al samples
clin = clin %>% filter(Study != "Routy_et_al")
OTU = OTU[,colnames(OTU) %in% clin$Patient_id]

# Convert to relative abundance
OTU = apply(OTU, 2, function(x){x/sum(x)})

# Filter out microbes that are not identified across at least 5 samples
OTU_filtered = OTU[apply(OTU, 1, function(x) sum(x > 0)) > 5,]
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
# Define metrics for clustering
drows1 <- "euclidean"
dcols1 <- "euclidean"
# Create annotations for heatmap
hm.rows = data.frame(Taxonomy = as.character(gsub(".*s__","",diffOTUs$taxonomy)))
hm.cols = data.frame(Response = gsub(".*_R","R",gsub(".*_NR","NR",colnames(OTU))))
hm.cols$Response = factor(hm.cols$Response,levels(hm.cols$Response)[c(2,1)])
rownames(hm.cols) = colnames(OTU_filtered)
# Create heatmap
# type "?pheatmap()" for more help
# CairoWin()
test = pheatmap(hm.data,
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

# Perform the analysis for all studies separately
studies = c("Gopalakrishnan","Matson","Frankel")
# Initialize empty plot list
p_list = list()
for (i in 1:length(studies)) {
  local({ #Sets local enviroment in order to not overwrite plot data
    # Extract three letter study ID
    study_id = substr(studies[i], start = 1, stop = 3)
    # Pre-process each individual study
    # Subset the full OTU matrix
    OTU_filtered_study = OTU[,grepl(study_id,colnames(OTU))]
    # Filter out microbes that are not identified across at least 3 samples
    OTU_filtered_study = OTU_filtered_study[apply(OTU_filtered_study, 1, function(x) sum(x > 0)) > 3,]
    # Sort according to most abundant OTUs
    OTU_filtered_study = OTU_filtered_study[order(rowSums(OTU_filtered_study),decreasing = TRUE),]
    # Initialise empty data frame for storing differentially abundant OTUs
    diffOTUs_study = data.frame(pvalue = 1:nrow(OTU_filtered_study), fold_change = 1:nrow(OTU_filtered_study), taxonomy = phylo_meta[rownames(OTU_filtered_study),1])
    rownames(diffOTUs_study) = rownames(OTU_filtered_study)
    # Loop over all OTUs
    for (j in 1:nrow(OTU_filtered_study)) {
      # Extract abundances of each OTU for both R and NR
      R_abundance = as.numeric(OTU_filtered_study[j,grepl("_R",colnames(OTU_filtered_study))])
      NR_abundance = as.numeric(OTU_filtered_study[j,grepl("_NR",colnames(OTU_filtered_study))])
      # Perform wilcoxon rank sum test between responders and non-responders
      diffOTUs_study$pvalue[j] = wilcox.test(R_abundance,NR_abundance, alternative = "two.sided")$p.value
      diffOTUs_study$fold_change[j] = log2((median(R_abundance) + 0.001)/(median(NR_abundance) + 0.001))
    }
    ## Adjust for multiple testing
    #diffOTUs_study$pvalue = p.adjust(diffOTUs_study$pvalue,method = "fdr")
    # Order according to most significant
    diffOTUs_study = diffOTUs_study[order(diffOTUs_study$pvalue),]
    # Keep only significant microbes (P < 0.05)
    diffOTUs_study = diffOTUs_study[diffOTUs_study$pvalue < 0.05,]
    # Remove patients with s
    # Lets make a heatmap of the top differentially abundant microbes
    hm.data = log10(OTU_filtered_study[rownames(diffOTUs_study),]+1e-10)
    # define metrics for clustering
    drows1 <- "euclidean"
    dcols1 <- "euclidean"
    # Create annotations for heatmap
    hm.rows = data.frame(Taxonomy = as.character(gsub(".*s__","",diffOTUs_study$taxonomy)))
    hm.cols = data.frame(Response = gsub(".*_R","R",gsub(".*_NR","NR",colnames(OTU_filtered_study))))
    hm.cols$Response = factor(hm.cols$Response,levels(hm.cols$Response)[c(2,1)])
    rownames(hm.cols) = colnames(OTU_filtered_study)
    # Create heatmap
    p_study = pheatmap(hm.data,
                    kmeans_k = NA,
                    show_rownames = TRUE, show_colnames = TRUE,
                    main = paste(studies[i]," et al.", sep = ""),
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
    
    # Save to plot list
    p_list[[i]] <<- p_study
  })
}

quartz()
p_list[[3]]

# ==== PLOT INDIVIDUAL OTUs ====
OTUforPlot = "ref_mOTU_v2_4880"
plotData = data.frame(Abundance = log10(as.numeric(OTU[OTUforPlot,])+1e-8), Response = gsub(".*_R","R",gsub(".*_NR","NR",colnames(OTU))))

# Call ggplot
p3 = ggplot(data = plotData, aes(x=Response,y=Abundance, color=Response)) + geom_boxplot() +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="red") +
  labs(y = "log10 normalized relative abundance", x = "") +
  labs(title = as.character(gsub(".*s__","",phylo_meta[OTUforPlot,1])))

p3

# Call ggplot
p4 = ggplot(data = plotData, aes(x=Response,y=Abundance, color=Response)) + geom_violin() + 
  labs(title=as.character(gsub(".*s__","",phylo_meta[OTUforPlot,1])), 
       x="",
       y="Relative abundance")
quartz()
p4
# ==== Phylum Level analysis ====

# ==== Diversity indices ====

# Extract phyla
a = OTU_filtered[1,2]
a = gsub(".*p__", "",a)
a = gsub("\\|c__.*","",a)

row.names(OTU_filtered) = OTU_filtered[,2]
OTU_filtered = OTU_filtered[,-c(1,2,3)]


