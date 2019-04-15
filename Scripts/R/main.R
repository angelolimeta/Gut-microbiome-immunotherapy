# == IMPORT DATA ==

# Remember to manually remove the "#" character in the heading of the .motus file
# This has already been done for the files that are imported below
Gop_OTU = read.delim("~/Documents/PhD/Gut-microbiome-immunotherapy/data/motus_files_merged/PRJEB22893.motus", stringsAsFactors = FALSE)
Mat_OTU = read.delim("~/Documents/PhD/Gut-microbiome-immunotherapy/data/motus_files_merged/PRJNA399742.motus", stringsAsFactors = FALSE)
Fra_OTU = read.delim("~/Documents/PhD/Gut-microbiome-immunotherapy/data/motus_files_merged/PRJNA397906.motus", stringsAsFactors = FALSE)

# == TIDY DATA ==

# Combine all data sets
OTU = cbind(Gop_OTU, Mat_OTU[,-c(1,2,3)], Fra_OTU[,-c(1,2,3)])

# Create phylogenetic metafile, which associates each mOTU to a specific taxon
phylo_meta = OTU[,c(1,2,3)]
rownames(phylo_meta) = phylo_meta[,1]
phylo_meta = phylo_meta[,-1]

# Only keep relative abundances in OTU data frame
rownames(OTU) = OTU[,1]
OTU = OTU[,-c(1,2,3)]

# == DATA PRE-PROCESSING ==

# Filter out microbes that are not identified across at least 3 samples
OTU_filtered = OTU[apply(OTU, 1, function(x) sum(x > 0)) > 3,]
# Sort according to most abundant OTUs
OTU_filtered = OTU_filtered[order(rowSums(OTU_filtered),decreasing = TRUE),]

# Extract phyla
a = OTU_filtered[1,2]
a = gsub(".*p__", "",a)
a = gsub("\\|c__.*","",a)

row.names(OTU_filtered) = OTU_filtered[,2]
OTU_filtered = OTU_filtered[,-c(1,2,3)]


