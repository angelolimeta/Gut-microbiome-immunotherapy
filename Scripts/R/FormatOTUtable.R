# == IMPORT DATA ==

# Remember to manually remove the "#" character in the heading of the .motus file
# This has already been done for the files that are imported below
Gop_OTU = read.delim("~/Documents/PhD/Gut-microbiome-immunotherapy/data/motus_files_merged/PRJEB22893.motus", stringsAsFactors = FALSE)
Mat_OTU = read.delim("~/Documents/PhD/Gut-microbiome-immunotherapy/data/motus_files_merged/PRJNA399742.motus", stringsAsFactors = FALSE)
Fra_OTU = read.delim("~/Documents/PhD/Gut-microbiome-immunotherapy/data/motus_files_merged/PRJNA397906.motus", stringsAsFactors = FALSE)

clin = readRDS(file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/clin.rds")

# == TIDY DATA ==

# Combine all data sets
OTU = cbind(Gop_OTU, Mat_OTU[,-c(1,2,3)], Fra_OTU[,-c(1,2,3)])

# Create phylogenetic metafile, which associates each mOTU to a specific taxon
phylo_meta = OTU[,c(1,2,3)]
rownames(phylo_meta) = phylo_meta[,1]
phylo_meta = phylo_meta[,-1]
phylo_meta = data.frame(lapply(phylo_meta, as.character), stringsAsFactors=FALSE)

# Only keep relative abundances in OTU data frame
rownames(OTU) = OTU[,1]
OTU = OTU[,-c(1,2,3)]

# Rename columns according to Patient_id, specified in the clin data frame
for (i in 1:length(colnames(OTU))) {
  if(sum(colnames(OTU)[i] == clin$Sample_id) == 1) {
    colnames(OTU)[i] = clin$Patient_id[colnames(OTU)[i] == clin$Sample_id] 
  }
}

# Remove the 5 repeat samples within 1 month of starting ICT from the Frankel et al data.
# See FormatMetadata.R for more info
OTU = OTU[,grepl("p",colnames(OTU))]

# == SAVE DATA ==

saveRDS(OTU, file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/OTU.rds")
saveRDS(phylo_meta, file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/phylo_meta.rds")