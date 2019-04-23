library(phyloseq)
library(ape)

# ==== IMPORT DATA ====

# Abundance table for each OTU
OTU = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/OTU.rds")
# Clinical metadata
clin = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/clin.rds")
# Pyholgenetic metadata
phylo_meta = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/phylo_meta.rds")
# Tree file
mOTUs_tree = read.tree(file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Tree_files/mOTUs.treefile")

# ==== CONVERT INTO PHYLOSEQ ====

# Create OTU table
OTU_phylo = otu_table(OTU, taxa_are_rows = TRUE)
# Create empty taxonomy matrix
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(OTU), ncol = 7)
rownames(taxmat) <- rownames(OTU)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Fill each values in the taxonomy matrix
for (i in 1:nrow(taxmat)) {
  full_taxonomy = phylo_meta[rownames(taxmat)[i],1] # Get the full taxonomy for each OTU
  
  taxmat[i,1] = gsub("\\|p.*","",gsub(".*k__","",full_taxonomy)) # Domain
  taxmat[i,2] = gsub("\\|c.*","",gsub(".*p__","",full_taxonomy)) # Phylum
  taxmat[i,3] = gsub("\\|o.*","",gsub(".*c__","",full_taxonomy)) # Class
  taxmat[i,4] = gsub("\\|f.*","",gsub(".*o__","",full_taxonomy)) # Order
  taxmat[i,5] = gsub("\\|g.*","",gsub(".*f__","",full_taxonomy)) # Family
  taxmat[i,6] = gsub("\\|s.*","",gsub(".*g__","",full_taxonomy)) # Genus
  taxmat[i,7] = gsub(".*s__","",full_taxonomy) # Species
}
# Convert into taxonomy table
TAX_phylo = tax_table(taxmat)
# Create phyloseq object
physeq = phyloseq(OTU_phylo, TAX_phylo)
# Add phylogenetic tree info
physeq = merge_phyloseq(physeq,mOTUs_tree)
# Add clinical data
rownames(clin) = clin$Patient_id
clin = clin[,-2]
sampledata = sample_data(clin)
physeq = merge_phyloseq(physeq,sampledata)


# ==== SAVE DATA ====

saveRDS(physeq, file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/physeq.rds")
