# ==== IMPORT DATA ====
library("phyloseq")
# Import Phyloseq object for all studies
physeq = readRDS(file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/physeq.rds")

physeq  = transform_sample_counts(physeq, function(x) x / sum(x) )
physeq = filter_taxa(physeq, function(x) mean(x) > 1e-5, TRUE)
plot_bar(physeq,x="Response", fill = "Phylum")
