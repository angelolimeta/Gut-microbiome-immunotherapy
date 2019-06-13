library("ggplot2")
library("phyloseq")
library("ape")
library("gridExtra")
library("svglite") # Requires Cairo. Run "brew install cairo" if on a mac

## if not installed, quickly add it as follows:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq", version = "3.8")


# ==== IMPORT DATA ====

# Phyloseq object (Both abundance and taxonomic info)
physeq = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/Abundance_tables/physeq.rds")

# ==== DATA PRE-PROCESSING ====

# Filter out microbes that are not identified
physeq_f = filter_taxa(physeq, function(x) sum(x > 0) > 0, TRUE)

# ==== ANALYSE ALPHA DIVERSITY ====

# Plot alpha diversity
p_full = plot_richness(physeq_f,x = "Response", measures = "Shannon", color = "Response") + 
  geom_boxplot() +
  ylab("") +
  xlab("")
p_full
#Calculate alpha diversity
alpha_div = estimate_richness(physeq_f, measures = "Shannon") 

# Extract diversities for responders and non responders
R_div = as.numeric(alpha_div[grepl("_R",rownames(alpha_div)),1])
NR_div = as.numeric(alpha_div[grepl("_NR",rownames(alpha_div)),1])
# Perform wilcoxon rank sum test
wilcox.test(R_div, NR_div, alternative = "two.sided")$p.value

# Perform the analysis for all studies separately
studies = c("Gopalakrishnan","Matson","Frankel")

# Define diversity indeces to be used
# - Shannon: Estimator of species richness and species evenness, more weight on species richness.
# - InvSimpson: Estimator of species richness and species evenness, more weight on species evenness.
# - ACE: Abundance-based Coverage Estimator of species richness.
# - Chao1: Abundance-based estimator of species richness.

div_index = c("Shannon","InvSimpson","ACE","Chao1")

# Create empty matrix for storing pvalues for significance testing between R and NR
pvalue_mat = data.frame(matrix(data = NA,nrow = length(studies), ncol = length(div_index)))
rownames(pvalue_mat) = studies
colnames(pvalue_mat) = div_index
# Initialize empty plot list
p_list = list()
for (i in 1:length(studies)) {
  local({ #Sets local enviroment in order to not overwrite plot data
    # Extract three letter study ID
    study_id = substr(studies[i], start = 1, stop = 3)
    # Subset OTU_filtered matrix
    physeq_f_study = prune_samples(samples = grepl(study_id,sample_names(physeq_f)), physeq_f)
    for (j in 1:length(div_index)) {
      #Calculate alpha diversity
      alpha_div_study = estimate_richness(physeq_f_study, measures = div_index[j]) 
      # Extract diversities for responders and non responders
      R_div_study = as.numeric(alpha_div_study[grepl("_R",rownames(alpha_div_study)),1])
      NR_div_study = as.numeric(alpha_div_study[grepl("_NR",rownames(alpha_div_study)),1])
      # Perform wilcoxon rank sum test
      pvalue_mat[i,j] <<- wilcox.test(R_div_study, NR_div_study, alternative = "two.sided")$p.value
    }
    # Create y-axis label for grid.arrange
    if(i==1){label = "Index value"}
    else{label = ""}
    # Plot alpha diversity (Shannon only)
    p_study = plot_richness(physeq_f_study,x = "Response", measures = "Shannon", color = "Response") +
      #geom_boxplot() +
      theme(legend.position = "none") +
      ylab(label) +
      xlab("")
    # Save to plot list
    p_list[[i]] <<- p_study
  })
}

p_list[[4]] = p_full 

p_all = grid.arrange(grobs = p_list,
                     layout_matrix = rbind(c(1, 2, 3, 4),
                                           c(1, 2, 3, 4))
)

# ==== ANALYSE ALPHA DIVERSITY ====

# ==== SAVE DATA ====
ggsave(filename = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Figures/Diversity.svg",
       plot = p_all,
       width = 8, 
       height = 4, 
       pointsize = 10)
