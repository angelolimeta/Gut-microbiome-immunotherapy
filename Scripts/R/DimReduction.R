library("ggplot2")
library("Rtsne")
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

# Filter out microbes that are not identified across at least 5 samples
physeq_f = filter_taxa(physeq, function(x) sum(x > 0) > 5, TRUE)

# ==== PCoA ====

# Calculate weighted UniFrac distances between all filtered samples
distance.matrix = UniFrac(physeq_f,weighted = TRUE)
# Perform PCoA
PCoA = cmdscale(distance.matrix, eig = TRUE, x.ret = TRUE)
# Extract variance information for each PCoA axis
PCoA.variance = round(PCoA$eig/sum(PCoA$eig)*100, 1)
# Extract values for each component
PCoA.values = PCoA$points
# Format data for ggplot
PCoA.data = data.frame(Sample = rownames(PCoA.values),
                       X = PCoA.values[,1],
                       Y = PCoA.values[,2],
                       Response = gsub(".*_R","R",gsub(".*_NR","NR",rownames(PCoA.values))),
                       Study = gsub(".*_Gop_.*","Gopalakrishnan",gsub(".*_Mat_.*","Matson",gsub(".*_Fra_.*","Frankel",rownames(PCoA.values)))))

# Call ggplot
p_MDS = ggplot(data = PCoA.data, aes(x = X, y = Y)) +
  geom_point(size = 3, aes(color = Response, shape = Study)) +
  stat_ellipse(aes(X,Y,color=Response, group=Response)) +
  xlab(paste("MDS1 - ", PCoA.variance[1],"%", sep = "")) +
  ylab(paste("MDS2 - ", PCoA.variance[2],"%", sep = "")) +
  ggtitle("All samples") +
  theme(plot.title = element_text(face="bold")) +
  theme(text = element_text(family = "Helvetica"))

# Perform the analysis for all studies separately
studies = c("Gopalakrishnan","Matson","Frankel")
# Initialize empty plot list
p_list = list()
for (i in 1:length(studies)) {
  local({ #Sets local enviroment in order to not overwrite plot data
    
    # Extract three letter study ID
    study_id = substr(studies[i], start = 1, stop = 3)
    # Subset OTU_filtered matrix
    physeq_f_study = prune_samples(samples = grepl(study_id,sample_names(physeq_f)), physeq_f)
    # Calculate weighted UniFrac distances between all filtered samples
    distance.matrix = UniFrac(physeq_f_study,weighted = TRUE)
    # Perform PCoA
    PCoA = cmdscale(distance.matrix, eig = TRUE, x.ret = TRUE)
    # Extract variance information for each PCoA axis
    PCoA.variance = round(PCoA$eig/sum(PCoA$eig)*100, 1)
    # Extract values for each component
    PCoA.values = PCoA$points
    # Format data for ggplot
    PCoA.data = data.frame(Sample = rownames(PCoA.values),
                           X = PCoA.values[,1],
                           Y = PCoA.values[,2],
                           Response = gsub(".*_R","R",gsub(".*_NR","NR",rownames(PCoA.values))))
    # Call ggplot
    p_MDS_study = ggplot(data = PCoA.data, aes(x = X, y = Y)) +
      geom_point(size = 3, aes(color = Response), show.legend = FALSE) +
      stat_ellipse(aes(X,Y,color=Response, group=Response), show.legend = FALSE) +
      xlab(paste("MDS1 - ", PCoA.variance[1],"%", sep = "")) +
      ylab(paste("MDS2 - ", PCoA.variance[2],"%", sep = "")) +
      ggtitle(paste(studies[i]," et al.", sep = "")) +
      theme(plot.title = element_text(face="bold")) +
      theme(text = element_text(family = "Helvetica"))
    # Save to plot list
    p_list[[i]] <<- p_MDS_study
  })
}

p_list[[4]] = p_MDS

quartz(width = 8, height = 6, pointsize = 10)
grid.arrange(
  grobs = p_list,
  layout_matrix = rbind(c(NA, 4, 4),
                        c(1, 2, 3))
)

# ==== tSNE ====

## https://github.com/opisthokonta/tsnemicrobiota
#library("devtools")
#install_github("opisthokonta/tsnemicrobiota")
library("tsnemicrobiota")

#Perform tSNE on filtered samples
set.seed(9)

tsne_res = tsne_phyloseq(physeq_f, distance='wunifrac',
                     perplexity = 20, verbose=0, rng_seed = 9)

tsne.data = data.frame(Sample = sample_names(physeq_f),
                       X = tsne_res$tsne$par[,1],
                       Y = tsne_res$tsne$par[,2],
                       Response = gsub(".*_R","R",gsub(".*_NR","NR",sample_names(physeq_f))),
                       Study = gsub(".*_Gop_.*","Gopalakrishnan",gsub(".*_Mat_.*","Matson",gsub(".*_Fra_.*","Frankel",sample_names(physeq_f)))))
# Call ggplot
p_tsne <- ggplot(tsne.data, aes(x=X, y=Y)) +
  geom_point(size = 2.5, aes(color = Response, shape = Study)) +
  stat_ellipse(aes(X,Y,color=Response, group=Response)) +
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  ggtitle("All samples") +
  theme(plot.title = element_text(face="bold")) +
  theme(text = element_text(family = "Helvetica"))

# Perform the analysis for all studies separately
studies = c("Gopalakrishnan","Matson","Frankel")
# Initialize empty plot list
p_list_tsne = list()
for (i in 1:length(studies)) {
  local({ #Sets local enviroment in order to not overwrite plot data
    
    # Extract three letter study ID
    study_id = substr(studies[i], start = 1, stop = 3)
    # Subset OTU_filtered matrix
    physeq_f_study = prune_samples(samples = grepl(study_id,sample_names(physeq_f)), physeq_f)
    # Perform tSNE
    tsne_res = tsne_phyloseq(physeq_f_study, distance='wunifrac',
                             perplexity = 20, verbose=0, rng_seed = 9)
    # Format data for ggplot
    tsne.data = data.frame(Sample = sample_names(physeq_f_study),
                           X = tsne_res$tsne$par[,1],
                           Y = tsne_res$tsne$par[,2],
                           Response = gsub(".*_R","R",gsub(".*_NR","NR",sample_names(physeq_f_study))))
    # Call ggplot
    p_tsne_study <- ggplot(tsne.data, aes(x=X, y=Y)) +
      geom_point(size = 2.5, aes(color = Response), show.legend = FALSE) +
      stat_ellipse(aes(X,Y,color=Response, group=Response), show.legend = FALSE) +
      xlab("tSNE 1") +
      ylab("tSNE 2") +
      ggtitle(paste(studies[i]," et al.", sep = "")) +
      theme(plot.title = element_text(face="bold")) +
      theme(text = element_text(family = "Helvetica"))
    # Save to plot list
    p_list_tsne[[i]] <<- p_tsne_study
  })
}

p_list_tsne[[4]] = p_tsne

quartz(width = 8, height = 6, pointsize = 10)
grid.arrange(
  grobs = p_list_tsne,
  layout_matrix = rbind(c(NA, 4, 4),
                        c(1, 2, 3))
)

# ==== SAVE PLOTS ====

# MDS
p_full = grid.arrange(
  grobs = p_list,
  layout_matrix = rbind(c(NA, 4, 4),
                        c(1, 2, 3))
)
ggsave(filename = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Figures/MDS.svg",
       plot = p_full,
       width = 8, 
       height = 6, 
       pointsize = 10)

# tSNE
p_full_tsne = grid.arrange(
  grobs = p_list_tsne,
  layout_matrix = rbind(c(NA, 4, 4),
                        c(1, 2, 3))
)
ggsave(filename = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Figures/tsne.svg",
       plot = p_full_tsne,
       width = 8, 
       height = 6, 
       pointsize = 10)