# ==== LOAD DEPENDENCIES ====

library("tidyverse")
library("data.table")
## if not installed, quickly add it as follows:
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("tidyverse", "data.table"))

# ==== IMPORT DATA ====

# Gene abundance table, values are in RPK
geneAbundance = fread("~/Documents/PhD/Gut-microbiome-immunotherapy/data/humann2/all_genefamilies.tsv", sep= "\t") %>% as_tibble

# Pathway abundance table
pathAbundance = as.data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/data/humann2/all_pathabundance.tsv", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE))

# Pathway coverage table
pathCoverage = as.data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/data/humann2/all_pathcoverage.tsv", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE))

# ==== DATA PRE-PROCESSING ====

###################
## geneAbundance ##
###################

# Let's start by converting our gene matrix into long format.
geneAbundance %>% gather(key = "sample", value = "abundance", -`# Gene Family`) -> geneAbundance_long
# Remove # character from column name
geneAbundance_long %>% rename(gene_family = `# Gene Family`) -> geneAbundance_long
# Convert into CoPM
geneAbundance_long %>% group_by(sample) %>% mutate(CoPM = abundance/sum(abundance, na.rm = T)*1e6) -> geneAbundance_long

# SET PREVALENCE THREASHOLD PARAMS:
PREVALENCE_THR = 0.8 # gene should be present in 80% of all samples
ABUNDANCE_THR = 1e-5 # genes with an RPK value of above 1e-5 count as present

# This command runs the following:
# - Create a temporary column called isPresent that tells us whether a gene is present or not.
# - Group data frame according to each gene
# - Calculate prevalence across samples for each gene, using the isPresent column.
geneAbundance_long %>% mutate(isPresent = ifelse(abundance > ABUNDANCE_THR, 1, 0)) %>% 
  group_by(gene_family) %>% summarize(prevalence = sum(isPresent)/length(isPresent)) -> geneAbundance_prevalence

# Let's keep the prevalent genes and convert back into wide format
geneAbundance %>% semi_join(geneAbundance_prevalence %>% filter(prevalence > PREVALENCE_THR ), by= c("# Gene Family" = "gene_family")) -> geneAbundance_filtered

# The resulting data frame has 61152 genes
tail(geneAbundance_filtered)

# Rename the columns
colnames(geneAbundance_filtered) = gsub("se_","",colnames(geneAbundance_filtered))
colnames(geneAbundance_filtered) = gsub("\\_.*","",colnames(geneAbundance_filtered))

# Add rownames
geneAbundance_filtered = as.data.frame(geneAbundance_filtered)
rownames(geneAbundance_filtered) = geneAbundance_filtered[,1]
geneAbundance_filtered = geneAbundance_filtered[,-1]

###################
## pathAbundance ##
###################

# SET PREVALENCE THREASHOLD PARAMS:
PATH_PREVALENCE_THR = 0.2

rownames(pathAbundance) = pathAbundance$`# Pathway`
pathAbundance = pathAbundance[,-1]
colnames(pathAbundance) = gsub("se_","",colnames(pathAbundance))
colnames(pathAbundance) = gsub("\\_.*","",colnames(pathAbundance))

# Convert to CoPM
pathAbundance = apply(pathAbundance, 2, function(x){(x*1e6)/sum(x)})

# Filter out pathways not found in less than 20% of samples
pathAbundance_filtered = as.data.frame(pathAbundance[apply(pathAbundance, 1, function(x) sum(x > ABUNDANCE_THR)) > PATH_PREVALENCE_THR*ncol(pathAbundance),])

##################
## pathCoverage ##
##################

rownames(pathCoverage) = pathCoverage$`# Pathway`
pathCoverage = pathCoverage[,-1]
colnames(pathCoverage) = gsub("se_","",colnames(pathCoverage))
colnames(pathCoverage) = gsub("\\_.*","",colnames(pathCoverage))

# Filter out pathways not found in less than 20% of samples
pathCoverage_filtered = as.data.frame(pathCoverage[apply(pathCoverage, 1, function(x) sum(x > 0)) > PATH_PREVALENCE_THR*ncol(pathCoverage),])


# ==== ANNOTATE SAMPLES ====

# Remove the 5 repeat samples from Frankel et al data
SamplesToRemove = colnames(pathAbundance)[which(!(colnames(pathAbundance) %in% clin$Sample_id))]

geneAbundance_filtered = geneAbundance_filtered[,-which(colnames(geneAbundance_filtered) %in% SamplesToRemove)]
pathAbundance_filtered = pathAbundance_filtered[,-which(colnames(pathAbundance_filtered) %in% SamplesToRemove)]
pathCoverage_filtered = pathCoverage_filtered[,-which(colnames(pathCoverage_filtered) %in% SamplesToRemove)]

# Load clinical metadata
clin = readRDS("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/clin.rds")

SampleNames = 1:ncol(geneAbundance_filtered)

for (i in 1:length(SampleNames)) {
  patient_ID = clin$Patient_id[which(clin$Sample_id %in% colnames(geneAbundance_filtered)[i])] 
  SampleNames[i] = patient_ID
}

colnames(geneAbundance_filtered) = SampleNames
colnames(pathAbundance_filtered) = SampleNames
colnames(pathCoverage_filtered) = SampleNames

# ==== SAVE DATA ====
saveRDS(geneAbundance_filtered, file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/humann2/geneAbundance_filtered.rds")
saveRDS(pathAbundance_filtered, file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/humann2/pathAbundance_filtered.rds")
saveRDS(pathCoverage_filtered, file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/data/humann2/pathCoverage_filtered.rds")



