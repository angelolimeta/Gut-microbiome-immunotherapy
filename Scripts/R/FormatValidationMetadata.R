##################################################
## Project: Gut-microbiome-immunotherapy
## Script purpose: Tidy up the metadata from validation study
## Date: 2020-01-08
## Author: Angelo Limeta (angelol@chalmers.se)
##################################################

# == LOAD DATA ==
library(readr)

# Study metadata
val_clin <- data.frame(read_delim("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Pet_meta.tsv", 
                                       delim = ",",
                                       escape_double = FALSE,
                                       col_names = TRUE,
                                       trim_ws = TRUE),
                            stringsAsFactors = FALSE)
val_clin = val_clin[,-1]

val_sampInfo = data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/PRJNA541981.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE),
                          stringsAsFactors = FALSE)

# Link each patient with an individual sample
Sample_id = vector(mode="character", length=nrow(val_clin))
for (i in 1:nrow(val_clin)) {
    Sample_id[i] = val_sampInfo$run_accession[which(val_clin$BioSample[i] == val_sampInfo$sample_accession)]
}
val_clin = cbind(Sample_id,val_clin)

# Define response as PFS < 6 mo.
pfs = val_clin$months_to_progression
Response = ifelse(pfs > 6 | is.na(pfs),"R","NR")
val_clin = cbind(val_clin,Response)

# Remove non-baseline samples
val_clin = val_clin[val_clin$Time == "Baseline",]

# Create unique identifier for each patient, e.g p001_Gop_NR
Patient_id = vector(mode="character", length=nrow(val_clin))
for (i in 1:nrow(val_clin)) {
Patient_id[i] = paste("p",sprintf("%03d",i),"_Val_",val_clin$Response[i],sep = "")
}

val_clin = cbind(Patient_id,val_clin)
val_clin$Patient_id = as.character(val_clin$Patient_id)

# Save metadata
saveRDS(val_clin, file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/val_clin.rds")
