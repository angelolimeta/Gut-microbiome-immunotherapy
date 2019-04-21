##################################################
## Project: Gut-microbiome-immunotherapy
## Script purpose: Tidy up the metadata from the 4 studies
## Date: 2019-04-12
## Author: Angelo Limeta (angelol@chalmers.se)
##################################################

# == LOAD DATA ==

# Study metadata
library(readr)
Gop_meta_full <- data.frame(read_delim("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Gop_meta.csv", 
                       ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE),
                       stringsAsFactors = FALSE)

Mat_meta_full <- data.frame(read_csv("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Mat_meta.csv"),
                       stringsAsFactors = FALSE)

Fra_meta_full <- data.frame(read_delim("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Fra_meta.csv", 
                        "\t", escape_double = FALSE, trim_ws = TRUE),
                       stringsAsFactors = FALSE)
# Sample data
Gop_sampInfo = data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/PRJEB22893.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE),
                        stringsAsFactors = FALSE)

Mat_sampInfo = data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/PRJNA399742.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE),
                        stringsAsFactors = FALSE)

Fra_sampInfo = data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/PRJNA397906.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE),
                        stringsAsFactors = FALSE)

# == TIDY DATA ==

# Gop_meta
colnames(Gop_meta_full)[1] = "Age"
Gop_meta_full[,1] = gsub(".*=","",Gop_meta_full[,1])
for (i in 2:ncol(Gop_meta_full)) {
  colnames(Gop_meta_full)[i] = gsub("=.*","",Gop_meta_full[1,i])
  Gop_meta_full[,i] = gsub(".*=","",Gop_meta_full[,i])
}
# Extract common data
Gop_meta = Gop_meta_full[,c("subject_id","phenotype","treatment")]
Gop_meta = cbind(Gop_meta,rep("Gopalakrishnan et al",nrow(Gop_meta)))
colnames(Gop_meta) = c("Study_patient_id","Response","Treatment","Study")

# Mat_meta
Mat_meta = Mat_meta_full[,c("Sample","Response")]
Mat_meta = cbind(Mat_meta,rep("Anti-PD1",nrow(Mat_meta)),rep("Matson et al",nrow(Mat_meta)))
Mat_meta$Response = gsub("NonResponder","NR",Mat_meta$Response)
Mat_meta$Response = gsub("Responder","R",Mat_meta$Response)
colnames(Mat_meta) = c("Study_patient_id","Response","Treatment","Study")
# Patient 09, 41 and 42 seems to only have Amplicon data. Exclude them from the dataset
Mat_meta = Mat_meta[Mat_meta$Study_patient_id != "P09",]
Mat_meta = Mat_meta[Mat_meta$Study_patient_id != "P41",]
Mat_meta = Mat_meta[Mat_meta$Study_patient_id != "P42",]

# Fra_meta
Fra_meta = Fra_meta_full[,c("Patient.Identifier","RECIST.Category","ICT.Therapy")]
Fra_meta = cbind(Fra_meta,rep("Frankel et al",nrow(Fra_meta)))
colnames(Fra_meta) = c("Study_patient_id","Response","Treatment","Study")
Fra_meta$Response = gsub("Stable","NR",Fra_meta$Response)
Fra_meta$Response = gsub("Progression","NR",Fra_meta$Response)
Fra_meta$Response = gsub("Response","R",Fra_meta$Response)
Fra_meta$Treatment[which(Fra_meta$Treatment == "IN")] = "Anti-PD1 + Anti-CTLA4"
Fra_meta$Treatment[which(Fra_meta$Treatment == "I")] = "Anti-CTLA4"
Fra_meta$Treatment[which(Fra_meta$Treatment == "P")] = "Anti-PD1"
Fra_meta$Treatment[which(Fra_meta$Treatment == "N")] = "Anti-PD1"

# Combine all metadata
clin = rbind(Gop_meta,Mat_meta,Fra_meta)
rownames(clin) = NULL

# Create unique identifier for each patient, e.g p001_Gop_NR
Patient_id = vector(mode="character", length=nrow(clin))
for (i in 1:nrow(clin)) {
  if(clin$Study[i] == "Gopalakrishnan et al"){
    # sprintf("%03d",i) formats numbers as fixed width, with leading zeros
    Patient_id[i] = paste("p",sprintf("%03d",i),"_Gop_",clin$Response[i],sep = "")
  }
  if(clin$Study[i] == "Matson et al"){
    Patient_id[i] = paste("p",sprintf("%03d",i),"_Mat_",clin$Response[i],sep = "")
  }
  if(clin$Study[i] == "Frankel et al"){
    Patient_id[i] = paste("p",sprintf("%03d",i),"_Fra_",clin$Response[i],sep = "")
  }
  
}

clin = cbind(Patient_id,clin)

# Remove non-WGS samples from Matson et al data
Mat_sampInfo = Mat_sampInfo[Mat_sampInfo$library_strategy == "WGS",]
Mat_sampInfo = Mat_sampInfo[,-2]

# Discard the 5 repeat samples within 1 month of starting ICT from the Frankel et al data.
Fra_sampInfo$sample_alias = gsub("melanoma_","",Fra_sampInfo$sample_alias)
Fra_sampInfo = Fra_sampInfo[!grepl("_",Fra_sampInfo$sample_alias),]



# Link each patient with an individual sample
Sample_id = vector(mode="character", length=nrow(clin))
for (i in 1:nrow(clin)) {
  if(clin$Study[i] == "Gopalakrishnan et al"){
    Sample_id[i] = Gop_sampInfo$run_accession[which(clin$Study_patient_id[i] == Gop_sampInfo$sample_alias)]
  }
  if(clin$Study[i] == "Matson et al"){
    Sample_id[i] = Mat_sampInfo$run_accession[which(gsub("P","",clin$Study_patient_id[i]) == gsub(".*_","",Mat_sampInfo$sample_alias))]
  }
  if(clin$Study[i] == "Frankel et al"){
    Sample_id[i] = Fra_sampInfo$run_accession[which(clin$Study_patient_id[i] == Fra_sampInfo$sample_alias)]
  }
}
clin = cbind(Sample_id,clin)

# Convert all columns into character vectors
clin = data.frame(lapply(clin, as.character), stringsAsFactors=FALSE)
# Response can be kept as a factor
clin$Response = as.factor(clin$Response)

# == SAVE DATA ==

# Metadata (clinical) is saved as a .rds file (clin.rds)
saveRDS(clin, file = "/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Processed_metadata/clin.rds")
