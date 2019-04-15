##################################################
## Project: Gut-microbiome-immunotherapy
## Script purpose: Tidy up the metadata from the 4 studies
## Date: 2019-04-12
## Author: Angelo Limeta (angelol@chalmers.se)
##################################################

# == LOAD DATA ==

library(readr)
Gop_meta_full <- data.frame(read_delim("Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Gop_meta.csv", 
                       ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE),
                       stringsAsFactors = FALSE)

Mat_meta_full <- data.frame(read_csv("Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Mat_meta.csv"),
                       stringsAsFactors = FALSE)

Fra_meta_full <- data.frame(read_delim("Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Fra_meta.csv", 
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
colnames(Gop_meta) = c("Patient_id","Response","Treatment","Study")

# Mat_meta
Mat_meta = Mat_meta_full[,c("Sample","Response")]
Mat_meta = cbind(Mat_meta,rep("Anti-PD1",nrow(Mat_meta)),rep("Matson et al",nrow(Mat_meta)))
Mat_meta$Response = gsub("NonResponder","NR",Mat_meta$Response)
Mat_meta$Response = gsub("Responder","R",Mat_meta$Response)
colnames(Mat_meta) = c("Patient_id","Response","Treatment","Study")

# Fra_meta
Fra_meta = Fra_meta_full[,c("Patient.Identifier","RECIST.Category","ICT.Therapy")]
Fra_meta = cbind(Fra_meta,rep("Frankel et al",nrow(Fra_meta)))
colnames(Fra_meta) = c("Patient_id","Response","Treatment","Study")
Fra_meta$Response = gsub("Stable","NR",Fra_meta$Response)
Fra_meta$Response = gsub("Progression","NR",Fra_meta$Response)
Fra_meta$Response = gsub("Response","R",Fra_meta$Response)
Fra_meta$Treatment[which(Fra_meta$Treatment == "IN")] = "Anti-PD1 + Anti-CTLA4"
Fra_meta$Treatment[which(Fra_meta$Treatment == "I")] = "Anti-CTLA4"
Fra_meta$Treatment[which(Fra_meta$Treatment == "P")] = "Anti-PD1"
Fra_meta$Treatment[which(Fra_meta$Treatment == "N")] = "Anti-PD1"

