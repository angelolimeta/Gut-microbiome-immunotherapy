#PRJEB22863 data set

# import dowload links
library(readr)
links_all <- as.data.frame(read_csv("~/Documents/PhD/Gut-microbiome-immunotherapy/Scripts/download/PRJEB22863/PRJEB22863.txt", col_names = FALSE))

# import links in Dir
filesInDir <- as.data.frame(read_csv("~/Documents/PhD/Gut-microbiome-immunotherapy/Scripts/download/PRJEB22863/filesInDir.txt", col_names = FALSE))
filesInDir <- filesInDir[-179,]

# import sample annotation
data <- as.data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/Scripts/download/PRJEB22863/allData.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE))
data$run_accession = paste(data$run_accession,".fastq.gz",sep = "")

# Keep only T0 samples
data_t0 = data[grepl("T0",data$run_alias),]

# Store T1 & T2 samples
data_t1t2 = data[!grepl("T0",data$run_alias),]

# Samples to download
SamplesToDownload = setdiff(data_t0$run_accession,filesInDir)
length(SamplesToDownload)

LinksToDownload = data_t0[data_t0$run_accession %in% SamplesToDownload, 3]
LinksToDownload = paste("http://",LinksToDownload,sep = "")


# Samples to remove
SamplesToRemove = intersect(data_t1t2$run_accession,filesInDir)
length(SamplesToRemove)

# All links
AllLinks = data_t0$fastq_ftp
AllLinks = paste("http://",AllLinks,sep = "")

# Save data
setwd("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Scripts/download/PRJEB22863")
write.table(SamplesToRemove, file = "SamplesToRemove.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(LinksToDownload, file = "LinksToDownload.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(AllLinks, file = "AllLinks.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
