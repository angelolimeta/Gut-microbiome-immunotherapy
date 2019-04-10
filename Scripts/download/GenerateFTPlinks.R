library(readr)

#

# PRJEB22863
PRJEB22863 <- data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJEB22863_meta.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE))
View(PRJEB22863)

PRJEB22863[,1] = paste("http://",PRJEB22863[,1],sep = "")

PRJEB22863 = PRJEB22863[,1]
colnames(PRJEB22863) = NULL

writeLines(PRJEB22863, "~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJEB22863.txt")

# PRJNA399742

PRJNA399742 <- data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJNA399742_meta.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE))

PRJNA399742 = t(data.frame(strsplit(PRJNA399742[,1],split = ";")))

rownames(PRJNA399742) = NULL

PRJNA399742 = c(PRJNA399742[,1], PRJNA399742[,2])

PRJNA399742 = paste("http://",PRJNA399742,sep = "")


writeLines(PRJNA399742, "~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJNA399742.txt")

# PRJNA397906

PRJNA397906 <- data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJNA397906_meta.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE))

PRJNA397906 = t(data.frame(strsplit(PRJNA397906[,1],split = ";")))

rownames(PRJNA397906) = NULL

PRJNA397906 = c(PRJNA397906[,1], PRJNA397906[,2])

PRJNA397906 = paste("http://",PRJNA397906,sep = "")


writeLines(PRJNA397906, "~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJNA397906.txt")

# PRJEB22893

PRJEB22893 <- data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJEB22893_meta.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE))

PRJEB22893 = t(data.frame(strsplit(PRJEB22893[,1],split = ";")))

rownames(PRJEB22893) = NULL

PRJEB22893 = c(PRJEB22893[,1], PRJEB22893[,2])

PRJEB22893 = paste("http://",PRJEB22893,sep = "")


writeLines(PRJEB22893, "~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJEB22893.txt")