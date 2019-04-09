library(readr)
PRJEB22863 <- data.frame(read_delim("~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJEB22863.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE))
View(PRJEB22863)

PRJEB22863[,1] = paste("http://",PRJEB22863[,1],sep = "")

PRJEB22863 = PRJEB22863[,1]
colnames(PRJEB22863) = NULL

write.table(PRJEB22863, "~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/PRJEB22863.txt", sep="\n")
