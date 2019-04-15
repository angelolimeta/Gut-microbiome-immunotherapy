PRJNA399742 <- read.delim("~/Documents/PhD/Gut-microbiome-immunotherapy/Scripts/download/PRJNA399742.txt")

samplesToRemove = PRJNA399742[PRJNA399742$library_strategy == "AMPLICON",]

samplesToRemove = as.character(samplesToRemove$run_accession)

samplesToRemove_1 = paste(samplesToRemove,"_1.fastq.gz",sep = "")
samplesToRemove_2 = paste(samplesToRemove,"_2.fastq.gz",sep = "")

samplesToRemove = c(samplesToRemove_1, samplesToRemove_2)

writeLines(samplesToRemove, "~/Documents/PhD/Gut-microbiome-immunotherapy/files2download/samplesToRemove.txt")
