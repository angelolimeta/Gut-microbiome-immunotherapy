# This script converts metadata, obtained from NCBI Biosample in xml format, to a 
# more human-readable .csv format. This also allows for easier downstream 
# processing of the metadata.
function(filePath, save.tsv = FALSE, file.tsv = NULL){
  
  require(xml2)
  # Read xml file
  meta = read_xml("/Users/angelol/Documents/PhD/Gut-microbiome-immunotherapy/Metadata/Rou_meta.xml")
  # Convert to list
  meta_list = as_list(meta)
  
  # Initialize empty data frame
  nSamples = length(meta_list$BioSampleSet)
  nAttributes = 1 + length(meta_list$BioSampleSet$BioSample$Description) +
    length(meta_list$BioSampleSet$BioSample$Attributes)
  meta_df = data.frame(matrix(NA, nrow = nSamples, ncol = nAttributes))
  
  # Fill data frame with values from .xml file
  index = 1
  for (i in 1:nSamples) {
    # Store BioSample ID
    meta_df[i,index] = meta_list$BioSampleSet[[i]]$Ids[[1]][[1]]
    index = index + 1
    # Store Biosample Description
    for (j in 1:length(meta_list$BioSampleSet$BioSample$Description)) {
      meta_df[i,index] = meta_list$BioSampleSet[[i]]$Description[[j]][[1]][[1]]
      index = index + 1
    }
    # Store Biosample Attributes
    for (j in 1:length(meta_list$BioSampleSet$BioSample$Attributes)) {
      meta_df[i,index] = meta_list$BioSampleSet[[i]]$Attributes[[j]][[1]][[1]]
      index = index + 1
    }
    index = 1
  }
  
  # Add colnames
  index = 1
  colnames(meta_df)[index] = attributes(meta_list$BioSampleSet[[i]]$Ids[[1]])[[1]]
  index = index +1
  
  for (j in 1:length(meta_list$BioSampleSet$BioSample$Description)) {
    colnames(meta_df)[index] = names(meta_list$BioSampleSet[[i]]$Description)[[j]]
    index = index + 1
  }
  
  for (j in 1:length(meta_list$BioSampleSet$BioSample$Attributes)) {
    colnames(meta_df)[index] = attributes(meta_list$BioSampleSet[[i]]$Attributes[[j]])[[1]]
    index = index + 1
  }
  
  if (save.tsv = TRUE) {
    if (is.null(file.tsv)){
      return("Please specify path to .tsv file!")
      break
    }
    write.table(meta_df, file = file.tsv, quote = FALSE, sep = "\t")
  }
  
  return(meta_df)
}




