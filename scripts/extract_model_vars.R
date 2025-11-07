get_model_variables <- function(model) {
  
  # Function to extract the variables in the Fine-Gray model.
  
  inModel_num <- c("Age_scaled","TsizeMerged_scaled","Pnodes_scaled")
  
  # extract variables:
  variables <- colnames(model$df)[-(1:2)]
  variables <- sapply(strsplit(variables, '...', fixed = T),'[[',2)
  variables <- gsub('_.','_-',variables, fixed = T)
  
  # extract powers:
  
  powers <- sapply(inModel_num, function(var) {
    tmp <- strsplit(variables[grepl(var,variables)],'_',fixed=T)
    sapply(tmp,'[[', unique(lengths(tmp))) |> as.numeric()
    })
  
  return(variables)
}

get_powers <- function(variables) {
  
  # Function to extract the powers to each the numarical variables in the model should be raised.
  
  inModel_num <- c("Age_scaled","TsizeMerged_scaled","Pnodes_scaled")
  
  # extract powers:
  
  powers <- sapply(inModel_num, function(var) {
    tmp <- strsplit(variables[grepl(var,variables)],'_',fixed=T)
    sapply(tmp,'[[', unique(lengths(tmp))) |> as.numeric()
  })
  
  return(powers)
}
