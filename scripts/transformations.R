
# get_transformation <- function(power, variable_name, dataframe) {
#   
#   # Function to add new column with the variable ("variable_name") values rised to the power of "power".
#   
#   var <- dataframe[, variable_name]
#   var_trans <- var^power # transformation
#   
#   dataframe[,paste0(variable_name, '_', power)] <- var_trans
#   
#   return(dataframe)
#   
# }

scaling <- function(x, shift=0, scale=1) {
  
  # Function to scale and shift the numerical vector.
  
  y <- (x+shift)/scale
  
  return(y)
}

get_transformation <- function(power, variable) {
  
  # Function to rise the variable values to the power of "power".
  
    var_trans <- variable^power # transformation
  
  # dataframe[,paste0(variable_name, '_', power)] <- var_trans
  
  return(var_trans)
  
}

get_trans_range <- function(variable_name, power_range, dataframe) {
  
  # Function to get the set of transformations for the range of powers and given variable ("variable_name").
  
  power_range <- setdiff(power_range, c(0))
  
  var <- dataframe[, variable_name]
  
  trans <- sapply(power_range, get_transformation, var) |> as.data.frame.matrix()
  
  colnames(trans) <- paste(variable_name, power_range, sep = '_')
  
  return(trans)
}
