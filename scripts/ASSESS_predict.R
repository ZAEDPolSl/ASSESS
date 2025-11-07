ASSESS_predict <- function(Age, Tsize, Tgrade, Pnodes, HR, HER2, endpoint, mdl_all, HR_pool_all=NULL, max_time=80){
# Prediction of breast-cancer free and overall survival using ASSESS model
# Inputs:
# Age - patient age at diagnosis (20 - 85)
# Tsize - tumor size (0-200 mm)
# Tgrade - tumor grade (1-3)
# Pnodes - no. positive nodes (0-100)
# HR - hormone receptor status (P or N)
# HER2 - HER2 status (P or N)
# endpoint - survival endpoint ("OS" - overall survival or "BCSS" - breast cancer specific survival)
#

  require(fastcmprsk)

  # define model features
  inModel <- c("Age","TsizeMerged","Tgrade","Pnodes","ChemoT","Rad") # all variables
  inModel_cat <- c("Tgrade","ChemoT","Rad") # categorical variables
  inModel_num <- paste0(setdiff(inModel, inModel_cat),'_scaled') # numerical variables
  
  # define subtype
  if (HR == "P") {
    if (HER2 == "N") subtype = "HR+/HER2-" else if (HER2 == "P") subtype = "HR+/HER2+"
  } else if (HR == "N") {
    if (HER2 == "N") subtype = "TNBC" else if (HER2 == "P") subtype = "HR-/HER2+"
  }

  # scale continuous variables and shift
  Age_scaled <- scaling(Age, scale = 100)
  TsizeMerged_scaled <- scaling(Tsize, scale = 10)
  Pnodes_scaled <- scaling(Pnodes, shift = 1)
  
  # model selection
  mdl <- mdl_all[grepl(subtype, names(mdl_all), fixed = T)] # for subtype
  
  if (endpoint == 'BCSS') { # for endpoint
    mdl <- mdl[grepl('BCSS', names(mdl))]
  } else {
    mdl <- mdl[!grepl('BCSS', names(mdl))]
  }
  
  names(mdl) <- sapply(strsplit(names(mdl),'_',fixed = T), '[[', 1)
  
  # model variables
  vars <- get_model_variables(mdl[[1]])

  # apply subtype/endpoint-specific MFP transformations
  powers <- get_powers(vars)
  if (is.matrix(powers)) {powers <- as.list(as.data.frame(powers))}

  trans_age <- sapply(powers$Age_scaled, get_transformation, variable=Age_scaled)
  names(trans_age) <- paste0('Age_scaled_', powers$Age_scaled)
  
  trans_tsize <- sapply(powers$TsizeMerged_scaled, get_transformation, variable=TsizeMerged_scaled)
  names(trans_tsize) <- paste0('TsizeMerged_scaled_', powers$TsizeMerged_scaled)
  
  trans_pnodes <- sapply(powers$Pnodes_scaled, get_transformation, variable=Pnodes_scaled)
  names(trans_pnodes) <- paste0('Pnodes_scaled_', powers$Pnodes_scaled)
  
  # prepare categorical variables (model matrix)
  inModel_modmat <- unique(grep(paste(inModel_cat,collapse="|"), vars, value=TRUE)) # categorical variables as model matrix names
  model_mat <- rep(0, times=length(inModel_modmat)); names(model_mat) <- inModel_modmat
  model_mat[paste0('Tgrade',Tgrade)] <- 1
  
  # create data table
  data <- c(model_mat, trans_age, trans_pnodes, trans_tsize)
  data <- data[vars]

  # predict cumulative incidence
  pred_BC <- predict(mdl[[setdiff(1:2,grep('non', names(mdl)))]], data, getBootstrapVariance=F, tL=1)
  pred_nonBC <- predict(mdl[[grep('non', names(mdl))]], data, getBootstrapVariance=F, tL=1)
  
  Pred_final <- data.frame(1:max_time)
  colnames(Pred_final) <- "Time"
  
  if (endpoint == 'OS') {
    Pred_final$OS <- 1 - (spline(pred_BC$ftime,pred_BC$CIF,xout=Pred_final$Time)$y + 
                            spline(pred_nonBC$ftime,pred_nonBC$CIF,xout=Pred_final$Time)$y)
    Pred_final$OS[Pred_final$OS < 0] <- 0
    Pred_final$OS[Pred_final$OS > 1] <- 1
  } else if (endpoint == 'BCSS') {
    Pred_final$BCSS <- 1 - spline(pred_BC$ftime,pred_BC$CIF,xout=Pred_final$Time)$y
    Pred_final$BCSS[Pred_final$BCSS < 0] <- 0
    Pred_final$BCSS[Pred_final$BCSS > 1] <- 1
  }
  
  # Add chemotherapy effect
  # load(paste0(HRpreds,".RData"))
  if (!is.null(HR_pool_all)) {
    if (endpoint == 'OS') {
      HR <- HR_pool_all$OS[,subtype]
    } else if (endpoint == 'BCSS') {
      HR <- HR_pool_all$DFS[,subtype]
    }
    un_grp <- names(HR)
    
    for (a in 1:length(un_grp)){
      tmp <- Pred_final[,2] ^ HR[un_grp[a]]
      tmp[tmp < 0] <- 0
      tmp[tmp > 1] <- 1
      Pred_final[,paste0(endpoint,"_",un_grp[a])] <- tmp
    }
    
    # res <- list()
    # res$Curves <- Pred_final
    # res$Estimates <- Surv_final
  }

  return(Pred_final)
}