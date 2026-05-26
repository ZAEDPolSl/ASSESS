rm(list = ls())
gc()

library(fastcmprsk) 
library(ggplot2)

#!!! Please set your working directory
setwd("c:/YALE/Projects/Lajos/Assess/Run_ASSESS/v6_16_04_2026/")
source("scripts/ASSESS_predict.R")
source("scripts/extract_model_vars.R")
source("scripts/transformations.R")

# load clinical data for given patient
pat_name <- "Patient7_P_P"
pat_data <- read.delim(file=paste0("data/",pat_name,".txt"))

# load ASSESS model parameters
load("scripts/ASSESS_BREAST_2010_2020_mdls.RData")
load("scripts/Pooled_HR.RData")

# survival endpoint
surv_end <- c("OS","BCSS")

# run ASSESS prediction
# (Age, Tsize, Tgrade, Pnodes, HR, HER2, endpoint, mdl_all, HR_pool_all, max_time=80)
res_OS <- ASSESS_predict(pat_data$Age, pat_data$Tsize, pat_data$Tgrade, pat_data$Pnodes, 
                      pat_data$HR, pat_data$HER2, "OS", mdl_all, HR_pool_all, max_time=96)
res_BCSS <- ASSESS_predict(pat_data$Age, pat_data$Tsize, pat_data$Tgrade, pat_data$Pnodes, 
                           pat_data$HR, pat_data$HER2, "BCSS", mdl_all, HR_pool_all, max_time=96)
res <- cbind(res_OS, res_BCSS[,2:ncol(res_BCSS)])

# save prediction
write.table(res, file=paste0("res/",pat_name,".txt"),sep="\t",row.names=F)

