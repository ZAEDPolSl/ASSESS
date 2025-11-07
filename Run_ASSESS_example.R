rm(list = ls())
gc()

library(fastcmprsk) 
library(ggplot2)

#!!! Please set your working directory
setwd("c:/POLSL/Analysis/Software/ASSESS/")
source("scripts/ASSESS_predict.R")
source("scripts/extract_model_vars.R")
source("scripts/transformations.R")

# load clinical data for given patient
pat_name <- "Patient2_TNBC"
pat_data <- read.delim(file=paste0("data/",pat_name,".txt"))

# load ASSESS model parameters
load("scripts/ASSESS_BREAST_2010_2020_mdls.RData")
load("scripts/Pooled_HR.RData")

grp_col <- c("#888888","#000000","#156c11","#ee6c4d","#b399fd","#0070C0","#FFC000","#A58114")
names(grp_col) <- c("None","Black","Green","Red","Purple","Blue","Olap","Pembro")

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

# plot survival curves
data_plot <- reshape::melt(res, id.vars="Time")
colnames(data_plot) <- c("Time","Name","Probability")
data_plot$Name <- as.character(data_plot$Name)

# get Treatment and Survival endpoint from Name variable
tmp <- strsplit(data_plot$Name, "_")
data_plot$Treatment <- sapply(tmp,function(x){x[2]})
data_plot$Treatment[is.na(data_plot$Treatment)] <- "None"
data_plot$Surv_type <- sapply(tmp,function(x){x[1]})

# change to factors
data_plot$Treatment <- factor(data_plot$Treatment, levels=names(grp_col))
data_plot$Surv_type <- factor(data_plot$Surv_type, levels=c("BCSS","OS"))

p <- ggplot(data_plot, aes(x=Time, y=Probability, col=Treatment)) + 
  geom_line() + scale_color_manual(values=grp_col) +
  facet_grid(~Surv_type) + 
  geom_line() + labs(x="Survival time [months]") + theme_bw()
print(p)
pdf(file=paste0("res/Survival_",pat_name,".pdf"), width=8, height=6)
print(p)
dev.off()

cat("Patient info:")
print(pat_data)