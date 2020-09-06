# Brian Stock
# Aug 17 2020
# Simulation test WHAM
# SNEMAYT
# M

# Assumes you open R in project directory
# source(here::here("code","bias_correct_oepe","SNEMAYT_M_oepe","2_sim_data_SNEMAYT_M_breaks_merge.R"))

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
# library(tidyverse)

# res_dir <- here("results","SNEMAYT_M")
out_dir <- here("data","simdata","bias_correct_oepe","SNEMAYT_M_oepe")
# dir.create(out_dir, showWarnings=FALSE)

# 1: Fit models
# 2: Simulate operating models
# 3: Fit OMs to simulated datasets
# 4: Collect + plot results
# -----------------------------------------------------------------------
n.mods <- 3
n.b <- 4
n.sim <- 100
# ind.sim <- 1:(n.sim/n.b) + (b-1)*(n.sim/n.b)
for(om in 1:n.mods){
	simdata <- list()
	for(b in 1:n.b){
		simdata_b <- readRDS(file.path(out_dir,paste0("simdata_om",om,"_",b,".rds")))
		simdata <- c(simdata, simdata_b)
	}
	saveRDS(simdata, file=file.path(out_dir,paste0("simdata_om",om,".rds")))
}

