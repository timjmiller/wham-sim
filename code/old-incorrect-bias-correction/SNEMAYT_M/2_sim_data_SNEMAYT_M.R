# Brian Stock
# Aug 17 2020
# Simulation test WHAM
# SNEMAYT
# M

# Assumes you open R in project directory
# source(here::here("code","SNEMAYT_M","2_sim_data_SNEMAYT_M.R"))

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","SNEMAYT_M")
out_dir <- here("data","simdata","SNEMAYT_M")
dir.create(out_dir, showWarnings=FALSE)

# 1: Fit models
# 2: Simulate operating models
# 3: Fit OMs to simulated datasets
# 4: Collect + plot results
# -----------------------------------------------------------------------

# Load fit models
# mod.list <- file.path(res_dir,paste0("m",1:n.mods,".rds"))
mod.list <- list.files(res_dir, full.names = TRUE)
mods2 <- lapply(mod.list, readRDS)
mods1 <- lapply(mod.list, readRDS)
n.mods <- length(mod.list)
for(m in 1:n.mods) mods1[[m]]$env$data$simulate_state <- rep(0,4) # simulate_state = 0 (fixed NAA, obs error only)

# set.seed(12345)
# sim.seeds = sample(1:1000000, n.sim, replace = FALSE) # set random seeds (re-use for each model)
# saveRDS(sim.seeds, here("data","sim_seeds.rds"))
n.sim <- 100
sim.seeds <- readRDS(here("data","sim_seeds.rds"))
for(m in 1:n.mods){
	simdata <- vector("list",n.sim)
	for(i in 1:n.sim){
		print(paste0("Model: ",m," Sim: ", i))
		set.seed(sim.seeds[i])
		simdata[[i]][[1]] <- mods1[[m]]$simulate(par=mods2[[m]]$env$last.par.best, complete=TRUE)
		simdata[[i]][[2]] <- mods2[[m]]$simulate(par=mods2[[m]]$env$last.par.best, complete=TRUE)		
	}
	saveRDS(simdata, file=file.path(out_dir,paste0("simdata_SNEMAYT_M_om",m,".rds")))
}

# # -------------------------------------------------------
# # Check simulated data
# m=1 # no random effects
# simdata <- readRDS(file.path(out_dir,paste0("simdata_SNEMAYT_M_om",m,".rds")))
# simdata[[1]][[1]]$MAA # obs error only (MAA constant)
# simdata[[1]][[2]]$MAA # oe + pe (MAA should not be constant)

# simdata[[1]][[1]]$M_repars # M re pars
# simdata[[1]][[2]]$M_repars # M re pars

# simdata[[1]][[1]]$M_re # M re
# simdata[[1]][[2]]$M_re # M re

# mods2[[m]]$parList$M_repars
# mods2[[m]]$parList$M_re

# m=2 # iid M_re
# simdata <- readRDS(file.path(out_dir,paste0("simdata_SNEMAYT_M_om",m,".rds")))
# simdata[[1]][[1]]$MAA 
# simdata[[1]][[2]]$MAA 

# # M re pars will match whether data simulated with/without process error (fixed effects controlling random effects)
# mods2[[m]]$parList$M_repars
# simdata[[1]][[1]]$M_repars 
# simdata[[1]][[2]]$M_repars

# mods2[[m]]$parList$M_re
# simdata[[1]][[1]]$M_re # M re simulated without process error should match model fit M_re (above)
# simdata[[1]][[2]]$M_re # M re simulated with process error should differ

# m=3 # 2d ar1 M_re
# simdata <- readRDS(file.path(out_dir,paste0("simdata_SNEMAYT_M_om",m,".rds")))
# simdata[[1]][[1]]$MAA 
# simdata[[1]][[2]]$MAA 

# # M re pars will match whether data simulated with/without process error (fixed effects controlling random effects)
# mods2[[m]]$parList$M_repars
# simdata[[1]][[1]]$M_repars 
# simdata[[1]][[2]]$M_repars

# mods2[[m]]$parList$M_re
# simdata[[1]][[1]]$M_re # M re simulated without process error should match model fit M_re (above)
# simdata[[1]][[2]]$M_re # M re simulated with process error should differ

