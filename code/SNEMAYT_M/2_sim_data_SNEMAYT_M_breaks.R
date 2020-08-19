# Brian Stock
# Aug 17 2020
# Simulation test WHAM
# SNEMAYT
# M

# Assumes you open R in project directory
# source(here::here("code","SNEMAYT_M","2_sim_data_SNEMAYT_M_breaks.R"))

args = commandArgs(trailingOnly=TRUE)
b = as.integer(args[1]) # break
om = as.integer(args[2])

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
# library(here)
# library(tidyverse)

res_dir <- out_dir <- getwd()
# res_dir <- here("results","SNEMAYT_M")
# out_dir <- here("data","simdata","SNEMAYT_M")
# dir.create(out_dir, showWarnings=FALSE)

# 1: Fit models
# 2: Simulate operating models
# 3: Fit OMs to simulated datasets
# 4: Collect + plot results
# -----------------------------------------------------------------------

# Load fit models
# mod.list <- file.path(res_dir,paste0("m",1:n.mods,".rds"))
mod.list <- file.path(res_dir,paste0("m",1:3,".rds"))
# mod.list <- list.files(res_dir, full.names = TRUE)
mods2 <- lapply(mod.list, readRDS)
mods1 <- lapply(mod.list, readRDS)
n.mods <- length(mod.list)
for(m in 1:n.mods) mods1[[m]]$env$data$simulate_state <- rep(0,4) # simulate_state = 0 (fixed NAA, obs error only)

# set.seed(12345)
# sim.seeds = sample(1:1000000, n.sim, replace = FALSE) # set random seeds (re-use for each model)
# saveRDS(sim.seeds, here("data","sim_seeds.rds"))
n.b <- 2
n.sim <- 100
ind.sim <- 1:(n.sim/n.b) + (b-1)*(n.sim/n.b)
# sim.seeds <- readRDS(here("data","sim_seeds.rds"))
sim.seeds <- readRDS("sim_seeds.rds")
# for(m in 1:n.mods){
	# simdata <- vector("list",n.sim)
	simdata <- vector("list",length(ind.sim))
	for(i in 1:length(ind.sim)){
		print(paste0("Model: ",om," Sim: ", ind.sim[i]))
		set.seed(sim.seeds[ind.sim[i]])
		simdata[[i]][[1]] <- mods1[[om]]$simulate(par=mods2[[om]]$env$last.par.best, complete=TRUE)
		simdata[[i]][[2]] <- mods2[[om]]$simulate(par=mods2[[om]]$env$last.par.best, complete=TRUE)		
	}
	saveRDS(simdata, file=file.path(out_dir,paste0("simdata_SNEMAYT_M_om",m,"_",b,".rds")))
}

# # fit model with random effects possibly on NAA, MAA, selAA, and Ecov
# # input = list with $data, $par, $map, $random
# mod <- TMB::MakeADFun(input)
# mod$opt <- stats::nlminb(mod$par, mod$fn, mod$gr)
# mod$sdrep <- TMB::sdreport(mod)

# # simulate observation error only
# mod$env$data$simulate_state <- rep(0,4) # see wham_v0.cpp
# simdata_oe <- mod$simulate(par=mod$env$last.par.best, complete=TRUE)

# # simulate obs + process error
# mod$env$data$simulate_state <- rep(1,4) # turn on process error in NAA, MAA, selAA, and Ecov
# simdata_oepe <- mod$simulate(par=mod$env$last.par.best, complete=TRUE)

# # fit model to data simulated with obs + process error
# input2 <- input
# input2$data <- simdata_oepe
# mod2 <- fit_wham(input2)
