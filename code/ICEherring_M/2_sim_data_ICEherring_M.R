# Brian Stock
# June 15 2020
# Simulation test WHAM

# Assumes you open R in project directory
# source(here::here("code","ICEherring_M","2_sim_data_ICEherring_M.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","ICEherring_M")
out_dir <- here("data","simdata","ICEherring_M")
dir.create(out_dir, showWarnings=FALSE)

# Simulate data
n.mods <- 3
n.sim <- 100

# Load all models in list first
mod.list <- file.path(res_dir,paste0("m",1:n.mods,".rds"))
mods2 <- lapply(mod.list, readRDS)
mods1 <- lapply(mod.list, readRDS)
for(m in 1:n.mods) mods1[[m]]$env$data$simulate_state <- rep(0,4) # simulate_state = 0 (fixed NAA, obs error only)

# set.seed(12345)
# sim.seeds = sample(1:1000000, n.sim, replace = FALSE) # set random seeds (re-use for each model)
# saveRDS(sim.seeds, here("data","sim_seeds.rds"))
sim.seeds <- readRDS(here("data","sim_seeds.rds"))
for(m in 1:n.mods){
	simdata <- vector("list",n.sim)
	for(i in 1:n.sim){
		print(paste0("Model: ",m," Sim: ", i))
		set.seed(sim.seeds[i])
		simdata[[i]][[1]] <- mods1[[m]]$simulate(par=mods2[[m]]$env$last.par.best, complete=TRUE)
		simdata[[i]][[2]] <- mods2[[m]]$simulate(par=mods2[[m]]$env$last.par.best, complete=TRUE)		
	}
	saveRDS(simdata, file=file.path(out_dir,paste0("simdata_om",m,".rds")))
}


