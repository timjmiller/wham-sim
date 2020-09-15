# Brian Stock
# June 15 2020
# Simulation test WHAM

# Assumes you open R in project directory
# source(here::here("code","bias_correct_oepe","GBhaddock_sel_oepe","2_sim_data_GBhaddock_sel.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","bias_correct_oepe","GBhaddock_sel_oepe")
out_dir <- here("data","simdata","bias_correct_oepe","GBhaddock_sel_oepe")
dir.create(out_dir, showWarnings=FALSE)

# Simulate data
n.mods <- 3
n.sim <- 100

# Load all models in list first
mod.list <- file.path(res_dir,paste0("m",1:n.mods,".rds"))
mods2 <- lapply(mod.list, readRDS)
mods1 <- lapply(mod.list, readRDS)
for(m in 1:n.mods) mods1[[m]]$env$data$simulate_state <- rep(0,4) # simulate_state = 0 (fixed NAA, obs error only)

sim.seeds <- readRDS(here("data","sim_seeds.rds"))
for(m in 1:n.mods){
	simdata <- vector("list",n.sim)
	for(i in 1:n.sim){
		print(paste0("Model: ",m," Sim: ", i))
		set.seed(sim.seeds[i])
		simdata[[i]][[1]] <- mods1[[m]]$simulate(par=mods2[[m]]$env$last.par.best, complete=TRUE)
		simdata[[i]][[2]] <- mods2[[m]]$simulate(par=mods2[[m]]$env$last.par.best, complete=TRUE)		
		# simdata[[i]] <- mod$simulate(par=mod$env$last.par.best, complete=TRUE)	
	}
	saveRDS(simdata, file=file.path(out_dir,paste0("simdata_om",m,".rds")))
}

# # check m3 is simulating selAA 
# # fixed effect pars controlling re (should always be same)
# simdata[[1]][[1]]$sel_repars 
# simdata[[2]][[1]]$sel_repars 
# simdata[[1]][[2]]$sel_repars 
# simdata[[2]][[2]]$sel_repars 

# # re should be same for type 1 
# cbind(simdata[[1]][[1]]$selpars_re, simdata[[2]][[1]]$selpars_re)

# # re should be different for type 2
# cbind(simdata[[1]][[2]]$selpars_re, simdata[[2]][[2]]$selpars_re)

# simdata[[1]][[2]]$selAA[[1]][,1:4]
# simdata[[2]][[2]]$selAA[[1]][,1:4]

# # ----------------------------------------------------------------
# # check m2 is simulating selAA 
# simdata <- readRDS("/home/bstock/Documents/ms/wham-sim/data/simdata/bias_correct_oepe/GBhaddock_sel_oepe/simdata_om2.rds")

# # fixed effect pars controlling re (should always be same)
# simdata[[1]][[1]]$sel_repars 
# simdata[[2]][[1]]$sel_repars 
# simdata[[1]][[2]]$sel_repars 
# simdata[[2]][[2]]$sel_repars 

# # re should be same for type 1 
# cbind(simdata[[1]][[1]]$selpars_re, simdata[[2]][[1]]$selpars_re)

# # re should be different for type 2
# cbind(simdata[[1]][[2]]$selpars_re, simdata[[2]][[2]]$selpars_re)

# simdata[[1]][[2]]$selAA[[1]][,1:4]
# simdata[[2]][[2]]$selAA[[1]][,1:4]

# # ----------------------------------------------------------------
# # check m1 is NOT simulating different selAA 
# simdata <- readRDS("/home/bstock/Documents/ms/wham-sim/data/simdata/bias_correct_oepe/GBhaddock_sel_oepe/simdata_om1.rds")

# # re should be same for type 1 
# cbind(simdata[[1]][[1]]$selpars_re, simdata[[2]][[1]]$selpars_re)

# # re should be same for type 2
# cbind(simdata[[1]][[2]]$selpars_re, simdata[[2]][[2]]$selpars_re)

# simdata[[1]][[2]]$selAA[[1]][,1:4]
# simdata[[2]][[2]]$selAA[[1]][,1:4]
