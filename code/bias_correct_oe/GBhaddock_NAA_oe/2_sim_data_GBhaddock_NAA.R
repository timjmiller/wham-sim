# Brian Stock
# June 15 2020
# Simulation test WHAM

# Assumes you open R in project directory
# source(here::here("code","bias_correct_oe","GBhaddock_NAA_oe","2_sim_data_GBhaddock_NAA.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","bias_correct_oe","GBhaddock_NAA_oe")
out_dir <- here("data","simdata","bias_correct_oe","GBhaddock_NAA_oe")
dir.create(out_dir, showWarnings=FALSE)

# Step 1: Fit 4 NAA models to 2019 SNE-MA yellowtail flounder data

# Step 2: Simulate 4 NAA operating models:
#  1. rec 	iid
#  2. rec 	ar1_y
#  3. rec+1 iid
#  4. rec+1 2dar1
# For each OM, simulate using
#  a. obs error only (fix NAA at estimated values)
#  b. obs + process error (new NAA)

# Step 3: Fit OMs to simulated datasets
# Step 4: Collect + plot results
# -----------------------------------------------------------------------

# Simulate data
n.mods <- 4
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
		# simdata[[i]] <- mod$simulate(par=mod$env$last.par.best, complete=TRUE)	
	}
	saveRDS(simdata, file=file.path(out_dir,paste0("simdata_om",m,".rds")))
}

# errors due to NAA_devs out of bounds in simulations?
#   R: malloc.c:4023: _int_malloc: Assertion `(unsigned long) (size) >= (unsigned long) (nb)' failed.

#   corrupted size vs. prev_size
#   Aborted (core dumped)

#   free(): invalid next size (fast)
#   Aborted (core dumped)
