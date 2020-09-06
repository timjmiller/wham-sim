# collect results when saved one sim at a time

# source("/home/bstock/Documents/ms/wham-sim/code/SNEMAYT_M/4_collect_sims_ind.R")

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
# library(here)
# library(tidyverse)

# res_dir <- "/home/bstock/Documents/NRC/code/remoteR/SNEMAYT_M"
res_dir <- out_dir <- "/home/bstock/Documents/ms/wham-sim/results/SNEMAYT_M"

n.sim <- 100
n.types <- 2
simdata <- readRDS(file.path(res_dir,"simdata_om1.rds"))
n.years <- simdata[[1]][[1]][["n_years_model"]]

nested.list <- function(len) if(length(len) == 1) vector("list", len) else lapply(1:len[1], function(...) nested.list(len[-1]))
# sdreps <- nested.list(c(n.types, n.sim)) # sdreps with analytical bias correction
# reps <- nested.list(c(n.types, n.sim))
# res.colnames <- c("om","em","type","year","sim","F_fit","F_fit_bc","F_sim","relF_fit","relF_fit_bc","relF_sim","SSB_fit","SSB_fit_bc","SSB_sim","relB_fit","relB_fit_bc","relB_sim","catch_fit","catch_fit_bc","catch_sim",paste0("NAA",1:6),paste0("NAA",1:6,"_bc"))
# results <- rep(list(rep(list(matrix(NA, ncol = length(res.colnames), nrow = n.years)),n.sim)),n.types) # nested lists with preallocated matrices
# results <- matrix(NA, ncol = length(res.colnames), nrow = n.years)

for(om in 1:3){
	for(em in 1:3){
		sdreps <- nested.list(c(n.types, n.sim)) # sdreps with analytical bias correction
		reps <- nested.list(c(n.types, n.sim))
		res.colnames <- c("om","em","type","year","sim","F_fit","F_fit_bc","F_sim","relF_fit","relF_fit_bc","relF_sim","SSB_fit","SSB_fit_bc","SSB_sim","relB_fit","relB_fit_bc","relB_sim","catch_fit","catch_fit_bc","catch_sim",paste0("NAA",1:6),paste0("NAA",1:6,"_bc"))
		results <- rep(list(rep(list(matrix(NA, ncol = length(res.colnames), nrow = n.years)),n.sim)),n.types) # nested lists with preallocated matrices
		# results <- matrix(NA, ncol = length(res.colnames), nrow = n.years)
		# n.sim = 50
		# if(om == 3 & em == 2) n.sim = 1
		# if(om == 3 & em == 3) n.sim = 48
		for(i in 1:n.sim){
			print(paste0("Sim ",i))
			res <- readRDS(file.path(res_dir,paste0("results_SNEMAYT_M_om",om,"_em",em,"_sim",i,".rds")))
			sdrep <- readRDS(file.path(res_dir,paste0("sdreps_SNEMAYT_M_om",om,"_em",em,"_sim",i,".rds")))
			rep <- readRDS(file.path(res_dir,paste0("reps_SNEMAYT_M_om",om,"_em",em,"_sim",i,".rds")))
			for(ty in 1:n.types){
				results[[ty]][[i]] <- res[[ty]]
				sdreps[[ty]][[i]] <- sdrep[[ty]]
				reps[[ty]][[i]] <- rep[[ty]]
				# results <- rbind(results, res[[ty]])
			}
			rm(list=c('res','sdrep','rep'))
			# rm('res')
		}
		# results <- as.data.frame(results)
		# results <- sapply(results, as.numeric)
		# results <- as.data.frame(results)
		saveRDS(results, file=file.path(out_dir,paste0("results_SNEMAYT_M_om",om,"_em",em,".rds")))		
		saveRDS(sdreps, file=file.path(out_dir,paste0("sdreps_SNEMAYT_M_om",om,"_em",em,".rds")))		
		saveRDS(reps, file=file.path(out_dir,paste0("reps_SNEMAYT_M_om",om,"_em",em,".rds")))		
	}
}


# saveRDS(sdreps, file=file.path(out_dir,paste0("sdreps_SNEMAYT_M_om",om,"_em",em,".rds")))
# saveRDS(reps, file=file.path(out_dir,paste0("reps_SNEMAYT_M_om",om,"_em",em,".rds")))
