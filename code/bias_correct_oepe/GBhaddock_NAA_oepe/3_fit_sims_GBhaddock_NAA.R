# Brian Stock
# June 15 2020
# Simulation test WHAM
# Step 3: Fit OMs to simulated datasets
#   NEFSC server version 2 (error catching)

# the following files should be copied to working directory:
#   simdata_omXX.rds (XX from 1-4)
#   mXX_input.rds (XX from 1-4)
#   sim_seeds.rds

# will write 3 files to working directory:
#   results_omXX_emYY.rds
#   sdreps_omXX_emYY.rds
#   reps_omXX_emYY.rds

args = commandArgs(trailingOnly=TRUE)
om = as.integer(args[1]) # operating model
em = as.integer(args[2]) # estimating model

# -----------------------------------------------------------------------
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
# library(here)
# library(tidyverse)

res_dir <- simdata_dir <- getwd()
# res_dir <- here("results","NScod_NAA")
# simdata_dir <- here("data","simdata","NScod_NAA")

# functions to calculate F/F40 and B/B40
calc_relF <- function(mod, sdrep, type="fit", bias.cor=TRUE){ # sdrep is summary(sdreport)
	if(type == "fit"){
		ind.FXSPR <- which(rownames(sdrep) == "log_FXSPR")
		if(bias.cor) F.t <- sdrep[ind.FXSPR,3] else F.t <- sdrep[ind.FXSPR,1]
		ind.faa <- which(rownames(sdrep) == "log_FAA_tot")
		n.yrs <- length(mod$years_full)
		n.ages <- mod$env$data$n_ages
  		if(bias.cor) faa <- matrix(sdrep[ind.faa,3], n.yrs, n.ages) else faa <- matrix(sdrep[ind.faa,1], n.yrs, n.ages)	
	}
	if(type == "sim"){
		F.t <- mod$env$data$log_FXSPR
		faa <- log(mod$env$data$FAA_tot)			
	}
  	age.full.f <- apply(faa,1, function(x) max(which(x == max(x))))
  	full.f <- faa[cbind(seq_along(age.full.f),age.full.f)]
	rel.f <- exp(full.f - F.t)
	return(rel.f)
}
calc_relB <- function(mod, sdrep, type="fit", bias.cor=TRUE){ # sdrep is summary(sdreport)
	if(type == "fit"){
		ind.SSB.FXSPR <- which(rownames(sdrep) == "log_SSB_FXSPR")
		if(bias.cor) SSB.t <- exp(sdrep[ind.SSB.FXSPR,3]) else SSB.t <- exp(sdrep[ind.SSB.FXSPR,1])
		ind.ssb <- which(rownames(sdrep) == "log_SSB")
  		if(bias.cor) ssb <- exp(sdrep[ind.ssb,3]) else ssb <- exp(sdrep[ind.ssb,1])
	}
	if(type == "sim"){
		SSB.t <- exp(mod$env$data$log_SSB_FXSPR)
		ssb <- mod$env$data$SSB
	}
	rel.ssb <- ssb / SSB.t	
	return(rel.ssb)
}
sdrep_out <- function(mod, sdrep, var, bias.cor=TRUE){
	if(var == "log_NAA_rep"){
		ind <- which(rownames(sdrep) == var)
		n.yrs <- length(mod$years_full)
		n.ages <- mod$env$data$n_ages
		if(bias.cor) x <- matrix(sdrep[ind,3], n.yrs, n.ages) else x <- matrix(sdrep[ind,1], n.yrs, n.ages)
	} else {
		ind <- which(rownames(sdrep) == var)
		if(bias.cor) x <- sdrep[ind,3] else x <- sdrep[ind,1]
	}
	return(exp(x)) # assumed on log scale
}
calc_results <- function(om, em, type, sim, fit1, s1){
	n.ages <- fit1$env$data$n_ages
	df <- as.matrix(data.frame(om=om, em=em, type=type, year=fit1$years, sim=sim, 
					F_fit=sdrep_out(fit1, s1, "log_F", bias.cor=F), F_fit_bc=sdrep_out(fit1, s1, "log_F", bias.cor=T), F_sim=fit1$env$data$F[,1], 
					relF_fit=calc_relF(fit1, s1, type="fit", bias.cor=F), relF_fit_bc=calc_relF(fit1, s1, type="fit", bias.cor=T), relF_sim=calc_relF(fit1, s1, type="sim"), 
					SSB_fit=sdrep_out(fit1, s1, "log_SSB", bias.cor=F), SSB_fit_bc=sdrep_out(fit1, s1, "log_SSB", bias.cor=T), SSB_sim=fit1$env$data$SSB, 
					relB_fit=calc_relB(fit1, s1, type="fit", bias.cor=F), relB_fit_bc=calc_relB(fit1, s1, type="fit", bias.cor=T), relB_sim=calc_relB(fit1, s1, type="sim"),  
					catch_fit=sdrep_out(fit1, s1, "log_pred_catch", bias.cor=F), catch_fit_bc=sdrep_out(fit1, s1, "log_pred_catch", bias.cor=T), catch_sim=fit1$env$data$pred_catch[,1]))
	dfnaa <- sdrep_out(fit1, s1, "log_NAA_rep", bias.cor=F)
	colnames(dfnaa) <- paste0("NAA",1:n.ages)
	dfnaa.bc <- sdrep_out(fit1, s1, "log_NAA_rep", bias.cor=T)
	colnames(dfnaa.bc) <- paste0("NAA",1:n.ages,"_bc")
	res <- cbind(df, dfnaa, dfnaa.bc)
	return(res)
}

# n.mods <- 4
n.sim <- 100
n.types <- 2
simdata <- readRDS("simdata_om1.rds")
n.years <- simdata[[1]][[1]][["n_years_model"]]
n.ages <- simdata[[1]][[1]][["n_ages"]]
rm("simdata")
sim.seeds <- readRDS("sim_seeds.rds")

options(warn=-1) # suppress warning messages

nested.list <- function(len) if(length(len) == 1) vector("list", len) else lapply(1:len[1], function(...) nested.list(len[-1]))
sdreps <- nested.list(c(n.types, n.sim)) # sdreps with analytical bias correction
# sdreps.bc <- nested.list(c(n.types, n.mods, n.sim))
reps <- nested.list(c(n.types, n.sim))
res.colnames <- c("om","em","type","year","sim","F_fit","F_fit_bc","F_sim","relF_fit","relF_fit_bc","relF_sim","SSB_fit","SSB_fit_bc","SSB_sim","relB_fit","relB_fit_bc","relB_sim","catch_fit","catch_fit_bc","catch_sim",paste0("NAA",1:n.ages),paste0("NAA",1:n.ages,"_bc"))
results <- rep(list(rep(list(matrix(NA, ncol = length(res.colnames), nrow = n.years)),n.sim)),n.types) # nested lists with preallocated matrices
# colnames(results) <- res.colnames
# for(m in 1:n.mods){
	simdata <- readRDS(file.path(simdata_dir,paste0("simdata_om",om,".rds")))
	for(i in 1:n.sim){
		# print(paste0("Model: ",m," Sim: ", i))
		print(paste0("OM: ",om," EM: ",em," Sim: ", i))
		set.seed(sim.seeds[i])

		# a) obs error, keep all parameters (incl NAA) as in fit model, simulate catch + index data
		input1 <- readRDS(file.path(res_dir,paste0("m",em,"_input.rds")))
		n.data <- length(input1$data)
		input1$data <- simdata[[i]][[1]]
		input1$data$n_NAA_sigma = length(input1$par$log_NAA_sigma)
		if(em > 2) input1$data$NAA_sigma_pointers = c(1,rep(2,input1$data$n_ages-1)) else input1$data$NAA_sigma_pointers = rep(1,input1$data$n_ages)
		ind.save <- c(1:n.data, match(c("F","SSB","pred_catch","log_FXSPR","FAA_tot","log_SSB_FXSPR"), names(input1$data)))
		input1$data <- input1$data[ind.save]
		fit1 <- tryCatch(fit_wham(input1, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE),
					error = function(e) conditionMessage(e))
		if(exists("err")) rm("err") # need to clean this up
		if(!'err' %in% names(fit1) & class(fit1) != "character"){
			reps[[1]][[i]] <- fit1$rep
			fit1$sdrep <- tryCatch(TMB::sdreport(fit1, bias.correct=TRUE), # also do bias correction
							error = function(e) conditionMessage(e))
			if(class(fit1$sdrep) == "sdreport"){
				s1 <- summary(fit1$sdrep)
				sdreps[[1]][[i]] <- s1
				results[[1]][[i]] <- tryCatch(calc_results(om=om, em=em, type=1, sim=i, fit1, s1),
					error = function(e) conditionMessage(e))
			} else {
				results[[1]][[i]] <- "Error: sdreport failed, no results to calculate"
				sdreps[[1]][[i]] <- fit1$sdrep # error message
			}
		} else {
			results[[1]][[i]] <- "Error: model did not converge, no results to calculate"
			if(class(fit1) != "character") reps[[1]][[i]] <- fit1$err # error message
			if(class(fit1) == "character") reps[[1]][[i]] <- fit1
			sdreps[[1]][[i]] <- "Error: model did not converge, sdreport not attempted"
		}

		# # turn off analytical bias correction and refit
		# input1$data$bias_correct_oe = 0
		# input1$data$bias_correct_pe = 0
		# fit1 <- fit_wham(input1, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)
		# if(exists("err")) rm("err") # need to clean this up
		# fit1$sdrep = TMB::sdreport(fit1, bias.correct = TRUE)

		# b) process + obs error, keep parameters except NAA as in fit model, simulate NAA (process error) and catch + index data (obs error)
		input2 <- readRDS(file.path(res_dir,paste0("m",em,"_input.rds")))
		input2$data <- simdata[[i]][[2]]
		input2$data$n_NAA_sigma = length(input2$par$log_NAA_sigma)
		if(em > 2) input2$data$NAA_sigma_pointers = c(1,rep(2,input2$data$n_ages-1)) else input2$data$NAA_sigma_pointers = rep(1,input2$data$n_ages)
		ind.save <- c(1:n.data, match(c("F","SSB","pred_catch","log_FXSPR","FAA_tot","log_SSB_FXSPR"), names(input2$data)))	
		input2$data <- input2$data[ind.save]
		fit2 <- tryCatch(fit_wham(input2, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE),
					error = function(e) conditionMessage(e))
		if(exists("err")) rm("err") # need to clean this up
		if(!'err' %in% names(fit2) & class(fit2) != "character"){
			reps[[2]][[i]] <- fit2$rep
			fit2$sdrep <- tryCatch(TMB::sdreport(fit2, bias.correct=TRUE), # also do bias correction
							error = function(e) conditionMessage(e))
			if(class(fit2$sdrep) == "sdreport"){
				s2 <- summary(fit2$sdrep)
				sdreps[[2]][[i]] <- s2
				results[[2]][[i]] <- tryCatch(calc_results(om=om, em=em, type=2, sim=i, fit2, s2),
					error = function(e) conditionMessage(e))
			} else {
				results[[2]][[i]] <- "Error: sdreport failed, no results to calculate"
				sdreps[[2]][[i]] <- fit2$sdrep # error message
			}
		} else {
			results[[2]][[i]] <- "Error: model did not converge, no results to calculate"
			if(class(fit2) != "character") reps[[2]][[i]] <- fit2$err # error message
			if(class(fit2) == "character") reps[[2]][[i]] <- fit2
			sdreps[[2]][[i]] <- "Error: model did not converge, sdreport not attempted"
		}

		# saveRDS(results, file=here("results","NAA",paste0("res_om",m,"_sim_",i,".rds")))
		rm(list=c("input1","input2","fit1","fit2","df","dfnaa","dfnaa.bc"))
	}
	rm(list=c("simdata"))
# }

# Collect results into one matrix
# flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
# res2 <- do.call(rbind, flatten.nested.list(results)) %>% as.data.frame
# saveRDS(res2, file=file.path(res_dir,"results.rds"))
saveRDS(results, file=paste0("results_om",om,"_em",em,".rds"))	
saveRDS(sdreps, file=paste0("sdreps_om",om,"_em",em,".rds"))	
saveRDS(reps, file=paste0("reps_om",om,"_em",em,".rds"))	

# un-suppress warnings
options(warn=0)
