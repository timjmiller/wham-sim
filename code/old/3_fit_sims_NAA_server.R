# Brian Stock
# June 15 2020
# Simulation test WHAM
# Step 3: Fit OMs to simulated datasets
#   NEFSC server version

# the following files should be copied to working directory:
#   simdata_omXX.rds (XX from 1-4)
#   mXX_input.rds (XX from 1-4)
#   sim_seeds.rds

# will write 3 files to working directory:
#   results_omXX_emYY.rds
#   sdreps_omXX_emYY.rds
#   reps_omXX_emYY.rds

args = commandArgs(trailingOnly=TRUE)
om = args[1] # operating model
em = args[2] # estimating model

# -----------------------------------------------------------------------
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
# library(here)
# library(tidyverse)

res_dir <- simdata_dir <- getwd()
# res_dir <- here("results","OE_PE","NAA")
# simdata_dir <- here("data","simdata","SNEMAYT","NAA")

# functions to calculate F/F40 and B/B40
calc_relF <- function(mod, sdrep, type="fit", bias.cor=TRUE){ # sdrep is summary(sdreport)
	if(type == "fit"){
		ind.FXSPR <- which(rownames(sdrep) == "log_FXSPR")
		if(bias.cor) F.t <- sdrep[ind.FXSPR,3] else F.t <- sdrep[ind.FXSPR,1]
		ind.faa <- which(rownames(sdrep) == "log_FAA_tot")
		n.yrs <- length(mod$years_full)
		n.ages <- mod$env$data$n_ages
  		if(bias.cor) faa <- matrix(sdrep[ind.faa,3], n.yrs, n.ages) else faa <- matrix(sdrep[ind.faa,1], n.yrs, n.ages)
		# F.t <- mod$rep$log_FXSPR
		# faa <- log(mod$rep$FAA_tot)		
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
		# SSB.t <- exp(mod$rep$log_SSB_FXSPR)
		# ssb <- mod$rep$SSB
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

# Assumes you open R in project directory
# Load fit models from step 1
# mod.list <- here("results","NAA",paste0("m",1:4,".rds"))
# mods <- lapply(mod.list, readRDS)
# input.list <- here("results","NAA",paste0("m",1:4,"_input.rds"))
# inputs <- lapply(input.list, readRDS)

# n.mods <- 4
n.sim <- 100
n.types <- 2
simdata <- readRDS("simdata_om1.rds")
n.years <- simdata[[1]][[1]][["n_years_model"]]
rm("simdata")
sim.seeds <- readRDS("sim_seeds.rds")

options(warn=-1) # suppress warning messages

nested.list <- function(len) if(length(len) == 1) vector("list", len) else lapply(1:len[1], function(...) nested.list(len[-1]))
sdreps <- nested.list(c(n.types, n.sim)) # sdreps with analytical bias correction
# sdreps.bc <- nested.list(c(n.types, n.mods, n.sim))
reps <- nested.list(c(n.types, n.sim))
res.colnames <- c("om","em","type","year","sim","F_fit","F_fit_bc","F_sim","relF_fit","relF_fit_bc","relF_sim","SSB_fit","SSB_fit_bc","SSB_sim","relB_fit","relB_fit_bc","relB_sim","catch_fit","catch_fit_bc","catch_sim",paste0("NAA",1:6),paste0("NAA",1:6,"_bc"))
results <- rep(list(rep(list(matrix(NA, ncol = length(res.colnames), nrow = n.years)),n.sim)),n.types) # nested lists with preallocated matrices
# colnames(results) <- res.colnames
# for(m in 1:n.mods){
	simdata <- readRDS(file.path(simdata_dir,paste0("simdata_om",om,".rds")))
	for(i in 1:n.sim){
		# print(paste0("Model: ",m," Sim: ", i))
		print(paste0("Sim: ", i))
		set.seed(sim.seeds[i])

		# a) obs error, keep all parameters (incl NAA) as in fit model, simulate catch + index data
		input1 <- readRDS(file.path(res_dir,paste0("m",em,"_input.rds")))
		n.data <- length(input1$data)
		input1$data <- simdata[[i]][[1]]
		input1$data$n_NAA_sigma = length(input1$par$log_NAA_sigma)
		if(em > 2) input1$data$NAA_sigma_pointers = c(1,rep(2,input1$data$n_ages-1)) else input1$data$NAA_sigma_pointers = rep(1,input1$data$n_ages)
		ind.save <- c(1:n.data, match(c("F","SSB","pred_catch","log_FXSPR","FAA_tot","log_SSB_FXSPR"), names(input1$data)))
		input1$data <- input1$data[ind.save]
		fit1 <- fit_wham(input1, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)
		if(exists("err")) rm("err") # need to clean this up
		fit1$sdrep = TMB::sdreport(fit1, bias.correct=TRUE) # also do bias correction
		s1 <- summary(fit1$sdrep)

		# # turn off analytical bias correction and refit
		# input1$data$bias_correct_oe = 0
		# input1$data$bias_correct_pe = 0
		# fit1 <- fit_wham(input1, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)
		# if(exists("err")) rm("err") # need to clean this up
		# fit1$sdrep = TMB::sdreport(fit1, bias.correct = TRUE)

		df <- as.matrix(data.frame(om=om, em=em, type=1, year=fit1$years, sim=i, 
			F_fit=sdrep_out(fit1, s1, "log_F", bias.cor=F), F_fit_bc=sdrep_out(fit1, s1, "log_F", bias.cor=T), F_sim=fit1$env$data$F[,1], 
			relF_fit=calc_relF(fit1, s1, type="fit", bias.cor=F), relF_fit_bc=calc_relF(fit1, s1, type="fit", bias.cor=T), relF_sim=calc_relF(fit1, s1, type="sim"), 
			SSB_fit=sdrep_out(fit1, s1, "log_SSB", bias.cor=F), SSB_fit_bc=sdrep_out(fit1, s1, "log_SSB", bias.cor=T), SSB_sim=fit1$env$data$SSB, 
			relB_fit=calc_relB(fit1, s1, type="fit", bias.cor=F), relB_fit_bc=calc_relB(fit1, s1, type="fit", bias.cor=T), relB_sim=calc_relB(fit1, s1, type="sim"),  
			catch_fit=sdrep_out(fit1, s1, "log_pred_catch", bias.cor=F), catch_fit_bc=sdrep_out(fit1, s1, "log_pred_catch", bias.cor=T), catch_sim=fit1$env$data$pred_catch[,1]))
		dfnaa <- sdrep_out(fit1, s1, "log_NAA_rep", bias.cor=F); colnames(dfnaa) <- paste0("NAA",1:6)
		dfnaa.bc <- sdrep_out(fit1, s1, "log_NAA_rep", bias.cor=T); colnames(dfnaa.bc) <- paste0("NAA",1:6,"_bc")
		results[[1]][[i]] <- cbind(df, dfnaa, dfnaa.bc)
		sdreps[[1]][[i]] <- s1
		reps[[1]][[i]] <- fit1$rep

		# b) process + obs error, keep parameters except NAA as in fit model, simulate NAA (process error) and catch + index data (obs error)
		input2 <- readRDS(file.path(res_dir,paste0("m",em,"_input.rds")))
		input2$data <- simdata[[i]][[2]]
		input2$data$n_NAA_sigma = length(input2$par$log_NAA_sigma)
		if(em > 2) input2$data$NAA_sigma_pointers = c(1,rep(2,input2$data$n_ages-1)) else input2$data$NAA_sigma_pointers = rep(1,input2$data$n_ages)
		ind.save <- c(1:n.data, match(c("F","SSB","pred_catch","log_FXSPR","FAA_tot","log_SSB_FXSPR"), names(input2$data)))	
		input2$data <- input2$data[ind.save]
		fit2 <- fit_wham(input2, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)		
		if(exists("err")) rm("err") # need to clean this up
		fit2$sdrep = TMB::sdreport(fit2, bias.correct=TRUE) # also do bias correction
		s2 <- summary(fit2$sdrep)

		df <- as.matrix(data.frame(om=om, em=em, type=2, year=fit2$years, sim=i, 
			F_fit=sdrep_out(fit2, s2, "log_F", bias.cor=F), F_fit_bc=sdrep_out(fit2, s2, "log_F", bias.cor=T), F_sim=fit2$env$data$F[,1], 
			relF_fit=calc_relF(fit2, s2, type="fit", bias.cor=F), relF_fit_bc=calc_relF(fit2, s2, type="fit", bias.cor=T), relF_sim=calc_relF(fit2, s2, type="sim"), 
			SSB_fit=sdrep_out(fit2, s2, "log_SSB", bias.cor=F), SSB_fit_bc=sdrep_out(fit2, s2, "log_SSB", bias.cor=T), SSB_sim=fit2$env$data$SSB, 
			relB_fit=calc_relB(fit2, s2, type="fit", bias.cor=F), relB_fit_bc=calc_relB(fit2, s2, type="fit", bias.cor=T), relB_sim=calc_relB(fit2, s2, type="sim"),  
			catch_fit=sdrep_out(fit2, s2, "log_pred_catch", bias.cor=F), catch_fit_bc=sdrep_out(fit2, s2, "log_pred_catch", bias.cor=T), catch_sim=fit2$env$data$pred_catch[,1]))
		dfnaa <- sdrep_out(fit2, s2, "log_NAA_rep", bias.cor=F); colnames(dfnaa) <- paste0("NAA",1:6)
		dfnaa.bc <- sdrep_out(fit2, s2, "log_NAA_rep", bias.cor=T); colnames(dfnaa.bc) <- paste0("NAA",1:6,"_bc")
		results[[2]][[i]] <- cbind(df, dfnaa, dfnaa.bc)
		sdreps[[2]][[i]] <- s2
		reps[[2]][[i]] <- fit2$rep

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
