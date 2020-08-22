# Brian Stock
# June 15 2020
# Simulation test WHAM
# Step 4: Collect simulation fits
#   NEFSC server version 
#   separate results files for each om/em cross test, for XX/YY in 1-4
#     results_omXX_emYY.rds
#     sdreps_omXX_emYY.rds
#     reps_omXX_emYY.rds
#   assume they are copied to /home/bstock/Documents/ms/wham-sim/results/SNEMAYT/NAA

# source("/home/bstock/Documents/ms/wham-sim/code/NScod_M/4_results_NScod_M.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)

# functions to calculate F/F40 and B/B40
calc_relF <- function(simdata, sdrep, type="fit", bias.cor=TRUE){ # sdrep is summary(sdreport)
	if(type == "fit"){
		ind.FXSPR <- which(rownames(sdrep) == "log_FXSPR")
		if(bias.cor) F.t <- sdrep[ind.FXSPR,3] else F.t <- sdrep[ind.FXSPR,1]
		ind.faa <- which(rownames(sdrep) == "log_FAA_tot")
		n.yrs <- simdata$n_years_model
		n.ages <- simdata$n_ages
  		if(bias.cor) faa <- matrix(sdrep[ind.faa,3], n.yrs, n.ages) else faa <- matrix(sdrep[ind.faa,1], n.yrs, n.ages)		
	}
	if(type == "sim"){
		# F.t <- mod$env$data$log_FXSPR
		# faa <- log(mod$env$data$FAA_tot)			
		F.t <- simdata$log_FXSPR
		faa <- log(simdata$FAA_tot)
	}
  	age.full.f <- apply(faa,1, function(x) max(which(x == max(x))))
  	full.f <- faa[cbind(seq_along(age.full.f),age.full.f)]
	rel.f <- exp(full.f - F.t)
	return(rel.f)
}
calc_relB <- function(simdata, sdrep, type="fit", bias.cor=TRUE){ # sdrep is summary(sdreport)
	if(type == "fit"){
		ind.SSB.FXSPR <- which(rownames(sdrep) == "log_SSB_FXSPR")
		if(bias.cor) SSB.t <- exp(sdrep[ind.SSB.FXSPR,3]) else SSB.t <- exp(sdrep[ind.SSB.FXSPR,1])
		ind.ssb <- which(rownames(sdrep) == "log_SSB")
  		if(bias.cor) ssb <- exp(sdrep[ind.ssb,3]) else ssb <- exp(sdrep[ind.ssb,1])
	}
	if(type == "sim"){
		SSB.t <- exp(simdata$log_SSB_FXSPR)
		ssb <- simdata$SSB
	}
	rel.ssb <- ssb / SSB.t	
	return(rel.ssb)
}
sdrep_out <- function(simdata, sdrep, var, bias.cor=TRUE){
	if(var == "log_NAA_rep"){
		ind <- which(rownames(sdrep) == var)
		# n.yrs <- length(mod$years_full)
		# n.ages <- mod$env$data$n_ages
		n.yrs <- simdata$n_years_model
		n.ages <- simdata$n_ages
		if(bias.cor) x <- matrix(sdrep[ind,3], n.yrs, n.ages) else x <- matrix(sdrep[ind,1], n.yrs, n.ages)
	} else {
		ind <- which(rownames(sdrep) == var)
		if(bias.cor) x <- sdrep[ind,3] else x <- sdrep[ind,1]
	}
	return(exp(x)) # assumed on log scale
}
calc_results <- function(om, em, type, sim, simdata, s1){
	n.ages <- simdata$n_ages
	df <- as.matrix(data.frame(om=om, em=em, type=type, year=(1:simdata$n_years_model + simdata$year1_model - 1), sim=sim, 
					F_fit=sdrep_out(simdata, s1, "log_F", bias.cor=F), F_fit_bc=sdrep_out(simdata, s1, "log_F", bias.cor=T), F_sim=simdata$F[,1], 
					relF_fit=calc_relF(simdata, s1, type="fit", bias.cor=F), relF_fit_bc=calc_relF(simdata, s1, type="fit", bias.cor=T), relF_sim=calc_relF(simdata, s1, type="sim"), 
					SSB_fit=sdrep_out(simdata, s1, "log_SSB", bias.cor=F), SSB_fit_bc=sdrep_out(simdata, s1, "log_SSB", bias.cor=T), SSB_sim=simdata$SSB, 
					relB_fit=calc_relB(simdata, s1, type="fit", bias.cor=F), relB_fit_bc=calc_relB(simdata, s1, type="fit", bias.cor=T), relB_sim=calc_relB(simdata, s1, type="sim"),  
					catch_fit=sdrep_out(simdata, s1, "log_pred_catch", bias.cor=F), catch_fit_bc=sdrep_out(simdata, s1, "log_pred_catch", bias.cor=T), catch_sim=simdata$pred_catch[,1]))
	dfnaa <- sdrep_out(simdata, s1, "log_NAA_rep", bias.cor=F)
	colnames(dfnaa) <- paste0("NAA",1:n.ages)
	dfnaa.bc <- sdrep_out(simdata, s1, "log_NAA_rep", bias.cor=T)
	colnames(dfnaa.bc) <- paste0("NAA",1:n.ages,"_bc")
	res <- cbind(df, dfnaa, dfnaa.bc)
	return(res)
}

n.sim <- 100
n.types <- 2
res_dir <- here("results","NScod_M")
simdata_dir <- here("data","simdata","NScod_M")
plots_dir <- here("plots","NScod_M")
simdata <- readRDS(file.path(simdata_dir,paste0("simdata_om1.rds")))
n.years <- simdata[[1]][[1]][["n_years_model"]]
n.ages <- simdata[[1]][[1]][["n_ages"]]
rm("simdata")
res.colnames <- c("om","em","type","year","sim","F_fit","F_fit_bc","F_sim","relF_fit","relF_fit_bc","relF_sim","SSB_fit","SSB_fit_bc","SSB_sim","relB_fit","relB_fit_bc","relB_sim","catch_fit","catch_fit_bc","catch_sim",paste0("NAA",1:n.ages),paste0("NAA",1:n.ages,"_bc"))
res.list <- vector("list", 2*2)
il <- 1
for(om in 1:2){
	simdata <- readRDS(file.path(simdata_dir,paste0("simdata_om",om,".rds")))
	for(em in 1:2){
		results <- rep(list(rep(list(matrix(NA, ncol = length(res.colnames), nrow = n.years)),n.sim)),n.types) # nested lists with preallocated matrices
		reps <- readRDS(file.path(res_dir,paste0("reps_om",om,"_em",em,".rds")))
		sdreps <- readRDS(file.path(res_dir,paste0("sdreps_om",om,"_em",em,".rds")))
		for(i in 1:n.sim){
			for(ty in 1:n.types){
				results[[ty]][[i]] <- tryCatch(calc_results(om=om, em=em, type=ty, sim=i, simdata=simdata[[i]][[ty]], s1=sdreps[[ty]][[i]]),
					error = function(e) conditionMessage(e))
			}
		}
		res.list[[il]] <- results
		il <- il + 1
	}
}

# # get results into data frame
# res_dir <- here("results","NScod_M")
# plots_dir <- here("plots","NScod_M")
# res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
# res.list <- lapply(res.files, readRDS)
# for(i in 1:length(res.list)){
# 	for(j in 1:length(res.list[[i]])){
# 		for(k in 1:length(res.list[[i]][[j]])){
# 			if(class(res.list[[i]][[j]][[k]])=="character"){
# 				res.list[[i]][[j]][[k]] <- res.list[[1]][[1]][[1]][FALSE,]
# 			}
# 		}
# 	}
# }
# flatten <- function(x) {
#   if (!inherits(x, "list")) return(list(x))
#   else return(unlist(c(lapply(x, flatten)), recursive = FALSE))
# }
# results <- do.call(rbind, flatten(res.list))

flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
results <- do.call(rbind, flatten.nested.list(res.list)) %>% as.data.frame
results <- sapply(results, as.numeric)
results <- as.data.frame(results)

types <- c("OE","OEPE")
results$om <- factor(results$om, levels=1:2, labels=c("m1: none","m2: IID"))
results$em <- factor(results$em, levels=1:2, labels=c("m1: none","m2: IID"))
results$type <- factor(results$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
results$em.x <- fct_recode(results$em, m1="m1: none", m2="m2: IID")

# calculate relative error
results$SSB.rel = results$SSB_fit / results$SSB_sim
results$SSB.rel.bc = results$SSB_fit_bc / results$SSB_sim
results$F.rel = results$F_fit / results$F_sim
results$F.rel.bc = results$F_fit_bc / results$F_sim
results$relB.rel = results$relB_fit / results$relB_sim
results$relB.rel.bc = results$relB_fit_bc / results$relB_sim
results$relF.rel = results$relF_fit / results$relF_sim
results$relF.rel.bc = results$relF_fit_bc / results$relF_sim
results$catch.rel = results$catch_fit / results$catch_sim
results$catch.rel.bc = results$catch_fit_bc / results$catch_sim

# Fig 1. SSB (sim fit) / SSB (sim data)
for(ty in 1:length(types)){
	df.plot <- filter(results, type==levels(results$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=SSB.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
	png(file.path(plots_dir,paste0("1_ssb_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	# boxplots (collapse time series)
	png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],".png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=SSB.rel)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()

	# with TMB bias correction
	p <- ggplot(df.plot, aes(x=year, y=SSB.rel.bc)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

	png(file.path(plots_dir,paste0("1_ssb_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],"_bc.png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=SSB.rel.bc)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()	
}

# Fig 2. F (sim fit) / F (sim data)
for(ty in 1:length(types)){
	df.plot <- filter(results, type==levels(results$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=F.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

	png(file.path(plots_dir,paste0("2_F_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	# boxplots (collapse time series)
	png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],".png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=F.rel)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()

	# with TMB bias correction
	p <- ggplot(df.plot, aes(x=year, y=F.rel.bc)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

	png(file.path(plots_dir,paste0("2_F_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],"_bc.png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=F.rel.bc)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()		
}

# Fig 3. relSSB (sim data) / relSSB (true data)
# mistake in relB calculations... need to exp(log_SSB)
# fix_relB <- function(s1, bias.cor=TRUE){ # sdrep is summary(sdreport)
# 	ind.SSB.FXSPR <- which(rownames(s1) == "log_SSB_FXSPR")
# 	if(bias.cor) SSB.t <- exp(s1[ind.SSB.FXSPR,3]) else SSB.t <- exp(s1[ind.SSB.FXSPR,1])
# 	ind.ssb <- which(rownames(s1) == "log_SSB")
# 	if(bias.cor) ssb <- exp(s1[ind.ssb,3]) else ssb <- exp(s1[ind.ssb,1])
# 	rel.ssb <- ssb / SSB.t	
# 	return(rel.ssb)
# }
# for(om in 1:4){
# 	for(em in 1:4){
# 		res <- readRDS(file=file.path(res_dir,paste0("results_om",om,"_em",em,".rds")))			
# 		sdrep <- readRDS(file=file.path(res_dir,paste0("sdreps_om",om,"_em",em,".rds")))	
# 		for(ty in 1:2){
# 			for(i in 1:100){
# 				s1 <- sdrep[[ty]][[i]]
# 				if(class(s1) != "character"){
# 					tmp <- as.data.frame(sapply(as.data.frame(res[[ty]][[i]]), as.numeric))
# 					tmp$relB_fit = fix_relB(s1, bias.cor=F)
# 					tmp$relB_fit_bc = fix_relB(s1, bias.cor=T)
# 					res[[ty]][[i]] <- tmp					
# 				}
# 			}
# 		}
# 		saveRDS(res, file=file.path(res_dir,paste0("results_om",om,"_em",em,".rds")))	
# 	}
# }
for(ty in 1:length(types)){
	df.plot <- filter(results, type==levels(results$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=relB.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

	png(file.path(plots_dir,paste0("3_relB_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	# boxplots (collapse time series)
	png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],".png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=relB.rel)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		# scale_fill_discrete(name="Estimation model") +
		# stat_summary(fun = "mean", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()

	# with TMB bias correction
	p <- ggplot(df.plot, aes(x=year, y=relB.rel.bc)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

	png(file.path(plots_dir,paste0("3_relB_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],"_bc.png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=relB.rel.bc)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		# scale_fill_discrete(name="Estimation model") +
		# stat_summary(fun = "mean", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()		
}

# Fig 4. relF (sim data) / relF (true data)
for(ty in 1:length(types)){
	df.plot <- filter(results, type==levels(results$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=relF.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
	png(file.path(plots_dir,paste0("4_relF_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	# boxplots (collapse time series)
	png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],".png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=relF.rel)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()

	# with TMB bias correction
	p <- ggplot(df.plot, aes(x=year, y=relF.rel.bc)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
	png(file.path(plots_dir,paste0("4_relF_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],"_bc.png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=relF.rel.bc)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()		
}

# Fig 5. pred_catch (sim data) / pred_catch (true data)
for(ty in 1:length(types)){
	df.plot <- filter(results, type==levels(results$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=catch.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
	png(file.path(plots_dir,paste0("5_catch_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	# boxplots (collapse time series)
	png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],".png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=catch.rel)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()

	# with TMB bias correction
	p <- ggplot(df.plot, aes(x=year, y=catch.rel.bc)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
	png(file.path(plots_dir,paste0("5_catch_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],"_bc.png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=catch.rel.bc)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()		
}

# Fig 6. Recruitment (sim data) / Recruitment (true data)
simdata <- lapply(1:2, function(x) readRDS(here("data","simdata","NScod_M",paste0("simdata_om",x,".rds"))))
results <- results[complete.cases(results),]
res.R <- results %>% group_by(om, em, type, sim) %>%
	mutate(R.sim = simdata[[unique(om)]][[unique(sim)]][[unique(type)]]$NAA[,1],
		   R.rel = NAA1 / R.sim,
		   R.rel.bc = NAA1_bc / R.sim)
for(ty in 1:length(types)){
	df.plot <- filter(res.R, type==levels(res.R$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=R.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
	png(file.path(plots_dir,paste0("6_R_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	# boxplots (collapse time series)
	png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],".png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=R.rel)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()

	# with TMB bias correction
	p <- ggplot(df.plot, aes(x=year, y=R.rel.bc)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(em), cols=vars(om)) +
		theme_bw() +
		theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
	png(file.path(plots_dir,paste0("6_R_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],"_bc.png")), width=5, height=3, units='in',res=100)
	print(ggplot(df.plot, aes(x=em.x, y=R.rel.bc)) +
		geom_boxplot(aes(fill=em), outlier.shape = NA) +
		scale_fill_jco(name="Estimation model") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Estimation model") +
		ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
		labs(title="Operating model") +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_wrap(vars(om), nrow=1) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)))
	dev.off()		
}

