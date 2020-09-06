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

# source("/home/bstock/Documents/ms/wham-sim/code/bias_correct_oepe/SNEMAYT_NAA_oepe/4_results_SNEMAYT_NAA.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)

# get results into data frame
res_dir <- here("results","bias_correct_oepe","SNEMAYT_NAA_oepe")
plots_dir <- here("plots","bias_correct_oepe","SNEMAYT_NAA")
res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
res.list <- lapply(res.files, readRDS)
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
results$om <- factor(results$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
results$em <- factor(results$em, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
results$type <- factor(results$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
results$em.x <- fct_recode(results$em, m1="m1: SCAA (iid)", m2="m2: SCAA (AR1_y)", m3="m3: NAA (iid)", m4="m4: NAA (2D AR1)")

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
	png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
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

	png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
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
	png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
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

	png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
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
	png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
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

	png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
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
	png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
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

	png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
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
	png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
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

	png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
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
simdata <- lapply(1:4, function(x) readRDS(here("data","simdata","bias_correct_oepe","SNEMAYT_NAA_oepe",paste0("simdata_om",x,".rds"))))
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
	png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
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

	png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
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

