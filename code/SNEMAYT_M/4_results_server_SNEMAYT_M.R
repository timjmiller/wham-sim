# Brian Stock
# June 15 2020
# Simulation test WHAM
# Step 4: Collect simulation fits
#   NEFSC server version 

# source("/home/bstock/Documents/ms/wham-sim/code/SNEMAYT_M/4_results_server_SNEMAYT_M.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)

# get results into data frame
res_dir <- here("results","SNEMAYT_M")
plots_dir <- here("plots","SNEMAYT_M")
res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
om2em2 <- res.files[4]
res.files <- res.files[-4] # om2/em2 in diff format
res.list <- lapply(res.files, readRDS)
flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
results <- do.call(rbind, flatten.nested.list(res.list)) %>% as.data.frame
results <- sapply(results, as.numeric)
types <- c("OE","OEPE")
results <- as.data.frame(results)
results <- rbind(results, readRDS(om2em2))
# results %>% group_by(om, em) %>% tally()

results$om <- factor(results$om, levels=1:3, labels=c("m1: none","m2: IID","m3: 2D AR1"))
results$em <- factor(results$em, levels=1:3, labels=c("m1: none","m2: IID","m3: 2D AR1"))
results$em.x <- fct_recode(results$em, m1="m1: none", m2="m2: IID", m3="m3: 2D AR1")
results$type <- factor(results$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new M re)"))

# Fig 1. SSB (sim fit) / SSB (sim data)
res.ssb <- results %>% group_by(om, em, type, sim) %>%
	mutate(SSB.rel = SSB_fit / SSB_sim, SSB.rel.bc = SSB_fit_bc / SSB_sim)
for(ty in 1:length(types)){
	df.plot <- filter(res.ssb, type==levels(res.ssb$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=SSB.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
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
		# scale_fill_discrete(name="Estimation model") +
		# stat_summary(fun = "mean", geom = "line", color = "red") +
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
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
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
		# scale_fill_discrete(name="Estimation model") +
		# stat_summary(fun = "mean", geom = "line", color = "red") +
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
res.F <- results %>% group_by(om, em, type, sim) %>%
	mutate(F.rel = F_fit / F_sim, F.rel.bc = F_fit_bc/F_sim)
for(ty in 1:length(types)){
	df.plot <- filter(res.F, type==levels(res.F$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=F.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
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
		# scale_fill_discrete(name="Estimation model") +
		# stat_summary(fun = "mean", geom = "line", color = "red") +
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
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
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
		# scale_fill_discrete(name="Estimation model") +
		# stat_summary(fun = "mean", geom = "line", color = "red") +
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
res.relB <- results %>% group_by(om, em, type, sim) %>%
	mutate(relB.rel = relB_fit / relB_sim, relB.rel.bc = relB_fit_bc/relB_sim)
for(ty in 1:length(types)){
	df.plot <- filter(res.relB, type==levels(res.relB$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=relB.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
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
		facet_grid(rows=vars(om), cols=vars(em)) +
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
res.relF <- results %>% group_by(om, em, type, sim) %>%
	mutate(relF.rel = relF_fit / relF_sim, relF.rel.bc = relF_fit_bc / relF_sim)
for(ty in 1:length(types)){
	df.plot <- filter(res.relF, type==levels(res.relF$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=relF.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
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
		facet_grid(rows=vars(om), cols=vars(em)) +
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
res.catch <- results %>% group_by(om, em, type, sim) %>%
	mutate(catch.rel = catch_fit / catch_sim, catch.rel.bc = catch_fit_bc / catch_sim)
for(ty in 1:length(types)){
	df.plot <- filter(res.catch, type==levels(res.catch$type)[ty])
	p <- ggplot(df.plot, aes(x=year, y=catch.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
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
		facet_grid(rows=vars(om), cols=vars(em)) +
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
# simdata <- lapply(1:3, function(x) readRDS(here("data","simdata","SNEMAYT_M",paste0("simdata_SNEMAYT_M_om",x,".rds"))))
# results <- results[complete.cases(results),]
# res.R <- results %>% group_by(om, em, type, sim) %>%
# 	mutate(R.sim = simdata[[unique(om)]][[unique(sim)]][[unique(type)]]$NAA[,1],
# 		   R.rel = NAA1 / R.sim,
# 		   R.rel.bc = NAA1_bc / R.sim)

mlabs = c("m1: none","m2: IID","m3: 2D AR1")
tylabs = c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)")
n.mods <- length(mlabs)
simdata <- lapply(1:n.mods, function(x) readRDS(here("data","simdata","SNEMAYT_M",paste0("simdata_SNEMAYT_M_om",x,".rds"))))
results$R.sim = NA
for(om in 1:n.mods){
	for(em in 1:n.mods){
		for(i in 1:100){
			for(ty in 1:2){
				res.ind <- which(results$om == mlabs[om] & results$em == mlabs[em] & results$sim == i & results$ty == tylabs[ty])
				results$R.sim[res.ind] <- simdata[[om]][[i]][[ty]]$NAA[,1]
			}
		}
	}
}
results$R.rel <- results$NAA1 / results$R.sim
results$R.rel.bc <- results$NAA1_bc / results$R.sim
res.R <- results %>% group_by(om, em, type, sim) 
plot.bc = FALSE

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
		facet_grid(rows=vars(om), cols=vars(em)) +
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
	if(plot.bc){
	p <- ggplot(df.plot, aes(x=year, y=R.rel.bc)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
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
}

