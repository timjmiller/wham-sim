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

# source("/home/bstock/Documents/ms/wham-sim/code/4_results_NAA_server.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)

# get results into data frame
res_dir <- here("results","SNEMAYT_NAA")
plots_dir <- here("plots","SNEMAYT_NAA")
res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
res.list <- lapply(res.files, readRDS)
flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
results <- do.call(rbind, flatten.nested.list(res.list)) %>% as.data.frame
results <- sapply(results, as.numeric)
types <- c("OE","OEPE")
results <- as.data.frame(results)

# Fig 1. SSB (sim fit) / SSB (sim data)
res.ssb <- results %>% group_by(om, em, type, sim) %>%
	mutate(SSB.rel = SSB_fit / SSB_sim, SSB.rel.bc = SSB_fit_bc / SSB_sim)
	# mutate(SSB.fit = mods[[unique(om)]]$rep$SSB,
	# 	   SSB.rel = SSB / SSB.fit)
res.ssb$om <- factor(res.ssb$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.ssb$em <- factor(res.ssb$em, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.ssb$type <- factor(res.ssb$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
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

	png(file.path(plots_dir,paste0("1_ssb_",types[ty],".png")), width=7, height=7, units='in',res=300)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	# boxplots (collapse time series)
	df.plot$em.x <- fct_recode(df.plot$em, m1="m1: SCAA (iid)", m2="m2: SCAA (AR1_y)", m3="m3: NAA (iid)", m4="m4: NAA (2D AR1)")
	png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],".png")), width=8, height=3, units='in',res=300)
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

	png(file.path(plots_dir,paste0("1_ssb_",types[ty],"_bc.png")), width=7, height=7, units='in',res=300)
	print(p)
	grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
	grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
	dev.off()

	png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=300)
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
	mutate(F.rel = F_fit / F_sim)
	# mutate(F.fit = mods[[unique(om)]]$rep$F[,1],
	# 	   F.rel = F / F.fit)
res.F$om <- factor(res.F$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.F$em <- factor(res.F$em, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.F$type <- factor(res.F$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
for(ty in 1:length(types)){
	df.plot <- filter(res.F, type==ty)
	png(file.path(plots_dir,paste0("2_F_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(ggplot(df.plot, aes(x=year, y=F.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
		theme_bw())
	dev.off()
}

# Fig 3. relSSB (sim data) / relSSB (true data)
res.relB <- results %>% group_by(om, em, type, sim) %>%
	mutate(relB.rel = relB_fit / relB_sim)
	# mutate(relB.fit = calc_relB(mods[[unique(om)]]),
	# 	   relB.rel = relB / relB.fit)
res.relB$om <- factor(res.relB$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.relB$em <- factor(res.relB$em, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.relB$type <- factor(res.relB$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
for(ty in 1:length(types)){
	df.plot <- filter(res.relB, type==ty)
	png(file.path(plots_dir,paste0("3_relB_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(ggplot(df.plot, aes(x=year, y=relB.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
		theme_bw())
	dev.off()
}

# Fig 4. relF (sim data) / relF (true data)
res.relF <- results %>% group_by(om, em, type, sim) %>%
	mutate(relF.rel = relF_fit / relF_sim)
	# mutate(relF.fit = calc_relF(mods[[unique(om)]]),
	# 	   relF.rel = relF / relF.fit)
res.relF$om <- factor(res.relF$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.relF$em <- factor(res.relF$em, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.relF$type <- factor(res.relF$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
for(ty in 1:length(types)){
	df.plot <- filter(res.relF, type==ty)
	png(file.path(plots_dir,paste0("4_relF_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(ggplot(df.plot, aes(x=year, y=relF.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
		theme_bw())
	dev.off()
}

# Fig 5. pred_catch (sim data) / pred_catch (true data)
res.catch <- results %>% group_by(om, em, type, sim) %>%
	mutate(catch.rel = catch_fit / catch_sim)
	# mutate(catch.fit = mods[[unique(om)]]$rep$pred_catch[,1],
	# 	   catch.rel = pred_catch / catch.fit)
res.catch$om <- factor(res.catch$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.catch$em <- factor(res.catch$em, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.catch$type <- factor(res.catch$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
for(ty in 1:length(types)){
	df.plot <- filter(res.catch, type==ty)
	png(file.path(plots_dir,paste0("5_catch_",types[ty],".png")), width=7, height=7, units='in',res=100)
	print(ggplot(df.plot, aes(x=year, y=catch.rel)) +
	    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
	    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
		# geom_line(aes(group=sim),color='grey', alpha=.2) +
		stat_summary(fun = "median", geom = "line", color = "red") +
		coord_cartesian(ylim=c(0,2)) +
		xlab("Year") +
		ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
		geom_hline(yintercept = 1, linetype=2, color='black') +
		facet_grid(rows=vars(om), cols=vars(em)) +
		theme_bw())
	dev.off()	
}

# Fig 6. Recruitment (sim data) / Recruitment (true data)
simdata <- lapply(1:4, function(x) readRDS(here("data","simdata","SNEMAYT","NAA",paste0("simdata_om",x,".rds"))))
res.R <- results %>% group_by(om, em, type, sim) %>%
	mutate(R.sim = simdata[[unique(om)]][[unique(sim)]][[unique(type)]]$NAA[,1],
		   R.rel = NAA1 / R.sim)
res.R$om <- factor(res.R$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.R$type <- factor(res.R$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
png(here("plots","NAA","6_R.png"), width=7, height=7, units='in',res=100)
print(ggplot(res.R, aes(x=year, y=R.rel)) +
    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
	# geom_line(aes(group=sim),color='grey', alpha=.2) +
	stat_summary(fun = "median", geom = "line", color = "red") +
	coord_cartesian(ylim=c(0,2)) +
	xlab("Year") +
	ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
	geom_hline(yintercept = 1, linetype=2, color='black') +
	facet_grid(rows=vars(om), cols=vars(type)) +
	theme_bw())
dev.off()


