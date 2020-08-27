# Brian Stock
# June 15 2020
# Simulation test WHAM

# source("/home/bstock/Documents/ms/wham-sim/code/4_results_NAA.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "om_mode")
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)

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

sdreps <- readRDS(here("results","NAA","sdreps.rds"))
results <- readRDS(here("results","NAA","results.rds"))
mod.list <- here("results","NAA",paste0("m",1:4,".rds"))
mods <- lapply(mod.list, readRDS)
input.list <- here("results","NAA",paste0("m",1:4,"_input.rds"))
inputs <- lapply(input.list, readRDS)
n.mods <- length(mods)

# Fig 1. SSB (sim fit) / SSB (sim data)
res.ssb <- results %>% group_by(om, type, sim) %>%
	mutate(SSB.rel = SSB_fit / SSB_sim)
	# mutate(SSB.fit = mods[[unique(om)]]$rep$SSB,
	# 	   SSB.rel = SSB / SSB.fit)
res.ssb$om <- factor(res.ssb$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.ssb$type <- factor(res.ssb$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
png(here("plots","NAA","1_ssb.png"), width=7, height=7, units='in',res=100)
print(ggplot(res.ssb, aes(x=year, y=SSB.rel)) +
    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
	# geom_line(aes(group=sim),color='grey', alpha=.2) +
	stat_summary(fun = "median", geom = "line", color = "red") +
	coord_cartesian(ylim=c(0,2)) +
	xlab("Year") +
	ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
	geom_hline(yintercept = 1, linetype=2, color='black') +
	facet_grid(rows=vars(om), cols=vars(type)) +
	theme_bw())
dev.off()

# Fig 2. F (sim fit) / F (sim data)
res.F <- results %>% group_by(om, type, sim) %>%
	mutate(F.rel = F_fit / F_sim)
	# mutate(F.fit = mods[[unique(om)]]$rep$F[,1],
	# 	   F.rel = F / F.fit)
res.F$om <- factor(res.F$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.F$type <- factor(res.F$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
png(here("plots","NAA","2_F.png"), width=7, height=7, units='in',res=100)
print(ggplot(res.F, aes(x=year, y=F.rel)) +
    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
	# geom_line(aes(group=sim),color='grey', alpha=.2) +
	stat_summary(fun = "median", geom = "line", color = "red") +
	coord_cartesian(ylim=c(0,2)) +
	xlab("Year") +
	ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
	geom_hline(yintercept = 1, linetype=2, color='black') +
	facet_grid(rows=vars(om), cols=vars(type)) +
	theme_bw())
dev.off()

# Fig 3. relSSB (sim data) / relSSB (true data)
res.relB <- results %>% group_by(om, type, sim) %>%
	mutate(relB.rel = relB_fit / relB_sim)
	# mutate(relB.fit = calc_relB(mods[[unique(om)]]),
	# 	   relB.rel = relB / relB.fit)
res.relB$om <- factor(res.relB$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.relB$type <- factor(res.relB$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
png(here("plots","NAA","3_relB.png"), width=7, height=7, units='in',res=100)
print(ggplot(res.relB, aes(x=year, y=relB.rel)) +
    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
	# geom_line(aes(group=sim),color='grey', alpha=.2) +
	stat_summary(fun = "median", geom = "line", color = "red") +
	coord_cartesian(ylim=c(0,2)) +
	xlab("Year") +
	ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
	geom_hline(yintercept = 1, linetype=2, color='black') +
	facet_grid(rows=vars(om), cols=vars(type)) +
	theme_bw())
dev.off()

# Fig 4. relF (sim data) / relF (true data)
res.relF <- results %>% group_by(om, type, sim) %>%
	mutate(relF.rel = relF_fit / relF_sim)
	# mutate(relF.fit = calc_relF(mods[[unique(om)]]),
	# 	   relF.rel = relF / relF.fit)
res.relF$om <- factor(res.relF$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.relF$type <- factor(res.relF$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
png(here("plots","NAA","4_relF.png"), width=7, height=7, units='in',res=100)
print(ggplot(res.relF, aes(x=year, y=relF.rel)) +
    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
	# geom_line(aes(group=sim),color='grey', alpha=.2) +
	stat_summary(fun = "median", geom = "line", color = "red") +
	coord_cartesian(ylim=c(0,2)) +
	xlab("Year") +
	ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
	geom_hline(yintercept = 1, linetype=2, color='black') +
	facet_grid(rows=vars(om), cols=vars(type)) +
	theme_bw())
dev.off()

# Fig 5. pred_catch (sim data) / pred_catch (true data)
res.catch <- results %>% group_by(om, type, sim) %>%
	mutate(catch.rel = catch_fit / catch_sim)
	# mutate(catch.fit = mods[[unique(om)]]$rep$pred_catch[,1],
	# 	   catch.rel = pred_catch / catch.fit)
res.catch$om <- factor(res.catch$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
res.catch$type <- factor(res.catch$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
png(here("plots","NAA","5_catch.png"), width=7, height=7, units='in',res=100)
print(ggplot(res.catch, aes(x=year, y=catch.rel)) +
    stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
    stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
	# geom_line(aes(group=sim),color='grey', alpha=.2) +
	stat_summary(fun = "median", geom = "line", color = "red") +
	coord_cartesian(ylim=c(0,2)) +
	xlab("Year") +
	ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
	geom_hline(yintercept = 1, linetype=2, color='black') +
	facet_grid(rows=vars(om), cols=vars(type)) +
	theme_bw())
dev.off()

# Fig 6. Recruitment (sim data) / Recruitment (true data)
simdata <- lapply(1:4, function(x) readRDS(here("data","simdata","SNEMAYT","NAA",paste0("simdata_om",x,".rds"))))
res.R <- results %>% group_by(om, type, sim) %>%
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


