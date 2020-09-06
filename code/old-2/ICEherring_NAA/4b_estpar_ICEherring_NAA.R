# Brian Stock
# June 15 2020
# Simulation test WHAM
# Step 4: Collect simulation fits
#   NEFSC server version 
#   separate results files for each om/em cross test, for XX/YY in 1-4
#     results_omXX_emYY.rds
#     sdreps_omXX_emYY.rds
#     reps_omXX_emYY.rds
#   assume they are copied to /home/bstock/Documents/ms/wham-sim/results/ICEherring/NAA

# source("/home/bstock/Documents/ms/wham-sim/code/old-incorrect-bias-correction/ICEherring_NAA/4b_estpar_ICEherring_NAA.R")
id = "ICEherring_NAA"

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)

# get results into data frame
res_dir <- here("results","old-incorrect-bias-correction",id)
simdata_dir <- here("data","simdata","old-incorrect-bias-correction",id)
plots_dir <- here("plots",id)
inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)

res.colnames <- c("om","em","type","sim","par.true","par.est","par.lab")
df <- as.data.frame(matrix(NA, ncol = length(res.colnames), nrow = 0))
colnames(df) <- res.colnames

n.mods = 4
n.ty = 2
n.sim = 100
for(om in 1:n.mods){
	# simdata <- readRDS(file.path(simdata_dir, paste0("simdata_om",om,".rds")))
	mod <- readRDS(file.path(res_dir, paste0("m",om,".rds")))
	for(em in 1:n.mods){
		sdreps <- readRDS(file.path(res_dir, paste0("sdreps_om",om,"_em",em,".rds")))
		reps <- readRDS(file.path(res_dir, paste0("reps_om",om,"_em",em,".rds")))
		for(ty in 1:n.ty){
			for(i in 1:n.sim){
				s1 <- sdreps[[ty]][[i]]
				tmp <- data.frame(om = om, em = em, type = ty, sim = i,
								par.true = rep(NA, 4), par.est = rep(NA, 4),
								par.lab = c("sigR","sigA","rhoY","rhoA"))
				if(class(s1)[1] != "character"){
					# dat <- simdata[[i]][[ty]]
					tmp$par.est[1:2] <- exp(s1[rownames(s1) == "log_NAA_sigma",1])[1:2]
					if(em %in% c(2,4)) tmp$par.est[3] <- s1[rownames(s1) == "NAA_rho_y",1]
					if(em == 4) tmp$par.est[4] <- s1[rownames(s1) == "NAA_rho_a",1]
					tmp$par.true[1] <- exp(mod$parList$log_NAA_sigma)[1]
					if(om %in% 1:2) tmp$par.true[2] = 0 else tmp$par.true[2] = exp(mod$parList$log_NAA_sigma)[2]
					tmp$par.true[3:4] <- rev(inv.rho.trans(mod$parList$trans_NAA_rho))
				}
				df <- rbind(df, tmp)
			}
		}
	}
}
# df[df == 0] = NA

types <- c("OE","OEPE")
df$om <- factor(df$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
df$em <- factor(df$em, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
df$type <- factor(df$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
df$em.x <- fct_recode(df$em, m1="m1: SCAA (iid)", m2="m2: SCAA (AR1_y)", m3="m3: NAA (iid)", m4="m4: NAA (2D AR1)")
df$par.lab <- factor(df$par.lab, levels=c("sigR","sigA","rhoY","rhoA"))

# Fig 1. SSB (sim fit) / SSB (sim data)
for(ty in 1:length(types)){
	for(em in 1:n.mods){
		df.plot <- df[df$em==levels(df$em)[em] & df$type==levels(df$type)[ty],]
		p <- ggplot(df.plot, aes(x=par.est)) +
			geom_histogram() +
			xlab("Estimate") +
			ylab("Density") +
			geom_vline(aes(xintercept = par.true), linetype=2, color='red') +
			facet_grid(rows=vars(om), cols=vars(par.lab), scales='free') +
			theme_bw() +
			theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
		
		png(file.path(plots_dir,paste0("7_estpar_em",em,"_",types[ty],".png")), width=7, height=7, units='in',res=100)
		print(p)
		grid::grid.text(unit(0.98,"npc"),0.5, label = 'Operating model', rot = 270) # right
		grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Parameter', rot = 0)   # top)
		dev.off()
	}
}
