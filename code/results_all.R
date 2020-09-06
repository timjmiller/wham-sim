# Brian Stock
# June 15 2020
# Simulation test WHAM

# source("/home/bstock/Documents/ms/wham-sim/code/results_all.R")
ids = "SNEMAYT"
re = "NAA"

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)
library(cowplot)
inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
source("/home/bstock/Documents/ms/wham-sim/code/plot_NAA.R")
source("/home/bstock/Documents/ms/wham-sim/code/plot_NAA_pars.R")

for(j in 1:length(ids)){
	# bc.type     bias corrected obs only (= 1) or obs + process (= 2)
	# sim.types   simulated obs only (= 1) or obs + process (= 2)
	# n.mods      4 if all NAA models converged
	# multipanel  TRUE makes a 5-panel plot (B, F, relB, relF, recruit), FALSE makes individual plots
	# plot.eps    FALSE does not plot the TMB epsilon results, TRUE does
	if(re[j] == "NAA"){
		plot_NAA(stock.id=ids[j], bc.type=2, sim.types=1:2, n.mods=4, n.sim=100, multipanel=TRUE, plot.eps=FALSE)
		plot_NAA(stock.id=ids[j], bc.type=2, sim.types=1:2, n.mods=4, n.sim=100, multipanel=FALSE, plot.eps=FALSE)
		plot_NAA_pars(stock.id=ids[j], bc.type=2, sim.types=1:2, n.mods=4, n.sim=100, plot.eps=FALSE)
	}
}
