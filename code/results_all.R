# Brian Stock
# June 15 2020
# Simulation test WHAM

# source("/home/bstock/Documents/ms/wham-sim/code/results_all.R")

setwd("/home/bstock/Documents/ms/wham-sim")
ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring", "SNEMAYT","butterfish","NScod")
re = c(rep("NAA",5), rep("M",3))
# ids = c("SNEMAYT","butterfish","NScod")
# re = rep("M",3)
# bc.type = 2 # _oepe
bc.type = 1 # _oe

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(tidyverse)
library(ggplotFL)
library(ggsci)
library(cowplot)
inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
source("/home/bstock/Documents/ms/wham-sim/code/plot_rel_err.R")
source("/home/bstock/Documents/ms/wham-sim/code/plot_rel_err_pars.R")
source("/home/bstock/Documents/ms/wham-sim/code/plot_rel_err_pars_multipanel.R")
source("/home/bstock/Documents/ms/wham-sim/code/plot_conv.R")
source("/home/bstock/Documents/ms/wham-sim/code/get_conv.R")
source("/home/bstock/Documents/ms/wham-sim/code/get_results.R")

for(j in 1:length(ids)){
# for(j in 6:length(ids)){
  results <- get_results(stock.id=ids[j], re=re[j], bc.type=bc.type)
    
	# bc.type     bias corrected obs only (= 1) or obs + process (= 2)
	# sim.types   simulated obs only (= 1) and/or obs + process (= 2)
	# n.mods      4 if all NAA models converged
	# multipanel  TRUE makes a 5-panel plot (B, F, relB, relF, recruit), FALSE makes individual plots
	# plot.eps    FALSE does not plot the TMB epsilon results, TRUE does
	plot_rel_err(results, stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=1:2, multipanel=TRUE)
	plot_rel_err(results, stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=1:2, multipanel=FALSE, plot.eps=FALSE)
	plot_rel_err_pars(stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=1:2, 
	              n.mods=length(table(results$om)), n.sim=length(table(results$sim)), plot.eps=FALSE)
}

plot_rel_err_pars_multipanel(ids=ids, re=re, bc.type=bc.type, sim.types=1:2, plot.eps=FALSE)

df.colnames <- c("type","om","em","p.conv","id","bc.type","sim.type","re")
df.conv <- as.data.frame(matrix(NA, ncol = length(df.colnames), nrow = 0))
colnames(df.conv) <- df.colnames
for(j in 1:length(ids)){
  df.conv <- rbind(df.conv, get_conv(stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=1:2))
}
# df.conv$id <- sapply(strsplit(df.conv$id,"_"), first)
saveRDS(df.conv, file.path(getwd(),"results",c("bias_correct_oe","bias_correct_oepe")[bc.type],"conv.rds"))
plot_conv(df.conv, plots_dir = file.path(getwd(),"plots",c("bias_correct_oe","bias_correct_oepe")[bc.type]))




