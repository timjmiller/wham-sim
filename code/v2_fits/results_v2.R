# Brian Stock
# June 15 2020
# Simulation test WHAM
# revision / v2
#   re-label models
#   grey-scale

# source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/results_v2.R")

setwd("/home/bstock/Documents/ms/wham-sim")
ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring", "SNEMAYT","butterfish","NScod","SNEMAYT","GBhaddock")
re = c(rep("NAA",5), rep("M",3),"Ecov2","sel")
bc.type = 2 # _oepe

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(tidyverse)
library(ggplotFL)
library(ggsci)
library(cowplot)
library(data.table)
inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
rho_trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/get_results.R")
source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/plot_rel_err.R")
# source("/home/bstock/Documents/ms/wham-sim/code/plot_rel_err_pars.R")
source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/plot_rel_err_pars_multipanel.R")
# source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/plot_conv.R")
# source("/home/bstock/Documents/ms/wham-sim/code/get_conv.R")
source("/home/bstock/Documents/ms/wham-sim/code/get_aic.R")
source("/home/bstock/Documents/ms/wham-sim/code/plot_aic_cross.R")
# source("/home/bstock/Documents/ms/wham-sim/code/plot_daic.R")
source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/plot_3panel_SSB_F_R_trends.R")

# ---------------------------------------------------------
# Relative error plots
for(j in 1:length(ids)){
  results <- get_results(stock.id=ids[j], re=re[j], bc.type=bc.type)
  plot_rel_err(results, stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=2, multipanel=TRUE)
	# plot_rel_err(results, stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=1:2, multipanel=FALSE, plot.eps=FALSE)
	# plot_rel_err_pars(stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=1:2, plot.eps=FALSE)
}
plot_rel_err_pars_multipanel(ids=ids, re=re, bc.type=bc.type, sim.types=2, plot.eps=FALSE)

# --------------------------------------------------------
# Convergence plots
# df.colnames <- c("type","om","em","p.conv","id","bc.type","sim.type","re")
# df.conv <- as.data.frame(matrix(NA, ncol = length(df.colnames), nrow = 0))
# colnames(df.conv) <- df.colnames
# for(j in 1:length(ids)){
#   df.conv <- rbind(df.conv, get_conv(stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=1:2))
# }
# # df.conv$id <- sapply(strsplit(df.conv$id,"_"), first)
# saveRDS(df.conv, file.path(getwd(),"results",c("bias_correct_oe","bias_correct_oepe")[bc.type],"conv.rds"))
# df.conv <- readRDS(file.path(getwd(),"results",c("bias_correct_oe","bias_correct_oepe")[bc.type],"conv.rds"))
# plot_conv(df.conv, plots_dir = "/home/bstock/Documents/ms/wham-sim/plots/v2")

# ------------------------------------------------------
# AIC and dAIC plots
# ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring", "SNEMAYT","butterfish","SNEMAYT","GBhaddock")
# re = c(rep("NAA",5), rep("M",2),"Ecov","sel") # remove NScod M bc only 2/3 models fit
# ids = c("GBhaddock")
# re = c("sel")
# ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring", "SNEMAYT","NScod","butterfish","GBhaddock","SNEMAYT")
# re = c(rep("NAA",5), rep("M",3),"sel","Ecov2")
# plot_3panel_SSB_F_R_trends(ids, re)

# df.colnames <- c("om","em","sim","aic","id","bc.type","sim.type","re")
# df.aic <- as.data.frame(matrix(NA, ncol = length(df.colnames), nrow = 0))
# colnames(df.aic) <- df.colnames
# for(j in 1:length(ids)){
#   df.aic <- rbind(df.aic, get_aic(id.j=ids[j], re.j=re[j], bc.type=bc.type, sim.types=1:2))
# }
# saveRDS(df.aic, file.path(getwd(),"results",c("bias_correct_oe","bias_correct_oepe")[bc.type],"aic.rds"))
df.aic <- readRDS(file.path(getwd(),"results",c("bias_correct_oe","bias_correct_oepe")[bc.type],"aic.rds"))
plot_aic_cross(df.aic, plots_dir = file.path(getwd(),"plots",c("bias_correct_oe","bias_correct_oepe")[bc.type]), bystock=TRUE) # plot each stock individually
plot_aic_cross(df.aic, plots_dir = file.path(getwd(),"plots",c("bias_correct_oe","bias_correct_oepe")[bc.type]), bystock=FALSE) # aggregate across stocks

# plot_daic() # dAIC multipanel full fits (all stocks each re)

