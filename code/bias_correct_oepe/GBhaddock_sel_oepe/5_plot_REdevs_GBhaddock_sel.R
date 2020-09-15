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

# source("/home/bstock/Documents/ms/wham-sim/code/bias_correct_oepe/GBhaddock_sel_oepe/5_plot_REdevs_GBhaddock_sel.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)
library(viridis)

# get results into data frame
n.mods = 3
res_dir <- here("results","bias_correct_oepe","GBhaddock_sel_oepe")
plots_dir <- here("plots","bias_correct_oepe","GBhaddock_sel")
mod.list <- file.path(res_dir,paste0("m",1:n.mods,".rds"))
mods <- lapply(mod.list, readRDS)
selAA <- lapply(mods, function(x) x$rep$selAA[[1]])

mlabs = c("m1: none","m2: IID","m3: 2D AR1")
mlabs_expr = c(expression(paste("m1:")~paste("none")), 
               expression(paste("m2:")~paste("IID")), 
               expression(paste("m3:")~paste("2D")~paste("AR1")))
mlabs_short <- mlabs
names(mlabs_short) = paste0("m",1:3)

df.mods <- data.frame(Model = 1:n.mods, stringsAsFactors=FALSE)
n_ages <- mods[[1]]$env$data$n_ages
df.selAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
colnames(df.selAA) <- c(paste0("Age_",1:n_ages),"Year","Model")
for(m in 1:n.mods){
	tmp <- as.data.frame(selAA[[m]])
	tmp$Year <- mods[[1]]$years
	colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
	tmp$Model <- df.mods$Model[m]
	df.selAA <- rbind(df.selAA, tmp)
}
df.plot <- df.selAA %>% pivot_longer(-c(Year,Model),
				names_to = "Age", 
				names_prefix = "Age_",
				names_transform = list(Age = as.integer),
				values_to = "Selectivity")
df.plot$mlabs <- factor(df.plot$Model, levels=1:n.mods, labels=mlabs_expr)

png(file.path(plots_dir,"8_REdevs_GBhaddock_sel.png"), width=7, height=5, units='in', res=300)
print(ggplot(df.plot, ggplot2::aes(x=Year, y=Age, fill=Selectivity)) +
      geom_tile() +
      scale_x_continuous(expand=c(0,0), breaks=seq(1935, 2015, by=10)) +
      scale_y_continuous(expand=c(0,0), breaks=seq(2,8,by=2)) +
      theme_bw() +
      facet_wrap(vars(mlabs), ncol=1, labeller = label_parsed, strip.position='right') +
	  scale_fill_viridis(limits=c(0,1)))
dev.off()

df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F, do.print=F)$tab)
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL
#  Model  dAIC     AIC
# 1     1 453.9 -5691.1
# 2     2 177.1 -5967.9
# 3     3   0.0 -6145.0

