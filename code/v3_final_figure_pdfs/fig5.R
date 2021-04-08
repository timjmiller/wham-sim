# Brian Stock
# June 15 2020
# Simulation test WHAM

# source("/home/bstock/Documents/ms/wham-sim/code/v3_final_figure_pdfs/fig5.R")

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
res_dir <- here("results","v2_fits","SNEMAYT_M")
plots_dir = here("plots","v3","final_pdfs")
mod.list <- list.files(res_dir, full.names = TRUE)
mods <- lapply(mod.list, readRDS)

MAA <- lapply(mods, function(x) x$rep$MAA)
mlabs = c("M-1: none","M-2: IID","M-3: 2D AR1")
mlabs_expr = c(expression(paste("M-1:")~paste("None")), 
               expression(paste("M-2:")~paste("Indep.")), 
               expression(paste("M-3:")~paste("2D")~paste("AR1")))
mlabs_short <- mlabs
names(mlabs_short) = paste0("M-",1:3)

df.mods <- data.frame(Model = 1:n.mods, stringsAsFactors=FALSE)
n_ages <- mods[[1]]$env$data$n_ages
df.MAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
colnames(df.MAA) <- c(paste0("Age_",1:n_ages),"Year","Model")
for(m in 1:n.mods){
	tmp <- as.data.frame(MAA[[m]])
	tmp$Year <- mods[[1]]$years_full
	colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
	tmp$Model <- df.mods$Model[m]
	df.MAA <- rbind(df.MAA, tmp)
}
df.plot <- df.MAA %>% pivot_longer(-c(Year,Model),
				names_to = "Age", 
				names_prefix = "Age_",
				names_transform = list(Age = as.integer),
				values_to = "M")
df.plot$mlabs <- factor(df.plot$Model, levels=1:n.mods, labels=mlabs_expr)
df.plot$xint = tail(mods[[1]]$years,1)+0.5 # terminal year dashed line

# png(file.path(plots_dir,"logM_SNEMAYT.png"), width=6, height=5, units='in', res=300)
cairo_pdf(filename=file.path(plots_dir,"fig5.pdf"), width=6, height=5)
print(ggplot(df.plot, ggplot2::aes(x=Year, y=Age, fill=log(M))) +
      geom_tile() +
      geom_vline(aes(xintercept = xint), linetype=2) +
      # scale_x_continuous(expand=c(0,0), breaks=seq(1990, 2015, by=5)) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw() +
      facet_wrap(vars(mlabs), ncol=1, labeller = label_parsed, strip.position='right') +
	  scale_fill_viridis(name = ""))
dev.off()

df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F, do.print=F)$tab)
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL
#   Model dAIC     AIC
# 1     1 75.6 -1187.0
# 2     2 42.6 -1220.0
# 3     3  0.0 -1262.6

