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

# source("/home/bstock/Documents/ms/wham-sim/code/bias_correct_oepe/NScod_NAA_oepe/5_plot_REdevs_NScod_NAA.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)

# get results into data frame
n.mods = 4
res_dir <- here("results","bias_correct_oepe","NScod_NAA_oepe")
plots_dir <- here("plots","bias_correct_oepe","NScod_NAA")
mod.list <- file.path(res_dir,paste0("m",1:n.mods,".rds"))
mods <- lapply(mod.list, readRDS)

mlabs = c("m1: SCAA (IID)","m2: SCAA (AR1_y)","m3: NAA (IID)","m4: NAA (2D AR1)")
mlabs_expr = c(expression(paste("m1:")~paste("SCAA")~paste("(IID)")), 
               expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
               expression(paste("m3:")~paste("NAA")~paste("(IID)")), 
               expression(paste("m4:")~paste("NAA")~paste("(2D AR1)")))
mlabs_short <- mlabs
names(mlabs_short) = paste0("m",1:n.mods)

df.mods <- data.frame(NAA_cor = c('iid','ar1_y','iid','2dar1'),
                      NAA_sigma = c('rec','rec','rec+1','rec+1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- 1:n.mods
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods$mlabs <- factor(df.mods$Model, levels=1:n.mods, labels=mlabs_expr)

n_ages <- mods[[1]]$env$data$n_ages
df.NAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
colnames(df.NAA) <- c(paste0("Age_",1:n_ages),"Year","Model")
for(i in 1:length(mods)){
  tmp = as.data.frame(mods[[i]]$rep$NAA_devs)[1:(length(mods[[i]]$years_full)-1),]
  tmp$Year <- mods[[i]]$years_full[-1]
  colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
  tmp$Model = df.mods$mlabs[i]
  df.NAA <- rbind(df.NAA, tmp)
}
df.plot <- df.NAA %>% tidyr::pivot_longer(-c(Year,Model),
          names_to = "Age",
          names_prefix = "Age_",
          names_transform = list(Age = as.integer),
          values_to = "NAA_re")
df.plot$NAA_re[df.plot$Model %in% c("paste(\"m1:\") ~ paste(\"SCAA\") ~ paste(\"(IID)\")", "paste(\"m2:\") ~ paste(\"SCAA (\") * AR1[y] * paste(\")\")") & df.plot$Age > 1] = 0
df.mods$xint = tail(mods[[1]]$years,1) # terminal year dashed line

png(file.path(plots_dir,"8_REdevs_NScod_NAA.png"), width=7, height=7, units='in', res=300)
print(ggplot(df.plot, ggplot2::aes(x=Year, y=Age)) +
      geom_tile(aes(fill=NAA_re)) +
      # geom_vline(data=df.mods, mapping=aes(xintercept = xint), linetype=2, size=.4) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      facet_wrap(vars(Model), ncol=1, labeller = label_parsed, strip.position='right') +
      scale_fill_gradient2(name = "NAA devs", low = scales::muted("blue"), mid = "white", high = scales::muted("red")))
dev.off()

df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F, do.print=F)$tab)
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL
#   xint  dAIC     AIC
# 1 2016 147.2 -2203.4
# 2 2016 119.0 -2231.6
# 3 2016  15.8 -2334.8
# 4 2016   0.0 -2350.6

