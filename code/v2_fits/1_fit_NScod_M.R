# Brian Stock
# June 15 2020
# Simulation test WHAM

# source(here::here("code","v2_fits","1_fit_NScod_M.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","v2_fits","NScod_M")
dir.create(res_dir, showWarnings=FALSE)

# North Sea cod
# see /home/bstock/Documents/wg_MGWG/state-space/NScod/WHAM/fit_wham_models.r
# Fbar.ages = 2:4
# selectivity blocks:
#     	 			est_ages 	ages_fix_0 	ages_fix_1
#   1 logistic 
#   2 age-specific 	1-4 		6 			5 
#   3 age-specific  1-2 		5-6 		3-4
# age comps: 7 (logistic normal, treat 0 obs as missing)
asap3 <- read_asap3_dat(here("data","NScod_ASAP.dat"))

# Fit 3 M operating models (best models from 2D-AR1-survival paper
#  1. None
#  2. IID
#  3. 2D AR1
df.mods <- data.frame(M_re = c('none','iid','2dar1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("M-",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods

# Fit models
for(m in 1:n.mods){
  input <- prepare_wham_input(asap3, recruit_model = 2,
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor="iid", sigma="rec"), # m1/base NAA model
                              M = list(re=df.mods$M_re[m]),
                              selectivity=list(model=c("logistic","age-specific","age-specific"),
                                 initial_pars=list(c(2,2), c(.5,.5,.5,.5,1,0), c(.8,.8,1,1,0,0)),
                                 fix_pars=list(NULL, 5:6, 3:6)),
                              age_comp = "logistic-normal-miss0")
  input$data$Fbar_ages = 2:4

  # Fit model
  mod <- fit_wham(input, do.sdrep=T, do.retro=T, do.osa=F, do.proj=T, proj.opts=list(proj.F=rep(0.001, 3)))

  # Save model
  saveRDS(mod, file=file.path(res_dir, paste0(df.mods$Model[m],".rds")))
}

# check that all models converged, pdHess, and bias correction succeeded
mod.list <- file.path(res_dir,paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)
pdHess = sapply(mods, function(x) if(x$sdrep$pdHess) TRUE else FALSE)
conv = sapply(mods, function(x) if(x$opt$convergence == 0) TRUE else FALSE)
conv
pdHess

