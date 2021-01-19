# Brian Stock
# June 15 2020
# Simulation test WHAM

# source(here::here("code","v2_fits","1_fit_ICEherring_NAA.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","v2_fits","ICEherring_NAA")
dir.create(res_dir, showWarnings=FALSE)

# Iceland herring
# see /home/bstock/Documents/wg_MGWG/state-space/ICEherring/WHAM/fit_wham_models.r
# Fbar.ages = 5:10 - 2 #first age is 3 in the model
# 11 ages
# 2 selectivity blocks
#   logistic (default)
#   age-specific fix 5:8 at 1, 9:11 at 0
# age comps: 7 (logistic normal, treat 0 obs as missing)
asap3 <- read_asap3_dat(here("data","ICEherring_ASAP.dat"))

# Fit 5 NAA operating models:
#  Base  rec fixed effects
#  NAA-1 rec  iid
#  NAA-2 rec  ar1_y
#  NAA-3 rec+1 iid
#  NAA-4 rec+1 2dar1
df.mods <- data.frame(NAA_cor = c('none','iid','ar1_y','iid','2dar1'),
                      NAA_sigma = c('none','rec','rec','rec+1','rec+1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- c("Base",paste0("NAA-",1:(n.mods-1)))
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods

# Fit models
for(m in 1:n.mods){
  if(m == 1){ # no NAA re
    input <- prepare_wham_input(asap3, recruit_model = 2,
                                model_name = df.mods$Model[m],                         
                                selectivity=list(model=c("logistic","age-specific"),
                                   initial_pars=list(c(3,3), c(.5,.5,.5,.5,1,1,1,1,0,0,0)),
                                   fix_pars=list(NULL, 5:11)),
                                age_comp = "logistic-normal-miss0")
  } else {
    input <- prepare_wham_input(asap3, recruit_model = 2,
                                model_name = df.mods$Model[m],                         
                                NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]),
                                selectivity=list(model=c("logistic","age-specific"),
                                   initial_pars=list(c(3,3), c(.5,.5,.5,.5,1,1,1,1,0,0,0)),
                                   fix_pars=list(NULL, 5:11)),
                                age_comp = "logistic-normal-miss0")
  }
  input$data$Fbar_ages = 5:10 - 2

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

