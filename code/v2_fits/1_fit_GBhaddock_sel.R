# Brian Stock
# June 15 2020
# Simulation test WHAM

# source(here::here("code","v2_fits","1_fit_GBhaddock_sel.R"))

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
# devtools::load_all("/home/bstock/Documents/wham")
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","v2_fits","GBhaddock_sel")
dir.create(res_dir, showWarnings=FALSE)

# Georges Bank haddock
# see /home/bstock/Documents/wg_MGWG/state-space/GBhaddock/WHAM/fit_wham_models.r
# Fbar.ages = 5:7
# 9 ages
# 4 selectivity blocks, all logistic (default, specified in ASAP file)
# age comps: 7 (logistic normal, treat 0 obs as missing)
asap3 <- read_asap3_dat(here("data","GBhaddock_ASAP.dat"))

# Fit 3 sel re models:
#  1. none
#  2. iid
#  3. 2dar1
df.mods <- data.frame(sel_re = c('none','iid','2dar1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("Sel-",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods

# Fit models
for(m in 1:n.mods){
  sel_inits <- c(floor(asap3$dat$n_ages/2),0.5) 
  input <- prepare_wham_input(asap3, recruit_model = 2,
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor='iid', sigma='rec+1'),
                              selectivity=list(model=rep("logistic",4), re=c(df.mods$sel_re[m], 'none', 'none', 'none'),
                                 initial_pars=rep(list(sel_inits), 4),
                                 fix_pars=list(NULL, NULL, NULL, NULL)),
                              age_comp = "logistic-normal-miss0")
  input$data$Fbar_ages = 5:7

  # Fit model
  mod <- fit_wham(input, do.sdrep=T, do.retro=T, do.osa=F, do.proj=T, proj.opts=list(proj.F=rep(0.001, 3)))

  # Save model
  saveRDS(mod, file=file.path(res_dir, paste0(df.mods$Model[m],".rds")))
}

# check that all models converged, pdHess, and bias correction succeeded
mod.list <- file.path(res_dir,paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)
conv = sapply(mods, function(x) if(x$sdrep$pdHess) TRUE else FALSE)
conv

