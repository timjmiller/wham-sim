# Brian Stock
# Aug 17 2020
# Simulation test WHAM
# SNEMAYT
# M

# source(here::here("code","v2_fits","1_fit_SNEMAYT_M.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="misreport")
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","v2_fits","SNEMAYT_M")
dir.create(res_dir, showWarnings=FALSE)

# Fit SNEMA yellowtail data
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))

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
                              selectivity=list(model=c("age-specific","logistic","age-specific","logistic","logistic","logistic","logistic","logistic","age-specific"),
                                 initial_pars=list(c(.01,1,1,1,1,1), c(3,3), c(.01,1,1,1,1,1), c(3,3), c(3,3), c(3,3), c(3,3), c(3,3), c(.01,.25,1,1,1,1)),
                                 fix_pars=list(2:6, NULL, 2:6, NULL, NULL, NULL, NULL, NULL, 3:6)),
                              age_comp = "logistic-normal-pool0") 

  # Fit model
  # mod <- fit_wham(input, do.retro=T, do.osa=F, do.proj=F)  
  mod <- fit_wham(input, do.retro=T, do.osa=F, do.proj=T, proj.opts=list(proj.F=rep(0.001, 3)) 
  # mod <- fit_wham(input, do.retro=T, do.osa=F, do.proj=T) 
  if(exists("err")) rm("err") # need to clean this up

  # Save model
  saveRDS(mod, file=file.path(res_dir, paste0(df.mods$Model[m],".rds")))

  # # mod_proj <- project_wham(mod)
  # mod_proj <- project_wham(mod, proj.opts=list(proj.F=rep(0.001, 3)))
  # saveRDS(mod_proj, file=file.path(res_dir, paste0(df.mods$Model[m],"_proj.rds")))  
}

# # check that all models converged, pdHess, and bias correction succeeded
mod.list <- file.path(res_dir,paste0(df.mods$Model,".rds"))
# mod.list <- list.files(path=res_dir,full.names=TRUE)
mods <- lapply(mod.list, readRDS)
pdHess = sapply(mods, function(x) if(x$sdrep$pdHess) TRUE else FALSE)
conv = sapply(mods, function(x) if(x$opt$convergence == 0) TRUE else FALSE)
conv
pdHess

# # m3 projections failed... constant F too high?
# mod <- readRDS("/home/bstock/Documents/ms/wham-sim/results/v2_fits/SNEMAYT_M/M-3.rds")
# mod$err_retro = NULL
# mod$err_proj = NULL
# mod2 <- project_wham(mod, proj.opts=list(proj.F=rep(0.001, 3)))

# # err in m3 retro
# as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=T)$tab)
