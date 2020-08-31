# Brian Stock
# June 15 2020
# Simulation test WHAM

# source(here::here("code","GBhaddock_M","1_fit_GBhaddock_M.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","GBhaddock_M")
dir.create(res_dir, showWarnings=FALSE)
ins_dir <- here("data","simdata","GBhaddock_M")

# Georges Bank haddock
# see /home/bstock/Documents/wg_MGWG/state-space/GBhaddock/WHAM/fit_wham_models.r
# Fbar.ages = 5:7
# 9 ages
# 4 selectivity blocks, all logistic (default, specified in ASAP file)
# age comps: 7 (logistic normal, treat 0 obs as missing)
asap3 <- read_asap3_dat(here("data","GBhaddock_ASAP.dat"))

# Fit 3 M operating models (best models from 2D-AR1-survival paper
#  1. None
#  2. IID
#  3. 2D AR1
df.mods <- data.frame(M_re = c('none','iid','2dar1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods

# Fit models
for(m in 1:n.mods){
  input <- prepare_wham_input(asap3, recruit_model = 2,
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor="iid", sigma="rec"), # m1/base NAA model
                              M = list(re=df.mods$M_re[m]))
  # input <- prepare_wham_input(asap3, recruit_model = 2,
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]),
  #                             selectivity=list(model=c("logistic","age-specific","age-specific","logistic"),
  #                                initial_pars=list(c(3,3), c(.5,.5,.5,1,.5,.5,1,1,0), c(.5,1,.5,.5,1,1,1,1,1), c(3,3)),
  #                                fix_pars=list(NULL, c(4,8,9), c(2,5), NULL)))

  # age comp = 7, logistic normal, treat 0 obs as missing, 1 par
  input$data$age_comp_model_indices = rep(7, input$data$n_indices)
  input$data$age_comp_model_fleets = rep(7, input$data$n_fleets)
  input$data$n_age_comp_pars_indices = rep(1, input$data$n_indices)
  input$data$n_age_comp_pars_fleets = rep(1, input$data$n_fleets)
  input$par$index_paa_pars = rep(0, input$data$n_indices)
  input$par$catch_paa_pars = rep(0, input$data$n_fleets)
  input$map = input$map[!(names(input$map) %in% c("index_paa_pars", "catch_paa_pars"))]

  # analytical bias correction, both obs and process error
  input$data$bias_correct_oe = 1
  input$data$bias_correct_pe = 1

  # input$data$Fbar_ages = 5:7
  input$data$Fbar_ages = seq(asap3$dat$Frep_ages[1], asap3$dat$Frep_ages[2])
  input$par$log_N1_pars = log(asap3$dat$N1_ini)
  input$par$log_NAA = matrix(input$par$log_N1_pars, ncol=asap3$dat$n_ages, nrow=asap3$dat$n_years-1, byrow=TRUE)
  input$par$log_F1 = log(asap3$dat$F1_ini)  

  # if(m > 1){
  #   m1 <- readRDS(file.path(res_dir, "m1.rds"))
  #   input$par$logit_selpars <- m1$parList$logit_selpars
  #   input$map$logit_selpars <- factor(rep(NA, length(input$par$logit_selpars)))
  # }

  # Fit model
  # mod <- fit_wham(input, do.retro=F, do.osa=F, do.proj=F) # no TMB bias correction
  mod <- fit_wham(input, do.sdrep=F, do.retro=F, do.osa=F, do.proj=F)  
  if(exists("err")) rm("err") # need to clean this up
  mod$sdrep = TMB::sdreport(mod, bias.correct=TRUE) # also do bias correction
  # simdata <- mod$simulate(par=mod$env$last.par.best, complete=TRUE)

  # Save model
  if(exists("err")) rm("err") # need to clean this up
  saveRDS(mod, file=file.path(res_dir, paste0(df.mods$Model[m],".rds")))
  saveRDS(input, file=file.path(ins_dir, paste0(df.mods$Model[m],"_input.rds")))
}

# check that all models converged, pdHess, and bias correction succeeded
mod.list <- file.path(res_dir,paste0("m",1:n.mods,".rds"))
mods <- lapply(mod.list, readRDS)
conv = sapply(mods, function(x) if(x$sdrep$pdHess) TRUE else FALSE)
conv

# a50 par is -15 for 2nd block... try age-specific?
