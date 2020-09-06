# Brian Stock
# Aug 17 2020
# Simulation test WHAM
# SNEMAYT
# M

# source(here::here("code","bias_correct_oepe","SNEMAYT_M_oepe","1_fit_SNEMAYT_M.R"))

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","bias_correct_oepe","SNEMAYT_M_oepe")
dir.create(res_dir, showWarnings=FALSE)
ins_dir <- here("data","simdata","bias_correct_oepe","SNEMAYT_M_oepe")

# Fit SNEMA yellowtail data
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))

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
  # blocks 1, 3, 9 have issues, try age-specific
  input <- prepare_wham_input(asap3, recruit_model = 2, # Bev Holt recruitment
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor="iid", sigma="rec"), # m1/base NAA model
                              M = list(re=df.mods$M_re[m]),
                              selectivity=list(model=c("age-specific","logistic","age-specific","logistic","logistic","logistic","logistic","logistic","age-specific"),
                                 initial_pars=list(c(.01,1,1,1,1,1), c(3,3), c(.01,1,1,1,1,1), c(3,3), c(3,3), c(3,3), c(3,3), c(3,3), c(.01,.25,1,1,1,1)),
                                 fix_pars=list(2:6, NULL, 2:6, NULL, NULL, NULL, NULL, NULL, 3:6))) 

  # age comp logistic normal pool obs (not multinomial, the default)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # analytical bias correction, both obs and process error
  input$data$bias_correct_oe = 1
  input$data$bias_correct_pe = 1

  # Fit model
  # mod <- fit_wham(input, do.retro=F, do.osa=F, do.proj=F) # no TMB bias correction
  mod <- fit_wham(input, do.sdrep=F, do.retro=F, do.osa=F, do.proj=F)  
  if(exists("err")) rm("err") # need to clean this up
  mod$sdrep = TMB::sdreport(mod, bias.correct=TRUE) # also do bias correction

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

