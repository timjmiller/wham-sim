# Brian Stock
# June 15 2020
# Simulation test WHAM

# source(here::here("code","bias_correct_oe","NScod_M_oe","1_fit_NScod_M.R"))

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","bias_correct_oe","NScod_M_oe")
dir.create(res_dir, showWarnings=FALSE)
ins_dir <- here("data","simdata","bias_correct_oe","NScod_M_oe")

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
df.mods$Model <- paste0("m",1:n.mods)
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
                                 fix_pars=list(NULL, 5:6, 3:6)))
  # input <- prepare_wham_input(asap3, recruit_model = 2,
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = list(cor="iid", sigma="rec"), # m1/base NAA model
  #                             M = list(re=df.mods$M_re[m]),
  #                             selectivity=list(model=c("logistic","logistic","age-specific"),
  #                                initial_pars=list(c(2,2), c(3,3), c(.8,.8,1,1,0,0)),
  #                                fix_pars=list(NULL, NULL, 3:6)))  

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
  input$data$bias_correct_pe = 0

  input$data$Fbar_ages = 2:4

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

# mods[[3]]$sdrep

# # m3 = 2dar1 didn't converge
# # fix logit_selpars at average from m1 and m2 (similar)
# my.list <- list(mods[[1]]$parList$logit_selpars, mods[[2]]$parList$logit_selpars)
# input$par$logit_selpars <- apply(simplify2array(my.list), 1:2, mean)
# input$map$logit_selpars <- factor(rep(NA, length(input$par$logit_selpars)))
# mod <- fit_wham(input, do.sdrep=F, do.retro=F, do.osa=F, do.proj=F, n.newton=7)  
# if(exists("err")) rm("err") # need to clean this up
# mod$sdrep = TMB::sdreport(mod, bias.correct=TRUE) # also do bias correction

#   saveRDS(mod, file=file.path(res_dir, "mod.rds"))
#   saveRDS(input, file=file.path(ins_dir, "mod_input.rds"))


# mod <- fit_wham(input, do.fit=FALSE)
# mod$opt <- stats::nlminb(mod$par, mod$fn, mod$gr, control = list(iter.max = 10000, eval.max = 10000))
# mod$sdrep = TMB::sdreport(mod, bias.correct=TRUE) # also do bias correction
