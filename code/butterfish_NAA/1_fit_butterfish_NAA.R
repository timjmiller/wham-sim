# Brian Stock
# June 15 2020
# Simulation test WHAM

# source(here::here("code","butterfish_NAA","1_fit_butterfish_NAA.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","butterfish_NAA")
dir.create(res_dir, showWarnings=FALSE)
ins_dir <- here("data","simdata","butterfish_NAA")

# Butterfish

# 28 years
# 5 ages
# 4 selectivity blocks
asap3 <- read_asap3_dat(here("data","butterfish.dat")) # more recent assessment uses asap4... use 2017 asap3 file for now

# Fit 4 NAA operating models:
#  1. rec 	iid
#  2. rec 	ar1_y
#  3. rec+1 iid
#  4. rec+1 2dar1
df.mods <- data.frame(NAA_cor = c('iid','ar1_y','iid','2dar1'),
                      NAA_sigma = c('rec','rec','rec+1','rec+1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods

# Fit models
for(m in 1:n.mods){
  # input <- prepare_wham_input(asap3, recruit_model = 2,
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]))

  # # # as specified in butterfish.dat
  # input <- prepare_wham_input(asap3, recruit_model = 2,
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]),
  #                             selectivity=list(model=rep("age-specific",4),
  #                                initial_pars=list(c(1,1,1,1,1), c(1,0.58,0.632,0.632,0.632), c(1,0.461,0.657,0.349,0.349), c(1,1,0.298,0.298,0.298)),
  #                                fix_pars=list(3:5, c(1,4,5), c(1,5), c(1,4,5))))

  input <- prepare_wham_input(asap3, recruit_model = 2,
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]),
                              selectivity=list(model=rep("age-specific",4),
                                 initial_pars=list(c(1,1,1,1,1), c(1,1,1,0.632,0.632), c(1,1,0.657,0.349,0.349), c(1,1,1,0.298,0.298)),
                                 fix_pars=list(3:5, 1:3, c(1,2,5), c(1,2,3))))

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

  input$data$Fbar_ages = seq(asap3$dat$Frep_ages[1], asap3$dat$Frep_ages[2])
  input$par$log_N1_pars = log(asap3$dat$N1_ini)
  input$par$log_F1 = log(asap3$dat$F1_ini)

  # Fit model
  # mod <- fit_wham(input, do.retro=F, do.osa=F, do.proj=F) # no TMB bias correction
  mod <- fit_wham(input, do.sdrep=F, do.retro=F, do.osa=F, do.proj=F, do.check=T)  
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

