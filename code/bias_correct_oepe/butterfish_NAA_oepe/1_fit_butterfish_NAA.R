# Brian Stock
# June 15 2020
# Simulation test WHAM

# source(here::here("code","bias_correct_oepe","butterfish_NAA_oepe","1_fit_butterfish_NAA.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","bias_correct_oepe","butterfish_NAA_oepe")
dir.create(res_dir, showWarnings=FALSE)
ins_dir <- here("data","simdata","bias_correct_oepe","butterfish_NAA_oepe")

# Butterfish

# 28 years
# 5 ages
# 4 selectivity blocks
# asap3 <- read_asap3_dat(here("data","butterfish.dat")) # more recent assessment uses asap4... use 2017 asap3 file for now
asap3 <- read_asap3_dat(here("data","butterfish_2020_asap3.dat")) # asap3 format of 2020 data

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

  # input <- prepare_wham_input(asap3, recruit_model = 2,
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]),
  #                             selectivity=list(model=rep("age-specific",4),
  #                                initial_pars=list(c(1,1,1,1,1), c(1,1,1,0.632,0.632), c(1,1,1,0.349,0.349), c(1,1,1,0.298,0.298)),
  #                                fix_pars=list(3:5, 1:3, c(1,2,3,5), c(1,2,4,5))))

  input <- prepare_wham_input(asap3, recruit_model = 2,
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]),
                              selectivity=list(model=rep("age-specific",4),
                                 initial_pars=list(c(0.5,0.5,1,1,1), c(1,1,1,0.632,0.632), c(1,1,1,0.349,0.349), c(1,1,0.5,0.298,0.298)),
                                 fix_pars=list(3:5, 1:3, c(1,2,3,5), c(1,2,4,5))))  

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
  input$par$log_NAA = matrix(input$par$log_N1_pars, ncol=asap3$dat$n_ages, nrow=asap3$dat$n_years-1, byrow=TRUE)
  input$par$log_F1 = log(asap3$dat$F1_ini)  
  input$par$M0 <- log(1.28762) # value estimated in 2020 assessment (see data/butterfish2020/asap4_rev_grad.rep)

  inv.logit <- function(x, low, upp) return(low + (upp - low)/(1 + exp(-x)))
  logit <- function(x, low, upp){
     p <- (x - low) / (upp - low)
     return(log(p / (1 - p)))
  }
  q_ini <- c(0.123764, 0.0158975, 0.17078) # from data/butterfish2020/asap4_rev_grad.rep
  input$par$logit_q <- logit(q_ini, input$data$q_lower, input$data$q_upper)
  input$map$logit_q <- factor(c(NA, 1, 2)) # 2020 assessment fixed q for index 1

  # # m4 doesn't converge bc sigma_R -> 0. fix at value from m3
  # if(m == 4){
  #   input$par$log_NAA_sigma[1] <- mod$parList$log_NAA_sigma[1]
  #   input$map$log_NAA_sigma <- factor(c(NA,1))
  # }

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
pdHess = sapply(mods, function(x) if(x$sdrep$pdHess) TRUE else FALSE)
conv = sapply(mods, function(x) if(x$opt$convergence == 0) TRUE else FALSE)
conv
pdHess
