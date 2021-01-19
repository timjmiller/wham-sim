# Brian Stock
# June 15 2020
# Simulation test WHAM

# source(here::here("code","v2_fits","1_fit_butterfish_M.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","v2_fits","butterfish_M")
dir.create(res_dir, showWarnings=FALSE)

# Butterfish

# 28 years
# 5 ages
# 4 selectivity blocks
# asap3 <- read_asap3_dat(here("data","butterfish.dat")) # more recent assessment uses asap4... use 2017 asap3 file for now
asap3 <- read_asap3_dat(here("data","butterfish_2020_asap3.dat")) # asap3 format of 2020 data

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
                              NAA_re = list(cor="iid", sigma="rec"), # NAA-1 model
                              M = list(re=df.mods$M_re[m]),
                              selectivity=list(model=rep("age-specific",4),
                                 initial_pars=list(c(0.5,0.5,1,1,1), c(1,1,1,0.632,0.632), c(1,1,1,0.349,0.349), c(1,1,0.5,0.298,0.298)),
                                 fix_pars=list(3:5, 1:3, c(1,2,3,5), c(1,2,4,5))),
                              age_comp = "logistic-normal-miss0")  

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

