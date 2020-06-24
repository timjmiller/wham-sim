# Brian Stock
# June 15 2020
# Simulation test WHAM

# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "om_mode")
library(wham)
library(here)
library(tidyverse)

# Fit SNEMA yellowtail data
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))

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
  # blocks 1, 3, 9 have issues, try age-specific
  input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]),
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

  # 5/29/20, bias correct CAA and IAA
  input$data$bias_correct_oe = 1

  # Fit model with projections:
  #  - 3 years
  #  - F = 0
  #  - WAA, MAA, maturity fixed at terminal year (2011)
  mod <- fit_wham(input, do.retro=F, do.osa=F, do.proj=F)  
  # mod <- fit_wham(input, do.retro=F, do.osa=F, do.proj=T, proj.opts=list(n.yrs=3, use.last.F=FALSE, use.avg.F=FALSE, use.FXSPR=FALSE,
  #                                             proj.F=rep(0,3), proj.catch=NULL, avg.yrs=NULL,
  #                                             cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL, cont.Mre=NULL))

  # Save model
  if(exists("err")) rm("err") # need to clean this up
  saveRDS(mod, file=here("results","NAA",paste0(df.mods$Model[m],".rds")))
  saveRDS(input, file=here("results","NAA",paste0(df.mods$Model[m],"_input.rds")))
}

