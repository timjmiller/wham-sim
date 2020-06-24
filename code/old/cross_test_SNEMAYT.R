# Brian Stock
# June 15 2020
# Simulation test WHAM

# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "om_mode")
library(wham)
library(here)
library(tidyverse)

# Fit SNEMA yellowtail data
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))

# Vary NAA models:
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
mod.list <- here("results","NAA",paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)
input.list <- here("results","NAA",paste0(df.mods$Model,"_input.rds"))
inputs <- lapply(input.list, readRDS)

# version 1: keep all parameters (incl NAA) as in fit model, simulate catch + index data (obs error only)
tmp <- inputs[[2]]
tmp$data$simulate_state <- rep(0,4)
mod2_om_fixNAA = fit_wham(tmp, do.fit = FALSE)
tmp$par <- mods[[2]]$env$parList(par=mods[[2]]$env$last.par.best) # get fit pars - correct!
tmp$data <- mod2_om_fixNAA$simulate(par=mods[[2]]$env$last.par.best, complete=TRUE) # simulate data using fit pars
# tmp$data$NAA[,1]
# mods[[2]]$rep$NAA[,1]
newfit1 <- fit_wham(tmp, do.osa=F, do.retro=F, do.proj=F) # refit model with simulated data (start pars at fit opt)
newfit1 <- fit_wham(tmp, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F) # refit model with simulated data (start pars at fit opt)

# version 2: keep parameters except NAA as in fit model, simulate NAA (process error) and catch + index data (obs error)
tmp <- inputs[[2]]
tmp$par <- mods[[2]]$env$parList(par=mods[[2]]$env$last.par.best) # get fit pars - correct!
tmp$data <- mods[[2]]$simulate(par=mods[[2]]$env$last.par.best, complete=TRUE) # simulate data using fit pars
# tmp$data$NAA[,1]
# mods[[2]]$rep$NAA[,1]
newfit2 <- fit_wham(tmp, do.osa=F, do.retro=F, do.proj=F) # refit model with simulated data (start pars at fit opt)

# Now simulate different time-series of NAA
n.sim <- 10
sim.seeds = matrix(sample(1:1000000, n.sim*n.mods, replace = FALSE), n.mods, n.sim) # set random seeds (mod x sim)

dat <- mods[[3]]$simulate(complete=TRUE)
obj2 <- TMB::MakeADFun(dat, mods[[3]]$parList, DLL="wham", random = input$random, map = input$map)
mod2 <- nlminb(obj2$par, obj2$fn, obj2$gr)

mod <- TMB::MakeADFun(input$data, input$par, DLL = "wham", random = input$random, map = input$map)
simdata <- mod$simulate(complete=TRUE)
obj2 <- MakeADFun(data=list(Y=simdata$Y, n_ages=n.ages, n_years=n.years), 
                parameters=params, map=map, random="NAA_devs", DLL="NAA_2Dar1", silent=TRUE)
par.est <- nlminb(obj2$par, obj2$fn, obj2$gr)$par
c(exp(par.est[1]), exp(par.est[2]), exp(par.est[3]), f(par.est[4]), f(par.est[5]))


# om1 = recruitment devs correlated, rho = 0.7
input1 = prepare_wham_om_input(pop1, recruit_model = recruit_model, selectivity=sel.list, NAA_re = NAA.list1)
om1 = fit_wham(input1, do.fit = FALSE)
# default parameter values: obj$env$last.par
dat1 <- om1$simulate(complete=TRUE)

# om2 = recruitment devs IID, rho = 0
input2 = prepare_wham_om_input(pop2, recruit_model = recruit_model, selectivity=sel.list, NAA_re = NAA.list2)
om2 = fit_wham(input2, do.fit = FALSE)
dat2 <- om2$simulate(complete=TRUE)

# par(mfrow=c(2,2))
# plot(dat1$NAA_devs[,1],type='b', main="Rho(age-1) = 0.7")
# plot(dat2$NAA_devs[,1],type='b', main="Rho(age-1) = 0")
# acf(dat1$NAA_devs[,1])
# acf(dat2$NAA_devs[,1])

# ------------------------------------------------------------
# Fit models
fits <- list()

# self-test 1
tmp <- om1$input
tmp$data = dat1
fits[[1]] <- fit_wham(tmp, do.osa=F, do.retro=F)

# self-test 2
tmp <- om2$input
tmp$data = dat2
fits[[2]] <- fit_wham(tmp, do.osa=F, do.retro=F)

# cross-test fit model 1 to simulated data from om2
tmp <- om1$input
tmp$data = dat2
fits[[3]] <- fit_wham(tmp, do.osa=F, do.retro=F)

# cross-test fit model 2 to simulated data from om1
tmp <- om2$input
tmp$data = dat1
fits[[4]] <- fit_wham(tmp, do.osa=F, do.retro=F)


# get parameter estimates and CI
extract_estimate <- function(model, par, alpha, digits=4){
	tmp <- summary(model$sdrep)
	est = tmp[rownames(tmp) == par,1]
	se = tmp[rownames(tmp) == par,2]
	x <- data.frame(par = par,
			   est = round(est, digits), 
			   low = round(est - se*qnorm(1-alpha/2), digits), 
			   high = round(est + se*qnorm(1-alpha/2), digits))
	return(x)
}
# df <- data.frame(matrix(NA, nrow=0, ncol=4))
# colnames(df) <- c("par","est","low","high","true")
df <- rbind(extract_estimate(fits[[1]], "NAA_rho_y", .05), extract_estimate(fits[[1]], "NAA_sigma", .05))
df$true <- round(c(pop1$NAA_rho[2], pop1$NAA_sigma), 4)

df2 <- rbind(extract_estimate(fits[[2]], "NAA_rho_y", .05), extract_estimate(fits[[2]], "NAA_sigma", .05))
df2$true <- round(c(pop2$NAA_rho[2], pop2$NAA_sigma), 4)

df3 <- rbind(extract_estimate(fits[[3]], "NAA_rho_y", .05), extract_estimate(fits[[3]], "NAA_sigma", .05))
df3$true <- round(c(pop2$NAA_rho[2], pop2$NAA_sigma), 4)

df4 <- rbind(extract_estimate(fits[[4]], "NAA_rho_y", .05), extract_estimate(fits[[4]], "NAA_sigma", .05))
df4$true <- round(c(pop1$NAA_rho[2], pop1$NAA_sigma), 4)


scaa_input = prepare_wham_om_input(input, recruit_model = 2, selectivity=sel.list)
om_M_change = change_M_om(om, M_new_ratio = 3, n_ramp_years = 10, year_change = 2009) 
om_M_change_wham = fit_wham(om_M_change, do.fit = FALSE)
om_M_change_wham$report()$MAA
set.seed(123)
#set.seed(sim.seeds[669,15]) 
#fit scaa model with incorrect M
y = om_M_change_wham$simulate(complete= TRUE)
tinput = scaa_input
#tinput = om
tinput$data = y
tfit = fit_wham(tinput, do.osa = FALSE, do.sdrep = FALSE)
mohns_rho(tfit)

#fit state-space model with correct M
tinput = om_M_change
#tinput = om
tinput$data = y
tfit2 = fit_wham(tinput, do.osa = FALSE, do.sdrep = FALSE)
mohns_rho(tfit2)

#fit SCAA model with correct M
om_M_change_scaa_input = change_M_om(scaa_input, M_new_ratio = 3, n_ramp_years = 10, year_change = 2009) 
tinput = om_M_change_scaa_input
#tinput = om
tinput$data = y
tfit3 = fit_wham(tinput, do.osa = FALSE, do.sdrep = FALSE)
mohns_rho(tfit3)


