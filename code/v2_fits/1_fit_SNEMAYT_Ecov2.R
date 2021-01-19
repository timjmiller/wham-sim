# Brian Stock
# June 15 2020
# Simulation test WHAM
# Ecov-Recruitment link, self test only (one model)

# source(here::here("code","v2_fits","1_fit_SNEMAYT_Ecov2.R"))

# devtools::load_all("/home/bstock/Documents/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)

res_dir <- here("results","v2_fits","SNEMAYT_Ecov2")
dir.create(res_dir, showWarnings=FALSE)

# Fit SNEMA yellowtail data
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))

# created in ms/2D-AR1-survival/code/0_explore_data.R
cpi <- read.csv(here("data","CPI_v3.csv"), header=TRUE) 
colnames(cpi) <- c("Year","CPI","CPI_sigma")
cpi$use <- 1
cpi$use[is.nan(cpi$CPI)] = 0 # don't use 2017 (NaN, fall survey missing)

# Operating models
#  NAA: rec+1 iid (fastest, as in Miller 2016)
#  Ecov-Rec link: Bev-Holt limiting (as in Miller 2016)
#  Ecov obs model: fixed (as in Miller 2016)

# Five models
#  1. RW, none
#  2. RW, linear
#  3. RW, poly-2
#  4. AR1, linear
#  5. AR1, poly-2
df.mods <- data.frame(ecov_process = c('rw','rw','rw','ar1','ar1'), 
                      ecov_link = c('linear','linear','poly-2','linear','poly-2'),
                      ecov_how = c(0,2,2,2,2), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("Ecov-",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods

# Fit models
for(m in 1:n.mods){
# for(m in 4:n.mods){
  ecov <- list(
    label = "CPI",
    mean = as.matrix(cpi$CPI),
    logsigma = matrix(log(cpi$CPI_sigma), ncol=1, nrow=dim(cpi)[1]),
    year = cpi$Year,
    use_obs = matrix(cpi$use, ncol=1, nrow=dim(cpi)[1]), # use all obs except 2017 (missing)
    lag = 1, # CPI in year t affects Rec in year t + 1
    process_model = df.mods$ecov_process[m],
    where = "recruit", # CPI affects recruitment
    how = df.mods$ecov_how[m], # 0 = no effect (but still fit Ecov to compare AIC), 2 = limiting
    link_model = df.mods$ecov_link[m])

  input <- prepare_wham_input(asap3, recruit_model = 3, # 2 = random about mean, 3 = Bev Holt, 4 = Ricker
                              model_name = df.mods$Model[m],                         
                              ecov = ecov,
                              NAA_re = list(cor='iid', sigma='rec+1'),
                              selectivity=list(model=c("age-specific","logistic","age-specific","logistic","logistic","logistic","logistic","logistic","age-specific"),
                                 initial_pars=list(c(.01,1,1,1,1,1), c(3,3), c(.01,1,1,1,1,1), c(3,3), c(3,3), c(3,3), c(3,3), c(3,3), c(.01,.25,1,1,1,1)),
                                 fix_pars=list(2:6, NULL, 2:6, NULL, NULL, NULL, NULL, NULL, 3:6)),
                              age_comp = "logistic-normal-pool0") 

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

# mod=mods[[1]]
# mod <- project_wham(mod, proj.opts=list(proj.F=rep(0.001, 3)))
