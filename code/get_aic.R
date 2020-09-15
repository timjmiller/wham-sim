# setwd("/home/bstock/Documents/ms/wham-sim")
# id.j="SNEMAYT"; re.j="NAA"; bc.type=2; sim.types=1:2
# res_dir=file.path(getwd(),"results")
# simdata_dir=file.path(getwd(),"data","simdata")
# library(wham)
# library(tidyverse)

get_aic <- function(id.j, re.j, bc.type=2, sim.types=1:2,
                     res_dir=file.path(getwd(),"results"), 
                     simdata_dir=file.path(getwd(),"data","simdata")){ 
  id <- paste0(id.j,"_",re.j)
  if(bc.type == 1){
    bc <- "bias_correct_oe"
    id <- paste0(id,"_oe")    
  }
  if(bc.type == 2){
    bc <- "bias_correct_oepe"
    id <- paste0(id,"_oepe") 
  }
  res_dir <- file.path(res_dir, bc, id)
  simdata_dir <- file.path(simdata_dir, bc, id)
  if(re.j == "NAA") n.mods = 4
  if(re.j == "M") n.mods = 3
  if(re.j == "sel") n.mods = 3
  if(re.j == "Ecov") n.mods = 5
  mod.files <- file.path(res_dir, paste0("m",1:n.mods,".rds"))
  mods <- lapply(mod.files, readRDS)
  n.sim <- 100
  
  df.colnames <- c("om","em","sim","aic","id","bc.type","sim.type","re")
  df <- as.data.frame(matrix(NA, ncol = length(df.colnames), nrow = 0))
  colnames(df) <- df.colnames
  for(om in 1:n.mods){
    simdata <- readRDS(file.path(simdata_dir,paste0("simdata_om",om,".rds")))
    for(em in 1:n.mods){
      sdreps <- readRDS(file.path(res_dir, paste0("sdreps_om",om,"_em",em,".rds")))
      reps <- readRDS(file.path(res_dir, paste0("reps_om",om,"_em",em,".rds")))
      for(ty in sim.types){
        for(i in 1:n.sim){
          nfe = length(mods[[em]]$opt$par) # npar in em
          input <- readRDS(file.path(simdata_dir,paste0("m",em,"_input.rds")))
          n.data <- length(input$data)
          input$data <- simdata[[i]][[ty]][1:n.data]
          if(re.j == "NAA"){
            input$data$n_NAA_sigma = length(input$par$log_NAA_sigma)
            if(em > 2) input$data$NAA_sigma_pointers = c(1,rep(2,input$data$n_ages-1)) else input$data$NAA_sigma_pointers = rep(1,input$data$n_ages)
          }
          if(re.j == "M"){
            input$data$M_re_model <- c(1,2,5)[em]
          }
          if(re.j == "sel"){
            input$data$selblock_models_re[1] <- c(1,2,5)[em]
          }
          s1 <- sdreps[[ty]][[i]]
          if(class(s1)[1] != "character"){
            rnames <- rownames(s1)
            s1 <- as.data.frame(s1)
            s1$rname <- rnames
            parlist <- split(s1$Estimate, s1$rname)
            parnames <- unique(names(mods[[em]]$opt$par))
            lst = parlist[names(parlist) %in% parnames]
            optparlist <- unlist(lst[order(match(names(lst), parnames))])
            names(optparlist) <- names(mods[[em]]$opt$par)
            input$par <- mods[[em]]$env$parList(optparlist) 
            z = fit_wham(input, do.fit = F, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F)
            nll = z$fn()
            aic = 2*(nll + nfe)          
          } else {
            aic = NA
          }
          df <- rbind(df, data.frame(om=om, em=em, sim=i, aic=aic, id=id.j, bc.type=bc.type, sim.type=ty, re=re.j))
        }
      }
    }
  }
  return(df)
}
  