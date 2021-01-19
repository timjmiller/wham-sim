# Get results for one stock

get_results <- function(stock.id="SNEMAYT", re="NAA", bc.type=2, 
                        res_dir=file.path(getwd(),"results"), 
                        simdata_dir=file.path(getwd(),"data","simdata")){ 
  id <- paste0(stock.id,"_",re)
  if(bc.type == 1){
    bc <- "bias_correct_oe"
    # plots_dir <- file.path(plots_dir, bc, id)
    id <- paste0(id,"_oe")    
  }
  if(bc.type == 2){
    bc <- "bias_correct_oepe"
    # plots_dir <- file.path(plots_dir, bc, id)
    id <- paste0(id,"_oepe") 
  }
  # dir.create(plots_dir, showWarnings = FALSE)
  res_dir <- file.path(res_dir, bc, id)
  res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
  res.list <- lapply(res.files, readRDS)
  mod.files <- list.files(path=res_dir, pattern = "^m..rds", full.names = TRUE)
  n.mods <- length(mod.files)
  n.sim <- length(res.list[[1]][[1]])
  flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
  results <- do.call(rbind, flatten.nested.list(res.list)) %>% as.data.frame
  results <- sapply(results, as.numeric)
  results <- as.data.frame(results)
  types <- c("OE","OEPE") # simulation type, not bias correction type
  tylabs = c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)")
  results$type <- factor(results$type, levels=1:2, labels=tylabs)
  if(re == "NAA" & n.mods == 4){
    mlabs = paste0("NAA-",1:4)
    # mlabs = c("m1: SCAA (IID)","m2: SCAA (AR1_y)","m3: NAA (IID)","m4: NAA (2D AR1)")
    # mlabs_expr = c(expression(paste("m1:")~paste("SCAA")~paste("(IID)")), 
    #                expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
    #                expression(paste("m3:")~paste("NAA")~paste("(IID)")), 
    #                expression(paste("m4:")~paste("NAA")~paste("(2D AR1)")))
    mlabs_expr = c(expression(paste("NAA-1")), 
                   expression(paste("NAA-2")), 
                   expression(paste("NAA-3")), 
                   expression(paste("NAA-4")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:4)
  }
  if(re == "M" & n.mods == 3){ # M or selectivity
    mlabs = paste0("M-",1:3)
    mlabs_expr = c(expression(paste("M-1")), 
               expression(paste("M-2")), 
               expression(paste("M-3")))
    # mlabs = c("m1: none","m2: IID","m3: 2D AR1")
    # mlabs_expr = c(expression(paste("m1:")~paste("none")), 
    #                expression(paste("m2:")~paste("IID")), 
    #                expression(paste("m3:")~paste("2D")~paste("AR1")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:3)
  }
  if(re == "sel" & n.mods == 3){ # M or selectivity
    mlabs = paste0("Sel-",1:3)
    mlabs_expr = c(expression(paste("Sel-1")), 
               expression(paste("Sel-2")), 
               expression(paste("Sel-3")))    
    # mlabs = c("m1: none","m2: IID","m3: 2D AR1")
    # mlabs_expr = c(expression(paste("m1:")~paste("none")), 
    #                expression(paste("m2:")~paste("IID")), 
    #                expression(paste("m3:")~paste("2D")~paste("AR1")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:3)
  }  
  if(re == "M" & n.mods == 2){ # NScod no m3
    mlabs = paste0("M-",1:2)
    mlabs_expr = c(expression(paste("M-1")), 
               expression(paste("M-2")))    
    # mlabs = c("m1: none","m2: IID")
    # mlabs_expr = c(expression(paste("m1:")~paste("none")), 
    #                expression(paste("m2:")~paste("IID")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:2)
  }
  if(re == "Ecov" & n.mods == 1){
    mlabs = c("m1: CPI-Recruitment")
    mlabs_expr = c(expression(paste("m1:")~paste("CPI-Recruitment")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1)
  }  
  if(re == "Ecov" & n.mods == 2){
    mlabs = c("m1: none","m2: CPI-Recruitment")
    mlabs_expr = c(expression(paste("m1:")~paste("none")), 
                   expression(paste("m2:")~paste("CPI-Recruitment")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:2)
  }
  if(re == "Ecov2" & n.mods == 5){
    mlabs = paste0("Ecov-",1:5)
    mlabs_expr = c(expression(paste("Ecov-1")), 
                   expression(paste("Ecov-2")), 
                   expression(paste("Ecov-3")), 
                   expression(paste("Ecov-4")),
                   expression(paste("Ecov-5")))    
    # mlabs = c("m1: RW-none","m2: RW-linear","m3: RW-poly","m4: AR1-linear","m5: AR1-poly")
    # mlabs_expr = c(expression(paste("m1:")~paste("RW-none")), 
    #                expression(paste("m2:")~paste("RW-linear")),
    #                expression(paste("m3:")~paste("RW-poly")),
    #                expression(paste("m4:")~paste("AR1-linear")),
    #                expression(paste("m5:")~paste("AR1-poly")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:n.mods)
  }
  results$om <- factor(results$om, levels=1:n.mods, labels=mlabs)
  results$em <- factor(results$em, levels=1:n.mods, labels=mlabs)
  results$em.x <- results$em
  results$em2 <- results$em
  # results$om2 <- results$om
  # results$em.x <- fct_recode(results$em, !!!mlabs_short)
  results$om2 <- factor(results$om, labels=mlabs_expr)
  # results$em2 <- factor(results$em, labels=mlabs_expr)	
  
  # calculate relative error
  results$SSB.rel = results$SSB_fit / results$SSB_sim
  results$SSB.rel.bc = results$SSB_fit_bc / results$SSB_sim
  results$F.rel = results$F_fit / results$F_sim
  results$F.rel.bc = results$F_fit_bc / results$F_sim
  results$relB.rel = results$relB_fit / results$relB_sim
  results$relB.rel.bc = results$relB_fit_bc / results$relB_sim
  results$relF.rel = results$relF_fit / results$relF_sim
  results$relF.rel.bc = results$relF_fit_bc / results$relF_sim
  results$catch.rel = results$catch_fit / results$catch_sim
  results$catch.rel.bc = results$catch_fit_bc / results$catch_sim
  
  simdata <- lapply(1:n.mods, function(x) readRDS(file.path(simdata_dir, bc, id, paste0("simdata_om",x,".rds"))))
  results$R.sim = NA
  for(om in 1:n.mods){
    for(em in 1:n.mods){
      for(i in 1:n.sim){
        for(ty in 1:2){
          res.ind <- which(results$om == mlabs[om] & results$em == mlabs[em] & results$sim == i & results$ty == tylabs[ty])
          results$R.sim[res.ind] <- simdata[[om]][[i]][[ty]]$NAA[,1]
        }
      }
    }
  }
  results$R.rel <- results$NAA1 / results$R.sim
  results$R.rel.bc <- results$NAA1_bc / results$R.sim
  
  return(results)
}