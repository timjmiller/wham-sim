# convergence rates + % selected by AIC
#   identical for NAA and M

# stock.id="SNEMAYT"; re="Ecov2"; bc.type=2; sim.types=1:2
# stock.id="butterfish"; re="NAA"; bc.type=2; sim.types=1:2
# # stock.id="GBhaddock"; re="sel"; bc.type=2; sim.types=1:2
# res_dir=file.path(getwd(),"results")
# simdata_dir=file.path(getwd(),"data","simdata")
# plots_dir=file.path(getwd(),"plots")
# library(wham)
# library(tidyverse)

get_conv <- function(stock.id, re, bc.type=2, sim.types=1:2,
                     res_dir=file.path(getwd(),"results"), 
                     simdata_dir=file.path(getwd(),"data","simdata"),
                     plots_dir=file.path(getwd(),"plots")){ 
  id <- paste0(stock.id,"_",re)
  if(bc.type == 1){
    bc <- "bias_correct_oe"
    plots_dir <- file.path(plots_dir, bc, id)
    id <- paste0(id,"_oe")    
  }
  if(bc.type == 2){
    bc <- "bias_correct_oepe"
    plots_dir <- file.path(plots_dir, bc, id)
    id <- paste0(id,"_oepe") 
  }
  dir.create(plots_dir, showWarnings = FALSE)
  res_dir <- file.path(res_dir, bc, id)
  res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
  res.list <- lapply(res.files, readRDS)
  types <- c("OE","OEPE")

  # n.mods <- sqrt(length(res.list))
  flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X) 
  n.sim <- length(res.list[[1]][[1]])
  df.colnames <- c("om","em","type","sim","conv")
  df <- as.data.frame(matrix(NA, ncol = length(df.colnames), nrow = 0))
  colnames(df) <- df.colnames
  for(m in 1:length(res.list)){
    results <- do.call(rbind, flatten.nested.list(res.list[[m]])) %>% as.data.frame
    results <- sapply(results, as.numeric)
    results <- as.data.frame(results)
    for(ty in 1:2){
      for(i in 1:n.sim){
        tmp <- data.frame(om = unique(results$om)[!is.na(unique(results$om))], em = unique(results$em)[!is.na(unique(results$em))], type=ty, sim=i, conv=NA)
        if(class(res.list[[m]][[ty]][[i]])[1] != 'character') tmp$conv <- 1 else tmp$conv <- 0
        df <- rbind(df, tmp)
      }
    }
  }
  # df$em[df$em == 0] = n.mods
  
  # # results as a list of data frames by type
  # res <- list()
  # for(ty in sim.types){
  #   res[[ty]] <- df %>% filter(type == ty) %>% group_by(om, em) %>% 
  #     summarize(p.conv=sum(conv)/n.sim) %>% pivot_wider(names_from=om, values_from=p.conv) %>% select(-1) %>% as.data.frame
  #   colnames(res[[ty]]) <- paste0("om",1:n.mods)
  #   rownames(res[[ty]]) <- paste0("em", 1:n.mods)
  # }
  # names(res) <- types[sim.types]
  
  # results as a data frame with type as column (for plotting all stocks together)
  res <- df %>% group_by(type, om, em) %>% 
    summarize(p.conv=sum(conv)/n.sim) %>%
    mutate(id=stock.id)
  res$bc.type = bc.type
  res$sim.type = types[res$type]
  res$re = re
  
  return(as.data.frame(res))
}
  