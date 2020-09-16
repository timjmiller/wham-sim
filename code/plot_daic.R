# setwd("/home/bstock/Documents/ms/wham-sim")
# bc.type=2; sim.types=2
# plots_dir = file.path(getwd(),"plots","bias_correct_oepe")
# res_dir=file.path(getwd(),"results")
# simdata_dir=file.path(getwd(),"data","simdata")
# library(wham)
# library(tidyverse)

plot_daic <- function(plots_dir = file.path(getwd(),"plots","bias_correct_oepe"), 
                      res_dir=file.path(getwd(),"results"),
                      simdata_dir=file.path(getwd(),"data","simdata")){ 
  ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring", "SNEMAYT","butterfish","GBhaddock","SNEMAYT")
  re = c(rep("NAA",5), rep("M",2),"sel","Ecov2")
  id <- paste0(ids,"_",re)
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
  n.mods <- c(4,3,3,5)[match(re,c("NAA","M","sel","Ecov2"))]
  mod.files <- lapply(seq_along(res_dir), function(i) file.path(res_dir[i], paste0("m",1:n.mods[i],".rds")))
  mods <- lapply(mod.files, function(x) lapply(x, readRDS))

  get_daic <- function(mods){
    x <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F, do.print=F)$tab)
    minAIC <- min(x$AIC, na.rm=T)
    x$dAIC <- round(x$AIC - minAIC,1)
    return(x$dAIC)
  }
  daic <- list()
  for(i in 1:length(mods)){
    daic[[i]] <- get_daic(mods[[i]])
  }
  df <- data.frame(m = unlist(lapply(n.mods, function(x) 1:x)),
                   daic = unlist(daic),
                   Stock = rep(ids, n.mods),
                   re = rep(re, n.mods))
  df$m[df$re %in% c("M","sel")] = df$m[df$re %in% c("M","sel")] + 4
  df$m[df$re == "Ecov2"] = df$m[df$re == "Ecov2"] + 7
  # mlabs_expr = c(expression(paste("m1:")~paste("SCAA")~paste("(IID)")), 
  #                     expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
  #                     expression(paste("m3:")~paste("NAA")~paste("(IID)")), 
  #                     expression(paste("m4:")~paste("NAA")~paste("(2D AR1)")),
  #                   expression(paste("m1:")~paste("none")), 
  #                     expression(paste("m2:")~paste("IID")), 
  #                     expression(paste("m3:")~paste("2D")~paste("AR1")))
  mlabs_expr2 = c(expression(atop("SCAA","IID")),
                 expression(atop("SCAA",AR1[y])),
                 expression(atop("NAA","IID")),
                 expression(atop("NAA",paste("2D AR1"))),
                 expression("none"),
                 expression("IID"),
                 expression(paste("2D")~paste("AR1")),
                 expression(atop("RW","none")),
                 expression(atop("RW","linear")),
                 expression(atop("RW","poly")),
                 expression(atop("AR1","linear")),
                 expression(atop("AR1","poly")))
  # mlabs_expr2 = c(expression(atop("SCAA","IID")), 
  #                 expression(atop("SCAA",AR1[y])), 
  #                 expression(atop("NAA","IID")), 
  #                 expression(atop("NAA",paste("2D AR1"))),
  #                 expression("none"), 
  #                 expression("IID"), 
  #                 expression(paste("2D")~paste("AR1")),
  #                 expression(paste("RW\nnone")), 
  #                 expression(paste("RW\nlinear")),
  #                 expression(paste("RW\npoly")),
  #                 expression(paste("AR1\nlinear")),
  #                 expression(paste("AR1\npoly")))
  # mlabs_expr3 = c(expression("SCAA\nIID"), 
  #                 expression(paste("SCAA\n",AR1[y])), 
  #                 expression("NAA\nIID"), 
  #                 expression(paste("NAA\n",paste("2D AR1"))),
  #                 expression("none"), 
  #                 expression("IID"), 
  #                 expression(paste("2D")~paste("AR1")))
  df$m <- factor(df$m, levels=1:length(mlabs_expr2), labels=paste0("m",1:length(mlabs_expr2)))
  # df$Model <- factor(df$m, levels=levels(df$m), labels=mlabs_expr)
  df$Model <- factor(df$m, levels=levels(df$m), labels=mlabs_expr2)
  # df$Model <- factor(df$m, levels=levels(df$m), labels=mlabs_expr3)
  df$relab <- factor(df$re, levels=c("NAA","M","sel","Ecov2"), labels=c("NAA","M","Selectivity","CPI-Recruitment"))

  png(file.path(plots_dir, "daic.png"), width=8, height=3, res=300, units='in')
  print(ggplot(df, aes(x=Model, y=daic, shape=Stock)) +
          geom_point(size=3) +
          ylab(expression(Delta*phantom()*AIC)) +
          # scale_y_continuous(expand=c(0.02,0.02)) +
          scale_x_discrete(labels=function(l) parse(text=l)) +
          facet_wrap(vars(relab), nrow=1, scale='free_x') +
          theme_bw() +
          # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          #       strip.text.x = element_text(size = 10)))
        theme(strip.text.x = element_text(size = 10), axis.text.x = element_text(size = 8), 
              axis.title.x = element_text(margin=margin(-10,0,0,0))))
  dev.off()
}
  