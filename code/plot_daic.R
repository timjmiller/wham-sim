setwd("/home/bstock/Documents/ms/wham-sim")
bc.type=2; sim.types=2
plots_dir = file.path(getwd(),"plots","bias_correct_oepe")
res_dir=file.path(getwd(),"results")
simdata_dir=file.path(getwd(),"data","simdata")
library(wham)
library(tidyverse)

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

#      m  daic      Stock    re                        Model           relab
# 1   m1 252.6    SNEMAYT   NAA          atop("SCAA", "IID")             NAA
# 2   m2 179.5    SNEMAYT   NAA         atop("SCAA", AR1[y])             NAA
# 3   m3  44.9    SNEMAYT   NAA           atop("NAA", "IID")             NAA
# 4   m4   0.0    SNEMAYT   NAA atop("NAA", paste("2D AR1"))             NAA
# 5   m1  57.6 butterfish   NAA          atop("SCAA", "IID")             NAA
# 6   m2  59.1 butterfish   NAA         atop("SCAA", AR1[y])             NAA
# 7   m3  16.0 butterfish   NAA           atop("NAA", "IID")             NAA
# 8   m4   0.0 butterfish   NAA atop("NAA", paste("2D AR1"))             NAA
# 9   m1 147.2      NScod   NAA          atop("SCAA", "IID")             NAA
# 10  m2 119.0      NScod   NAA         atop("SCAA", AR1[y])             NAA
# 11  m3  15.8      NScod   NAA           atop("NAA", "IID")             NAA
# 12  m4   0.0      NScod   NAA atop("NAA", paste("2D AR1"))             NAA
# 13  m1 401.6  GBhaddock   NAA          atop("SCAA", "IID")             NAA
# 14  m2 378.3  GBhaddock   NAA         atop("SCAA", AR1[y])             NAA
# 15  m3  53.2  GBhaddock   NAA           atop("NAA", "IID")             NAA
# 16  m4   0.0  GBhaddock   NAA atop("NAA", paste("2D AR1"))             NAA
# 17  m1 226.7 ICEherring   NAA          atop("SCAA", "IID")             NAA
# 18  m2 223.7 ICEherring   NAA         atop("SCAA", AR1[y])             NAA
# 19  m3  17.5 ICEherring   NAA           atop("NAA", "IID")             NAA
# 20  m4   0.0 ICEherring   NAA atop("NAA", paste("2D AR1"))             NAA
# 21  m5 233.3    SNEMAYT     M                         none               M
# 22  m6  11.1    SNEMAYT     M                          IID               M
# 23  m7   0.0    SNEMAYT     M   paste("2D") ~ paste("AR1")               M
# 24  m5  75.6 butterfish     M                         none               M
# 25  m6  42.6 butterfish     M                          IID               M
# 26  m7   0.0 butterfish     M   paste("2D") ~ paste("AR1")               M
# 27  m5 453.9  GBhaddock   sel                         none     Selectivity
# 28  m6 177.1  GBhaddock   sel                          IID     Selectivity
# 29  m7   0.0  GBhaddock   sel   paste("2D") ~ paste("AR1")     Selectivity
# 30  m8  33.0    SNEMAYT Ecov2           atop("RW", "none") CPI-Recruitment
# 31  m9  12.7    SNEMAYT Ecov2         atop("RW", "linear") CPI-Recruitment
# 32 m10  14.0    SNEMAYT Ecov2           atop("RW", "poly") CPI-Recruitment
# 33 m11   0.0    SNEMAYT Ecov2        atop("AR1", "linear") CPI-Recruitment
# 34 m12   1.3    SNEMAYT Ecov2          atop("AR1", "poly") CPI-Recruitment
#   