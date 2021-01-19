
# source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/plot_retro.R")

library(wham)
library(tidyverse)
plots_dir = file.path(getwd(),"plots","v2") 
res_dir=file.path(getwd(),"results","v2_fits")

  # NAA only 
  ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring")
  re = rep("NAA",5)
  id <- paste0(ids,"_",re)
  
  res_dir <- file.path(res_dir, id)
  # n.mods <- c(4,3,3,5)[match(re,c("NAA","M","sel","Ecov2"))]
  # n.mods[6:7] = 2 # only 2 M models for NScod
  mod.files <- lapply(seq_along(res_dir), function(i) list.files(res_dir[i], full.names=TRUE))
  # mod.files[[6]] = mod.files[[6]][1:2]
  # mod.files[[7]] = mod.files[[7]][1:2]
  mods <- lapply(mod.files, function(x) lapply(x, readRDS))

  get_retro <- function(mods){
    x <- as.data.frame(compare_wham_models(mods, sort=FALSE, do.print=F)$tab)
    rho <- round(x[,3:5],2)
    return(rho)
    # minAIC <- min(x$AIC, na.rm=T)
    # x$dAIC <- round(x$AIC - minAIC,1)
    # return(x$dAIC)
  }
  mrho <- list()
  for(i in 1:length(mods)){
    mrho[[i]] <- get_retro(mods[[i]])
    mrho[[i]]$model = rownames(mrho[[i]])
    # mrho[[i]]$stock = ids[i]
  }
  names(mrho) = ids
  df <- bind_rows(mrho, .id = "stock") %>% tidyr::pivot_longer(-c(stock,model), names_to = "rho_type", values_to = "rho_val")
  df.plot <- df %>% group_by(stock, rho_type) %>% mutate(rho_diff = abs(rho_val) - abs(rho_val[model == "m1"])) %>% filter(model %in% c("m3","m4"))
  df.plot$var <- factor(df.plot$rho_type, levels=c("rho_R","rho_SSB","rho_Fbar"), labels=c("`Mohn's`~rho[R]","`Mohn's`~rho[SSB]","`Mohn's`~rho[F]"))
  df.plot$mod <- factor(df.plot$model, levels=c("m3","m4"), labels=c("NAA-3","NAA-4"))

  library(RColorBrewer)
  cols <- brewer.pal(5,"Greys")
  cols <- cols[-length(cols)]
  f <- function(x) {
    r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  png(file.path(plots_dir, "mohns_rho_reduction_boxplot.png"), width=4.5, height=3, res=300, units='in')
  print(ggplot(df.plot, aes(x=mod, y=rho_diff)) +
          geom_hline(yintercept=0, linetype="dashed") +
          # geom_boxplot(fill='grey') +        
          stat_summary(aes(fill=mod), fun.data = f, geom="boxplot", width=.6) +
          scale_fill_manual(values = cols[3:4]) +
          ylab(expression(Change~"in"~"|"*"Mohn's"~rho*"|")) +
          coord_cartesian(ylim=c(-1.1, .1)) +
          xlab("Model") +
          facet_wrap(vars(var), nrow=1, strip.position = "top", labeller = label_parsed) +
          theme_bw() +
          theme(legend.position = "none", 
            strip.text.x = element_text(size = 10), axis.text.x = element_text(size = 8)))
  dev.off()

  df.plot$rho_diff[df.plot$rho_diff < -1.1] = -1.1
  png(file.path(plots_dir, "mohns_rho_reduction.png"), width=4.5, height=3, res=300, units='in')
  print(ggplot(df.plot, aes(x=mod, y=rho_diff, shape=stock)) +
          geom_point(size=2.5, position=position_dodge(width=0.3)) +
          geom_hline(yintercept=0, linetype="dashed") +
          ylab(expression(Change~"in"~"|"*"Mohn's"~rho*"|")) +
          ylim(-1.1, .1) +
          xlab("Model") +
          # scale_x_discrete(labels=function(l) parse(text=l)) +
          facet_wrap(vars(var), nrow=1, strip.position = "top", labeller = label_parsed) +
          theme_bw() +
          theme(legend.position = "none", 
            strip.text.x = element_text(size = 10), axis.text.x = element_text(size = 8)))
  dev.off()




# mrho
# $SNEMAYT_NAA
#    rho_R rho_SSB rho_Fbar
# m1  5.57    0.85    -0.32
# m2  6.42    1.02    -0.43
# m3  4.56    0.90    -0.38
# m4  0.81    0.17    -0.10
# m5  0.51    0.06     0.00

# $butterfish_NAA
#    rho_R rho_SSB rho_Fbar
# m1  0.18    0.25    -0.16
# m2  0.05    0.10    -0.07
# m3  0.07    0.11    -0.08
# m4  0.05    0.05     0.02
# m5  0.11    0.16    -0.12

# $NScod_NAA
#    rho_R rho_SSB rho_Fbar
# m1 -0.09    0.17    -0.15
# m2 -0.04    0.22    -0.19
# m3 -0.05    0.19    -0.16
# m4  0.02    0.10    -0.11
# m5  0.00    0.05    -0.05

# $GBhaddock_NAA
#    rho_R rho_SSB rho_Fbar
# m1  0.94    0.36    -0.31
# m2  0.89    0.38    -0.33
# m3  0.88    0.36    -0.31
# m4  0.61    0.11    -0.11
# m5  0.31    0.06    -0.07

# $ICEherring_NAA
#    rho_R rho_SSB rho_Fbar
# m1  0.51    0.63    -0.36
# m2 -0.01    0.02    -0.03
# m3  0.01    0.03    -0.04
# m4 -0.07   -0.05     0.06
# m5 -0.05   -0.03     0.03

# $SNEMAYT_M
#    rho_R rho_SSB rho_Fbar
# m1  6.42    1.02    -0.43
# m2  0.12    0.18    -0.11

# $NScod_M
#    rho_R rho_SSB rho_Fbar
# m1 -0.04    0.22    -0.19
# m2  0.31    0.30    -0.29

# $butterfish_M
#    rho_R rho_SSB rho_Fbar
# m1  0.05    0.10    -0.07
# m2  0.05    0.03     0.04
# m3  0.03    0.08    -0.04

# $GBhaddock_sel
#    rho_R rho_SSB rho_Fbar
# m1  0.60    0.12    -0.11
# m2  0.67    0.15    -0.20
# m3  0.45    0.19    -0.24

# $SNEMAYT_Ecov2
#    rho_R rho_SSB rho_Fbar
# m1  0.77    0.18    -0.11
# m2  0.65    0.16    -0.09
# m3  0.58    0.15    -0.09
# m4  0.65    0.16    -0.10
# m5  0.58    0.15    -0.09
