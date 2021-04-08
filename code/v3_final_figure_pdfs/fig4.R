# Brian Stock
# Jan 4 2021
# Plot NAA devs

# source("/home/bstock/Documents/ms/wham-sim/code/v3_final_figure_pdfs/fig4.R")

# mods is a list of wham model fits
plot_NAAdevs <- function(mods){
  n.mods = length(mods)
  mlabs = c("Base: SCAA","NAA-1: Indep. R_y","NAA-2: AR1 R_y","NAA-3: Indep. N_a,y","NAA-4: 2D AR1 N_a,y")
  mlabs_expr = c(expression(paste("Base:")~paste("SCAA")),
                 expression(paste("NAA-1:")~paste("Indep. ")*R[y]), 
                 expression(paste("NAA-2:")~paste("AR1 ")*R[y]), 
                 expression(paste("NAA-3:")~paste("Indep. ")*N["a,y"]), 
                 expression(paste("NAA-4:")~paste("2D AR1 ")*N["a,y"]))
  
  df.mods <- data.frame(NAA_cor = c('none','iid','ar1_y','iid','2dar1'),
                        NAA_sigma = c('none','rec','rec','rec+1','rec+1'), stringsAsFactors=FALSE)
  df.mods$Model <- 1:n.mods
  df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
  df.mods$mlabs <- factor(df.mods$Model, levels=1:n.mods, labels=mlabs_expr)

  n_ages <- mods[[1]]$env$data$n_ages
  n_years <- length(mods[[1]]$years)
  df.NAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
  colnames(df.NAA) <- c(paste0("Age_",1:n_ages),"Year","Model")
  for(i in 1:n.mods){
    tmp = as.data.frame(mods[[i]]$rep$NAA_devs)[1:(length(mods[[i]]$years_full)-1),]
    tmp$Year <- mods[[i]]$years_full[-1]
    colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
    tmp$Model = df.mods$mlabs[i]
    df.NAA <- rbind(df.NAA, tmp)
  }
  df.plot <- df.NAA %>% tidyr::pivot_longer(-c(Year,Model),
            names_to = "Age",
            names_prefix = "Age_",
            names_transform = list(Age = as.integer),
            values_to = "NAA_re")
  df.plot$NAA_re[df.plot$Model == "paste(\"Base:\") ~ paste(\"SCAA\")"] <- df.plot$NAA_re[df.plot$Model == "paste(\"Base:\") ~ paste(\"SCAA\")"] - mean(mods[[1]]$rep$NAA_devs[1:n_years,1])
  df.plot$NAA_re[df.plot$Model == "paste(\"Base:\") ~ paste(\"SCAA\")" & df.plot$Year > tail(mods[[1]]$years,1)] = 0
  df.plot$NAA_re[df.plot$Model %in% c("paste(\"Base:\") ~ paste(\"SCAA\")","paste(\"NAA-1:\") ~ paste(\"Indep. \") * R[y]", "paste(\"NAA-2:\") ~ paste(\"AR1 \") * R[y]") & df.plot$Age > 1] = 0
  df.plot$xint = tail(mods[[1]]$years,1)+0.5 # terminal year dashed line

  g <- ggplot(df.plot, ggplot2::aes(x=Year, y=Age)) +
        geom_tile(aes(fill=NAA_re)) +
        geom_vline(aes(xintercept = xint), linetype=2, size=.4) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        facet_wrap(vars(Model), ncol=1, labeller = label_parsed, strip.position='right') +
        scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white", high = scales::muted("red"))
  return(g)
}

library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)

plots_dir = here("plots","v3","final_pdfs")
ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring")
# for(j in 1:length(ids)){
j = 4
  res_dir <- here("results","v2_fits",paste0(ids[j],"_NAA"))
  mod.list <- list.files(res_dir, full.names = TRUE)
  mods <- lapply(mod.list, readRDS)
  g <- plot_NAAdevs(mods)

  cairo_pdf(file.path(plots_dir, "fig4_naadevs.pdf"), width=6, height=7)
  print(g)
  dev.off()
# }

