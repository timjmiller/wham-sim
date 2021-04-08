# Figure 11

# source("/home/bstock/Documents/ms/wham-sim/code/v3_final_figure_pdfs/fig11.R")

library(wham)
library(here)
library(tidyverse)
library(cowplot)

df.aic <- readRDS(here("results",c("bias_correct_oe","bias_correct_oepe")[2],"aic.rds"))
source("/home/bstock/Documents/ms/wham-sim/code/plot_aic_cross.R")
p <- plot_aic_cross(df.aic, bystock=FALSE) # aggregate across stocks

grDevices::cairo_pdf(filename=here("plots","v3","final_pdfs","fig11.pdf"), width=4, height=10)
print(p)
dev.off()

