# Brian Stock
# June 15 2020
# Simulation test WHAM
# revision / v2
#   re-label models
#   grey-scale

# source("/home/bstock/Documents/ms/wham-sim/code/v3_final_figure_pdfs/fig8_fig9_fig10.R")

library(wham)
library(here)
library(tidyverse)
plots_dir = here("plots","v3","final_pdfs")
res_dir=file.path(getwd(),"results","v2_fits")

ids = c("GBhaddock", "SNEMAYT", "GBhaddock")
re = c("NAA","M","sel")
bc.type = 2 # _oepe

library(wham)
library(tidyverse)
library(ggplotFL)
library(ggsci)
library(cowplot)
library(data.table)
inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
rho_trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/get_results.R")
source("/home/bstock/Documents/ms/wham-sim/code/v3_final_figure_pdfs/plot_rel_err.R")

pwidth <- c(6.5, 5, 5)
pheight <- c(8, 8, 8)
fign <- 8:10
for(j in 1:length(ids)){
  results <- get_results(stock.id=ids[j], re=re[j], bc.type=bc.type)
  p <- plot_rel_err(results, stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=2, multipanel=TRUE)
  cairo_pdf(file.path(plots_dir, paste0("fig",fign[j],".pdf")), width=pwidth[j], height=pheight[j])
  print(p)
  dev.off()
}

