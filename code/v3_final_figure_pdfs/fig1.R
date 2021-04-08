# Figure 1
# dAIC by process and stock

library(wham)
library(here)
library(tidyverse)

df <- readRDS(here("plots","old","bias_correct_oepe","daic.rds"))

df$re2 <- factor(df$re, levels=c("NAA","M","sel","Ecov2"), labels=c("NAA","M","Sel","Ecov"))
df$Model <- paste0(df$re2,"-",df$m)
df$relab <- factor(df$re, levels=c("NAA","M","sel","Ecov2"), labels=c("Numbers at age","Natural mortality","Selectivity","Ecov (CPI-Recruitment)"))

grDevices::cairo_pdf(filename=here("plots","v3","final_pdfs","fig1_daic.pdf"), width=9, height=3)
print(ggplot(df, aes(x=Model, y=daic, shape=Stock)) +
      geom_point(size=3) +
      ylab(expression(Delta*phantom()*AIC)) +
      facet_wrap(vars(relab), nrow=1, scales='free') +
      theme_bw() +
      theme(strip.text.x = element_text(size = 10), axis.text.x = element_text(size = 8,angle = 45, vjust = 1, hjust=1), 
          axis.title.x = element_text(margin=margin(-5,0,0,0))))
dev.off()

