# Figure 2
# Mohn's rho reduction relative to Base

# source("/home/bstock/Documents/ms/wham-sim/code/v3_final_figure_pdfs/fig2.R")

library(wham)
library(here)
library(tidyverse)
library(RColorBrewer)

mrho <- readRDS(here("results","v2_fits","mrho_naa.rds"))
df <- bind_rows(mrho, .id = "stock") %>% tidyr::pivot_longer(-c(stock,model), names_to = "rho_type", values_to = "rho_val")
df.plot <- df %>% group_by(stock, rho_type) %>% mutate(rho_diff = abs(rho_val) - abs(rho_val[model == "Base"])) %>% filter(model %in% c("NAA-3","NAA-4"))
df.plot$var <- factor(df.plot$rho_type, levels=c("rho_R","rho_SSB","rho_Fbar"), labels=c("`Mohn's`~rho[R]","`Mohn's`~rho[SSB]","`Mohn's`~rho[F]"))

cols <- brewer.pal(5,"Greys")
cols <- cols[-length(cols)]
f <- function(x) {
r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
r
}

grDevices::cairo_pdf(filename=here("plots","v3","final_pdfs","fig2_mrho.pdf"), width=4.5, height=3)
print(ggplot(df.plot, aes(x=model, y=rho_diff)) +
      geom_hline(yintercept=0, linetype="dashed") +
      stat_summary(aes(fill=model), fun.data = f, geom="boxplot", width=.6) +
      scale_fill_manual(values = cols[3:4]) +
      ylab(expression(Change~"in"~"|"*"Mohn's"~rho*"|")) +
      coord_cartesian(ylim=c(-1.1, .1)) +
      xlab("Model") +
      facet_wrap(vars(var), nrow=1, strip.position = "top", labeller = label_parsed) +
      theme_bw() +
      theme(legend.position = "none", 
        strip.text.x = element_text(size = 10), axis.text.x = element_text(size = 8)))
dev.off()

