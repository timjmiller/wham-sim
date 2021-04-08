# Brian Stock
# Jan 4 2021
# Plot mean % increase in CV for SSB, F, Recruitment 
#   state-space models (NAA-3, NAA-4) have larger uncertainty

# source("/home/bstock/Documents/ms/wham-sim/code/v3_final_figure_pdfs/fig3.R")

library(wham)
library(tidyverse)
plots_dir = here("plots","v3","final_pdfs")
res_dir=file.path(getwd(),"results","v2_fits")

  ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring")
  re = rep("NAA",5)
  id <- paste0(ids,"_",re)
  res_dir <- file.path(res_dir, id)
  mod.files <- lapply(res_dir, function(x) list.files(path=x, pattern="NAA-", full.names=TRUE))
  mods <- lapply(mod.files, function(x) lapply(x, readRDS))
  n.mods <- length(mods[[1]])

  # from plot.cv() in wham_plots_tables.R
  get_CV <- function(mod){
    years_full = mod$years_full
    nyrs <- length(years_full)
    std <- summary(mod$sdrep)
    ssb.ind <- which(rownames(std) == "log_SSB")
    log.ssb <- std[ssb.ind,1]
    ssb.cv <- std[ssb.ind,2]
    R.ind <- which(rownames(std) == "log_NAA_rep")[1:nyrs]
    log.R <- std[R.ind,1]
    R.cv <- std[R.ind,2]
    faa.ind <- which(rownames(std) == "log_FAA_tot")
    log.faa <- matrix(std[faa.ind,1], nyrs, mod$env$data$n_ages)
    faa.cv <- matrix(std[faa.ind,2], nyrs, mod$env$data$n_ages)
    full.f.ind <- cbind(1:nyrs, apply(log.faa,1, function(x) max(which(x == max(x)))))
    log.full.f <- log.faa[full.f.ind]
    full.f.cv <- faa.cv[full.f.ind]
    return(list(ssb.cv, R.cv, full.f.cv))
  }
  df <- data.frame(matrix(NA, nrow=0, ncol=5))
  colnames(df) <- c("var","proj","mean_rel_cv","Model","Stock")
  for(st in 1:length(mods)){ # for each stock
    cv1 <- get_CV(mods[[st]][[1]]) # NAA-1
    for(i in 3:length(mods[[st]])){ # calc CV for NAA-3 and NAA-4 relative to NAA-1
      tmp <- get_CV(mods[[st]][[i]])
      p.cv = purrr::map2(tmp, cv1, ~ (.x / .y)-1) # relative change in CV
      names(p.cv) = c("SSB","Recruitment","F")
      df.tmp <- bind_rows(p.cv) %>% mutate(Yr = row_number(), proj = !Yr %in% 1:length(mods[[st]][[i]]$years)) %>% select(-Yr) %>% 
        pivot_longer(-proj, names_to = "var", values_to = "p_cv") %>% 
        group_by(var, proj) %>%
        summarize(mean_rel_cv = mean(p_cv)) %>% 
        mutate(Model = paste0("NAA-",i), Stock=ids[st]) %>% as.data.frame
      df <- rbind(df, df.tmp)
    }
  }
  df <- df[!(df$var == "F" & df$proj),]
  df$proj <- factor(df$proj, levels=c(F,T), labels=c("Model years","Projection years"))

  library(RColorBrewer)
  cols <- brewer.pal(n.mods+1,"Greys")
  cols <- cols[-length(cols)]
  f <- function(x) {
    r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  cairo_pdf(file.path(plots_dir, "fig3_cv.pdf"), width=3.75, height=5)
  print(ggplot(df, aes(x=Model, y=mean_rel_cv)) +
          geom_hline(aes(yintercept = 0), linetype=2) +
          stat_summary(aes(fill=Model), fun.data = f, geom="boxplot", width=.6) +
          scale_fill_manual(values = cols[3:4]) +
          ylab("Mean change in CV (relative to NAA-1)") +
          facet_grid(cols=vars(proj), rows=vars(var), scales="free_y") +
          expand_limits(y = -.2) +
          theme_bw() +
          theme(legend.position = "none"))
  dev.off()
