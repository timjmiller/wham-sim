# plot NAA random effect parameters (sigR, sigA, rhoY, rhoA)

# # testing
# library(wham)
# library(tidyverse)
# library(ggplotFL)
# library(ggsci)
# library(cowplot)
# ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring", "SNEMAYT","butterfish","NScod")
# re = c(rep("NAA",5), rep("M",3))
# bc.type = 2 # _oepe
# sim.types = 1:2
# inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
# res_dir=file.path("/home/bstock/Documents/ms/wham-sim/results")
# simdata_dir=file.path("/home/bstock/Documents/ms/wham-sim/data/simdata")
# plots_dir=file.path("/home/bstock/Documents/ms/wham-sim/plots")

# plot.eps    FALSE does not plot the TMB epsilon results, TRUE does
plot_rel_err_pars_multipanel <- function(ids, re, bc.type=2, sim.types=1:2, plot.eps=FALSE,
                      res_dir=file.path(getwd(),"results"), 
                      simdata_dir=file.path(getwd(),"data","simdata"),
                      plots_dir=file.path(getwd(),"plots")){ 
  ind.naa <- which(re == "NAA")
  id.naa <- paste0(ids[ind.naa],"_NAA")
  ind.m <- which(re == "M")
  id.m <- paste0(ids[ind.m],"_M")
  id.list <- list(NAA = id.naa, M = id.m)
  n.sim = 100
  
  for(ni in 1:length(id.list)){
    n.stocks <- length(id.list[[ni]])
    
    if(names(id.list)[ni] == "NAA"){
      par.labs <- c("sigR","sigA","rhoY","rhoA")
      par.labs.expr <- c(expression(sigma[R]), expression(sigma[A]), expression(rho[Y]), expression(rho[A]))
      # names(par.labs.expr) <- par.labs
    }
    if(names(id.list)[ni] == "M"){
      # par.labs <- c("sigM","phiY","phiA")
      # par.labs.expr <- c(expression(sigma[M]), expression(phi[Y]), expression(phi[A]))
      par.labs <- c("sigR","sigM","phiY","phiA")
      par.labs.expr <- c(expression(sigma[R]), expression(sigma[M]), expression(phi[Y]), expression(phi[A]))
    }
    n.pars <- length(par.labs)
    
    if(names(id.list)[ni] == "NAA"){
      n.mods=4
      mlabs = c("m1: SCAA (IID)","m2: SCAA (AR1_y)","m3: NAA (IID)","m4: NAA (2D AR1)")
      mlabs_expr = c(expression(paste("m1:")~paste("SCAA")~paste("(IID)")), 
                     expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
                     expression(paste("m3:")~paste("NAA")~paste("(IID)")), 
                     expression(paste("m4:")~paste("NAA")~paste("(2D AR1)")))
      mlabs_short <- mlabs
      names(mlabs_short) = paste0("m",1:4)
    }
    if(names(id.list)[ni] == "M"){
      n.mods=3
      mlabs = c("m1: none","m2: IID","m3: 2D AR1")
      mlabs_expr = c(expression(paste("m1:")~paste("none")), 
                     expression(paste("m2:")~paste("IID")), 
                     expression(paste("m3:")~paste("2D")~paste("AR1")))
      mlabs_short <- mlabs
      names(mlabs_short) = paste0("m",1:3)
    }
    
    par.colnames <- c("om","em","type","sim","par.true","par.est","par.lab","stock")
    df.par.all <- as.data.frame(matrix(NA, ncol = length(par.colnames), nrow = 0))
    colnames(df.par.all) <- par.colnames
    for(st in 1:n.stocks){
      if(bc.type == 1){
        bc <- "bias_correct_oe"
        id <- paste0(id.list[[ni]][st],"_oe")    
      }
      if(bc.type == 2){
        bc <- "bias_correct_oepe"
        id <- paste0(id.list[[ni]][st],"_oepe") 
      }
      res_dir_i <- file.path(res_dir, bc, id)
      
      n.mods.tmp = n.mods
      if(id.list[[ni]][st] == "NScod_M") n.mods.tmp <- 2
      
      df.par <- as.data.frame(matrix(NA, ncol = length(par.colnames), nrow = 0))
      colnames(df.par) <- par.colnames
      for(om in 1:n.mods.tmp){
        mod <- readRDS(file.path(res_dir_i, paste0("m",om,".rds")))
        for(em in 1:n.mods.tmp){
          sdreps <- readRDS(file.path(res_dir_i, paste0("sdreps_om",om,"_em",em,".rds")))
          reps <- readRDS(file.path(res_dir_i, paste0("reps_om",om,"_em",em,".rds")))
          for(ty in 1:2){
            for(i in 1:n.sim){
              s1 <- sdreps[[ty]][[i]]
              tmp <- data.frame(om = om, em = em, type = ty, sim = i,
                                par.true = rep(NA, n.pars), par.est = rep(NA, n.pars),
                                par.lab = par.labs)
              if(class(s1)[1] != "character"){
                if(names(id.list)[ni] == "NAA"){
                  tmp$par.est[1:2] <- exp(s1[rownames(s1) == "log_NAA_sigma",1])[1:2]
                  if(em %in% c(2,4)) tmp$par.est[3] <- s1[rownames(s1) == "NAA_rho_y",1]
                  if(em == 4) tmp$par.est[4] <- s1[rownames(s1) == "NAA_rho_a",1]
                  tmp$par.true[1] <- exp(mod$parList$log_NAA_sigma)[1]
                  if(om %in% 1:2) tmp$par.true[2] = 0 else tmp$par.true[2] = exp(mod$parList$log_NAA_sigma)[2]
                  tmp$par.true[3:4] <- rev(inv.rho.trans(mod$parList$trans_NAA_rho))
                }
                if(names(id.list)[ni] == "M"){
                  if(n.pars == 3){ # M re pars only (no sigR)
                    if(em > 1){ # sigma_M, rho_M_a, rho_M_y aren't in sdrep if no M random effects
                      tmp$par.est[1] <- s1[rownames(s1) == "sigma_M",1]
                      if(em == 3){
                        tmp$par.est[2] <- s1[rownames(s1) == "rho_M_y",1]
                        tmp$par.est[3] <- s1[rownames(s1) == "rho_M_a",1]
                      }
                    }
                    if(om == 1) tmp$par.true <- 0
                    if(om > 1){
                      tmp$par.true[1] <- exp(mod$parList$M_repars[1])
                      tmp$par.true[2:3] <- rev(inv.rho.trans(mod$parList$M_repars[2:3]))
                    }
                  } else { # 4 pars (also plot sigR)
                    tmp$par.est[1] <- exp(s1[rownames(s1) == "log_NAA_sigma",1])[1] # sigR
                    if(em > 1){ # sigma_M, rho_M_a, rho_M_y aren't in sdrep if no M random effects
                      tmp$par.est[2] <- s1[rownames(s1) == "sigma_M",1]
                      if(em == 3){
                        tmp$par.est[3] <- s1[rownames(s1) == "rho_M_y",1]
                        tmp$par.est[4] <- s1[rownames(s1) == "rho_M_a",1]
                      }
                    }
                    tmp$par.true[1] <- exp(mod$parList$log_NAA_sigma)[1]
                    if(om == 1) tmp$par.true[2:4] <- 0
                    if(om > 1){
                      tmp$par.true[2] <- exp(mod$parList$M_repars[1])
                      tmp$par.true[3:4] <- rev(inv.rho.trans(mod$parList$M_repars[2:3]))
                    }
                  }
                }
              }
              df.par <- rbind(df.par, tmp)
            }
          }
        }
      }
      df.par$stock <- strsplit(id.list[[ni]][st], "_")[[1]][1]
      df.par.all <- rbind(df.par.all, df.par)
    }
    
    # only plot m2 and m3 for M (if not plotting sigR)
    plot.mods <- 1:n.mods
    if(names(id.list)[ni] == "M" & n.pars == 3){
      plot.mods <- 2:n.mods
      df.par.all <- df.par.all[df.par.all$om %in% plot.mods & df.par.all$em %in% plot.mods,]
      mlabs <- mlabs[plot.mods]
      mlabs_expr <- mlabs_expr[plot.mods]
      mlabs_short <- mlabs_short[plot.mods]
    }
    
    # now get data frame with all stocks ready for plotting
    types <- c("OE","OEPE")
    tylabs = c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)")
    df.par.all$type <- factor(df.par.all$type, levels=1:2, labels=tylabs)
    df.par.all$om <- factor(df.par.all$om, levels=plot.mods, labels=mlabs)
    df.par.all$em <- factor(df.par.all$em, levels=plot.mods, labels=mlabs)
    df.par.all$em.x <- fct_recode(df.par.all$em, !!!mlabs_short)
    df.par.all$om2 <- factor(df.par.all$om, labels=mlabs_expr)
    df.par.all$par.lab <- factor(df.par.all$par.lab, levels=par.labs)
    df.par.all$par.lab.expr <- factor(df.par.all$par.lab, labels=par.labs.expr)
    
    # plot (NAA or M)
    for(ty in sim.types){
      bnds <- qbinom(c(0.025,0.975), n.sim, 0.5)/n.sim # 95% CI bounds for median (sims in a given year
      df.boxplots <- df.par.all[df.par.all$em == df.par.all$om & df.par.all$type==levels(df.par.all$type)[ty],] %>%
        mutate(rel.err = par.est/par.true-1) %>%
        group_by(stock, om2, par.lab)
      df.plot <- df.par.all[df.par.all$em == df.par.all$om & df.par.all$type==levels(df.par.all$type)[ty],] %>%
        mutate(rel.err = par.est/par.true) %>%
        group_by(stock, om2, par.lab) %>%
        # summarize(med.rel.err = median(rel.err, na.rm=T)-1, # median
        #           rel_err_lo = quantile(rel.err, probs=bnds[1], na.rm=T)-1,
        #           rel_err_hi = quantile(rel.err, probs=bnds[2], na.rm=T)-1)
        summarize(med.rel.err = mean(rel.err, na.rm=T)-1, rel.err.se = sd(rel.err)/sqrt(n.sim), .groups='keep') %>% 
        mutate(rel_err_lo = med.rel.err - 1.96*rel.err.se,
               rel_err_hi = med.rel.err + 1.96*rel.err.se)
      
      # mean in red on top of boxplot
      p <- ggplot(df.boxplots, aes(x=par.lab)) +
        geom_hline(yintercept = 0, linetype=2) +
        geom_boxplot(aes(y=rel.err), fill='grey', outlier.color = NA) +
        geom_linerange(data=df.plot, mapping=aes(ymin=rel_err_lo, ymax=rel_err_hi), size=.7, color='red') +
        geom_point(data=df.plot, mapping=aes(y=med.rel.err), size=2, color='red') +
        scale_x_discrete(labels = par.labs.expr) +
        coord_cartesian(ylim=c(-.75,.75)) +
        xlab("Parameter") +
        ylab("Relative error") +
        facet_grid(cols=vars(om2), rows=vars(stock), labeller = label_parsed) +
        theme_bw() +
        # theme(axis.text.x = element_text(size=12))
      theme(axis.text.x = element_text(size=12), plot.margin = unit(c(0.3,0.1,0.1,0.1), "in"))
      
      plots_dir_i <- file.path(plots_dir, bc)
      if(n.pars == 4 & names(id.list)[ni] == "NAA") png(file.path(plots_dir_i, paste0("estpar_",names(id.list)[ni],"_",types[ty],".png")), width=7, height=7.5, units='in',res=300)
      if(n.pars == 4 & names(id.list)[ni] == "M") png(file.path(plots_dir_i, paste0("estpar_",names(id.list)[ni],"_4par_",types[ty],".png")), width=6, height=5, units='in',res=300)
      if(n.pars == 3) png(file.path(plots_dir_i, paste0("estpar_",names(id.list)[ni],"_3par_",types[ty],".png")), width=4, height=5, units='in',res=300)
      print(p)
      # grid::grid.text(unit(0.98,"npc"),0.5, label = 'Operating model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()
    }
  }
}
