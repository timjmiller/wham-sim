# plot NAA random effect parameters (sigR, sigA, rhoY, rhoA)

# # testing
# library(wham)
# library(tidyverse)
# library(ggplotFL)
# library(ggsci)
# library(cowplot)
# stock.id = "SNEMAYT"
# # stock.id = "butterfish"
# re="M"; n.mods=3
# # re="NAA"; n.mods=4
# n.sim=100
# bc.type = 2
# sim.types = 2
# inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
# res_dir=file.path("/home/bstock/Documents/ms/wham-sim/results")
# simdata_dir=file.path("/home/bstock/Documents/ms/wham-sim/data/simdata")
# plots_dir=file.path("/home/bstock/Documents/ms/wham-sim/plots")

# bc.type     bias corrected obs only (= 1) or obs + process (= 2)
# sim.types   simulated obs only (= 1) or obs + process (= 2)
# n.mods      4 if all NAA models converged
# multipanel  TRUE makes a 5-panel plot (B, F, relB, relF, recruit), FALSE makes individual plots
# plot.eps    FALSE does not plot the TMB epsilon results, TRUE does
plot_rel_err_pars <- function(stock.id="SNEMAYT", re="NAA", bc.type=2, sim.types=1:2, n.mods=4, n.sim=100, plot.eps=FALSE,
                      res_dir=file.path(getwd(),"results"), 
                      simdata_dir=file.path(getwd(),"data","simdata"),
                      plots_dir=file.path(getwd(),"plots")){ 
  id <- paste0(stock.id,"_",re)
  if(bc.type == 1){
    bc <- "bias_correct_oe"
    plots_dir <- file.path(plots_dir, bc, id)
    id <- paste0(id,"_oe")    
  }
  if(bc.type == 2){
    bc <- "bias_correct_oepe"
    plots_dir <- file.path(plots_dir, bc, id)
    id <- paste0(id,"_oepe") 
  }
  res_dir <- file.path(res_dir, bc, id)

  if(re == "NAA"){
    par.labs <- c("sigR","sigA","rhoY","rhoA")
    par.labs.expr <- c(expression(sigma[R]), expression(sigma[A]), expression(rho[Y]), expression(rho[A]))
    # names(par.labs.expr) <- par.labs
  }
  if(re == "M"){
    par.labs <- c("sigM","phiY","phiA")
    par.labs.expr <- c(expression(sigma[M]), expression(phi[Y]), expression(phi[A]))
  }
  n.pars <- length(par.labs)
  
  if(re == "NAA" & n.mods == 4){
    mlabs = c("m1: SCAA (IID)","m2: SCAA (AR1_y)","m3: NAA (IID)","m4: NAA (2D AR1)")
    mlabs_expr = c(expression(paste("m1:")~paste("SCAA")~paste("(IID)")), 
                   expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
                   expression(paste("m3:")~paste("NAA")~paste("(IID)")), 
                   expression(paste("m4:")~paste("NAA")~paste("(2D AR1)")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:4)
  }
  if(re == "M" & n.mods == 3){
    mlabs = c("m1: none","m2: IID","m3: 2D AR1")
    mlabs_expr = c(expression(paste("m1:")~paste("none")), 
                   expression(paste("m2:")~paste("IID")), 
                   expression(paste("m3:")~paste("2D")~paste("AR1")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:3)
  }
  if(re == "M" & n.mods == 2){ # NScod no m3
    mlabs = c("m1: none","m2: IID")
    mlabs_expr = c(expression("m1: none"), 
                   expression("m2: IID"))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:2)
  }
  if(re == "Ecov2" & n.mods == 5){
    mlabs = c("m1: RW-none","m2: RW-linear","m3: RW-poly","m4: AR1-linear","m5: AR1-poly")
    mlabs_expr = c(expression(paste("m1:")~paste("RW-none")), 
                   expression(paste("m2:")~paste("RW-linear")),
                   expression(paste("m3:")~paste("RW-poly")),
                   expression(paste("m4:")~paste("AR1-linear")),
                   expression(paste("m5:")~paste("AR1-poly")))
    mlabs_short <- mlabs
    names(mlabs_short) = paste0("m",1:n.mods)
  }  
  
  par.colnames <- c("om","em","type","sim","par.true","par.est","par.lab")
  df.par <- as.data.frame(matrix(NA, ncol = length(par.colnames), nrow = 0))
  colnames(df.par) <- par.colnames
  for(om in 1:n.mods){
    mod <- readRDS(file.path(res_dir, paste0("m",om,".rds")))
    for(em in 1:n.mods){
      sdreps <- readRDS(file.path(res_dir, paste0("sdreps_om",om,"_em",em,".rds")))
      reps <- readRDS(file.path(res_dir, paste0("reps_om",om,"_em",em,".rds")))
      for(ty in 1:2){
        for(i in 1:n.sim){
          s1 <- sdreps[[ty]][[i]]
          tmp <- data.frame(om = om, em = em, type = ty, sim = i,
                  par.true = rep(NA, n.pars), par.est = rep(NA, n.pars),
                  par.lab = par.labs)
          if(class(s1)[1] != "character"){
            if(re == "NAA"){
              tmp$par.est[1:2] <- exp(s1[rownames(s1) == "log_NAA_sigma",1])[1:2]
              if(em %in% c(2,4)) tmp$par.est[3] <- s1[rownames(s1) == "NAA_rho_y",1]
              if(em == 4) tmp$par.est[4] <- s1[rownames(s1) == "NAA_rho_a",1]
              tmp$par.true[1] <- exp(mod$parList$log_NAA_sigma)[1]
              if(om %in% 1:2) tmp$par.true[2] = 0 else tmp$par.true[2] = exp(mod$parList$log_NAA_sigma)[2]
              tmp$par.true[3:4] <- rev(inv.rho.trans(mod$parList$trans_NAA_rho))
            }
            if(re == "M"){
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
            }
          }
          df.par <- rbind(df.par, tmp)
        }
      }
    }
  }

  types <- c("OE","OEPE")
  tylabs = c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)")
  df.par$type <- factor(df.par$type, levels=1:2, labels=tylabs)
  df.par$om <- factor(df.par$om, levels=1:n.mods, labels=mlabs)
  df.par$em <- factor(df.par$em, levels=1:n.mods, labels=mlabs)
  df.par$em.x <- fct_recode(df.par$em, !!!mlabs_short)
  df.par$om2 <- factor(df.par$om, labels=mlabs_expr)
  df.par$par.lab <- factor(df.par$par.lab, levels=par.labs)
  df.par$par.lab.expr <- factor(df.par$par.lab, labels=par.labs.expr)

  if(re == "NAA") plot.mods <- 1:n.mods
  if(re == "M") plot.mods <- 2:n.mods
  for(ty in sim.types){
    for(em in plot.mods){
      df.plot <- df.par[df.par$em==levels(df.par$em)[em] & df.par$type==levels(df.par$type)[ty],]
      p <- ggplot(df.plot, aes(x=par.est)) +
        geom_histogram() +
        xlab("Parameter estimate") +
        ylab("Density") +
        geom_vline(aes(xintercept = par.true), linetype=2, color='red') +
        facet_grid(rows=vars(om2), cols=vars(par.lab.expr), scales='free', labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), strip.text.x = element_text(size = 12),
              plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      
      png(file.path(plots_dir,paste0("7_estpar_em",em,"_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Operating model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Parameter', rot = 0)   # top)
      dev.off()
    }
  }
  
  for(ty in sim.types){
      bnds <- qbinom(c(0.025,0.975), n.sim, 0.5)/n.sim # 95% CI bounds for median (sims in a given year
      df.boxplots <- df.par[df.par$em == df.par$om & df.par$type==levels(df.par$type)[ty],] %>%
        mutate(rel.err = par.est/par.true-1) %>%
        group_by(om2, par.lab)
      df.plot <- df.par[df.par$em == df.par$om & df.par$type==levels(df.par$type)[ty],] %>%
        mutate(rel.err = par.est/par.true) %>%
        group_by(om2, par.lab) %>%
        # summarize(med.rel.err = median(rel.err, na.rm=T)-1, # median
        #           rel_err_lo = quantile(rel.err, probs=bnds[1], na.rm=T)-1,
        #           rel_err_hi = quantile(rel.err, probs=bnds[2], na.rm=T)-1)
        summarize(med.rel.err = mean(rel.err, na.rm=T)-1, rel.err.se = sd(rel.err)/sqrt(n.sim)) %>% 
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
              facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
              theme_bw() +
              theme(axis.text.x = element_text(size=12))
        # theme(axis.text.x = element_text(size=12), plot.margin = unit(c(0.3,0.1,0.1,0.1), "in"))
        
      if(n.pars == 4) png(file.path(plots_dir,paste0("7_estpar_selftest_",types[ty],".png")), width=6, height=2.5, units='in',res=300)
      if(n.pars == 3) png(file.path(plots_dir,paste0("7_estpar_selftest_",types[ty],".png")), width=5, height=2.5, units='in',res=300)
      print(p)
      # grid::grid.text(unit(0.98,"npc"),0.5, label = 'Operating model', rot = 270) # right
      # grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()
    }
}
