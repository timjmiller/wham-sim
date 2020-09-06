# plot NAA random effect parameters (sigR, sigA, rhoY, rhoA)

# bc.type     bias corrected obs only (= 1) or obs + process (= 2)
# sim.types   simulated obs only (= 1) or obs + process (= 2)
# n.mods      4 if all NAA models converged
# multipanel  TRUE makes a 5-panel plot (B, F, relB, relF, recruit), FALSE makes individual plots
# plot.eps    FALSE does not plot the TMB epsilon results, TRUE does
plot_NAA_pars <- function(stock.id="SNEMAYT", bc.type=2, sim.types=1:2, n.mods=4, n.sim=100, plot.eps=FALSE,
                      res_dir=file.path(getwd(),"results"), 
                      simdata_dir=file.path(getwd(),"data","simdata"),
                      plots_dir=file.path(getwd(),"plots")){ 
  id <- paste0(stock.id,"_NAA")
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
                  par.true = rep(NA, 4), par.est = rep(NA, 4),
                  par.lab = c("sigR","sigA","rhoY","rhoA"))
          if(class(s1)[1] != "character"){
            # dat <- simdata[[i]][[ty]]
            tmp$par.est[1:2] <- exp(s1[rownames(s1) == "log_NAA_sigma",1])[1:2]
            if(em %in% c(2,4)) tmp$par.est[3] <- s1[rownames(s1) == "NAA_rho_y",1]
            if(em == 4) tmp$par.est[4] <- s1[rownames(s1) == "NAA_rho_a",1]
            tmp$par.true[1] <- exp(mod$parList$log_NAA_sigma)[1]
            if(om %in% 1:2) tmp$par.true[2] = 0 else tmp$par.true[2] = exp(mod$parList$log_NAA_sigma)[2]
            tmp$par.true[3:4] <- rev(inv.rho.trans(mod$parList$trans_NAA_rho))
          }
          df.par <- rbind(df.par, tmp)
        }
      }
    }
  }

  types <- c("OE","OEPE")
  mlabs = c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)")
  tylabs = c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)")
  df.par$om <- factor(df.par$om, levels=1:n.mods, labels=mlabs)
  df.par$em <- factor(df.par$em, levels=1:n.mods, labels=mlabs)
  df.par$type <- factor(df.par$type, levels=1:2, labels=tylabs)
  df.par$em.x <- fct_recode(df.par$em, m1="m1: SCAA (iid)", m2="m2: SCAA (AR1_y)", m3="m3: NAA (iid)", m4="m4: NAA (2D AR1)")
  df.par$par.lab <- factor(df.par$par.lab, levels=c("sigR","sigA","rhoY","rhoA"))

  for(ty in sim.types){
    for(em in 1:n.mods){
      df.plot <- df.par[df.par$em==levels(df.par$em)[em] & df.par$type==levels(df.par$type)[ty],]
      p <- ggplot(df.plot, aes(x=par.est)) +
        geom_histogram() +
        xlab("Estimate") +
        ylab("Density") +
        geom_vline(aes(xintercept = par.true), linetype=2, color='red') +
        facet_grid(rows=vars(om), cols=vars(par.lab), scales='free') +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      
      png(file.path(plots_dir,paste0("7_estpar_em",em,"_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Operating model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Parameter', rot = 0)   # top)
      dev.off()
    }
  }
}
