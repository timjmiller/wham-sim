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
# # inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
# res_dir=file.path("/home/bstock/Documents/ms/wham-sim/results")
# simdata_dir=file.path("/home/bstock/Documents/ms/wham-sim/data/simdata")
# plots_dir=file.path("/home/bstock/Documents/ms/wham-sim/plots/v2")

# plot.eps    FALSE does not plot the TMB epsilon results, TRUE does
plot_rel_err_pars_multipanel <- function(ids, re, bc.type=2, sim.types=1:2, plot.eps=FALSE,
                      res_dir=file.path(getwd(),"results"), 
                      simdata_dir=file.path(getwd(),"data","simdata"),
                      plots_dir=file.path(getwd(),"plots","v2")){ 
  ind.naa <- which(re == "NAA")
  if(length(ind.naa) == 0) id.naa <- NULL else id.naa <- paste0(ids[ind.naa],"_NAA")
  ind.m <- which(re == "M")
  if(length(ind.m) == 0) id.m <- NULL else id.m <- paste0(ids[ind.m],"_M")
  ind.ecov <- which(re == "Ecov2")
  if(length(ind.ecov) == 0) id.ecov <- NULL else id.ecov <- paste0(ids[ind.ecov],"_Ecov2")
  ind.sel <- which(re == "sel")
  if(length(ind.sel) == 0) id.sel <- NULL else id.sel <- paste0(ids[ind.sel],"_sel")
  id.list <- list(NAA = id.naa, M = id.m, Ecov2=id.ecov, sel=id.sel)
  n.sim = 100
  
  for(ni in 1:length(id.list)){
    n.stocks <- length(id.list[[ni]])
    if(n.stocks > 0){
      if(names(id.list)[ni] == "NAA"){
        par.labs <- c("sigR","sigA","rhoY","rhoA")
        par.labs.expr <- c(expression(sigma[R]), expression(sigma[a]), expression(rho[year]), expression(rho[age]))
        # names(par.labs.expr) <- par.labs
      }
      if(names(id.list)[ni] == "M"){
        # par.labs <- c("sigM","phiY","phiA")
        # par.labs.expr <- c(expression(sigma[M]), expression(phi[Y]), expression(phi[A]))
        par.labs <- c("sigR","sigM","phiY","phiA")
        par.labs.expr <- c(expression(sigma[R]), expression(sigma[M]), expression(phi[year]), expression(phi[age]))
      }
      if(names(id.list)[ni] == "Ecov2"){
        par.labs <- c("sigR","sigA","sigX","phiX","alpha","beta0","beta1","beta2")
        par.labs.expr <- c(expression(sigma[R]), expression(sigma[a]), expression(sigma[x]), expression(varphi[x]), expression(alpha), expression(beta[0]), expression(beta[1]), expression(beta[2]))
      }      
      if(names(id.list)[ni] == "sel"){
        par.labs <- c("sigR","sigA","sigS","phiY","phiA")
        par.labs.expr <- c(expression(sigma[R]), expression(sigma[a]), expression(sigma[Sel]), expression(varphi[year]), expression(varphi[par]))
      }
      n.pars <- length(par.labs)
      
      if(names(id.list)[ni] == "NAA"){
        n.mods=4
        mlabs = c("m1: SCAA (IID)","m2: SCAA (AR1_y)","m3: NAA (IID)","m4: NAA (2D AR1)")
        mlabs_expr = c(expression(paste("NAA-1:")~paste("Indep. ")*R[y]), 
                 expression(paste("NAA-2:")~paste("AR1 ")*R[y]), 
                 expression(paste("NAA-3:")~paste("Indep. ")*N["a,y"]), 
                 expression(paste("NAA-4:")~paste("2D AR1 ")*N["a,y"]))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1:4)
      }
      if(names(id.list)[ni] == "M"){
        n.mods=3
        mlabs = c("m1: none","m2: IID","m3: 2D AR1")
        mlabs_expr = c(expression(paste("M-1:")~paste("None")), 
               expression(paste("M-2:")~paste("Indep.")), 
               expression(paste("M-3:")~paste("2D")~paste("AR1")))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1:3)
      }
      if(names(id.list)[ni] == "Ecov"){
        n.mods=1
        mlabs = c("m1: CPI-Recruitment")
        mlabs_expr = c(expression(paste("m1:")~paste("CPI-Recruitment")))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1)
      }
      if(names(id.list)[ni] == "Ecov2"){
        n.mods=5
        mlabs = c("m1: RW-none","m2: RW-linear","m3: RW-poly","m4: AR1-linear","m5: AR1-poly")
        mlabs_expr = c(expression(paste("Ecov-1:")~paste("RW-none")), 
                       expression(paste("Ecov-2:")~paste("RW-linear")),
                       expression(paste("Ecov-3:")~paste("RW-poly")),
                       expression(paste("Ecov-4:")~paste("AR1-linear")),
                       expression(paste("Ecov-5:")~paste("AR1-poly")))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1:n.mods)
      }
      if(names(id.list)[ni] == "sel"){
        n.mods=3
        mlabs = c("m1: none","m2: IID","m3: 2D AR1")
        mlabs_expr = c(expression(paste("Sel-1:")~paste("None")), 
               expression(paste("Sel-2:")~paste("Indep.")), 
               expression(paste("Sel-3:")~paste("2D")~paste("AR1")))
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
          for(em in om){
          # for(em in 1:n.mods.tmp){
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
                    # tmp$par.true[3:4] <- rev(inv.rho.trans(mod$parList$trans_NAA_rho))
                    tmp$par.true[3:4] <- rev(rho_trans(mod$parList$trans_NAA_rho))
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
                        # tmp$par.true[2:3] <- rev(inv.rho.trans(mod$parList$M_repars[2:3]))
                        tmp$par.true[2:3] <- rev(rho_trans(mod$parList$M_repars[2:3]))
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
                        tmp$par.true[3:4] <- rev(rho_trans(mod$parList$M_repars[2:3]))
                      }
                    }
                  }
                  if(names(id.list)[ni] == "Ecov"){
                    tmp$par.est[1:2] <- exp(s1[rownames(s1) == "log_NAA_sigma",1])[1:2]
                    tmp$par.est[3] <- exp(s1[rownames(s1) == "Ecov_process_pars",1])[2] # Ecov process sigma
                    tmp$par.est[4] <- -1 + 2/(1 + exp(-s1[rownames(s1) == "Ecov_process_pars",1])[3]) # Ecov process cor
                    tmp$par.est[5] <- s1[rownames(s1) == "Ecov_beta",1] # Ecov-Rec link
                    
                    tmp$par.true[1] <- exp(mod$parList$log_NAA_sigma)[1]
                    tmp$par.true[2] <- exp(mod$parList$log_NAA_sigma)[2]
                    tmp$par.true[3] <- exp(mod$parList$Ecov_process_pars[2,1])
                    tmp$par.true[4] <- -1 + 2/(1 + exp(-mod$parList$Ecov_process_pars[3,1]))
                    tmp$par.true[5] <- mod$parList$Ecov_beta[1,1]
                  }
                  if(names(id.list)[ni] == "Ecov2"){ # c("sigR","sigA","sigX","phiX","alpha","beta0","beta1","beta2")
                    tmp$par.est[1:2] <- exp(s1[rownames(s1) == "log_NAA_sigma",1])[1:2]
                    tmp$par.est[3] <- exp(s1[rownames(s1) == "Ecov_process_pars",1])[2] # Ecov process sigma
                    if(em > 3) tmp$par.est[4] <- -1 + 2/(1 + exp(-s1[rownames(s1) == "Ecov_process_pars",1])[3]) # Ecov process cor
                    tmp$par.est[5] <- exp(s1[rownames(s1) == "mean_rec_pars",1])[1]
                    tmp$par.est[6] <- exp(s1[rownames(s1) == "mean_rec_pars",1])[2]
                    if(em > 1) tmp$par.est[7] <- s1[rownames(s1) == "Ecov_beta",1][1] # Ecov-Rec link(s)
                    if(em %in% c(3,5)) tmp$par.est[8] <- s1[rownames(s1) == "Ecov_beta",1][2]
                    
                    tmp$par.true[1] <- exp(mod$parList$log_NAA_sigma)[1]
                    tmp$par.true[2] <- exp(mod$parList$log_NAA_sigma)[2]
                    tmp$par.true[3] <- exp(mod$parList$Ecov_process_pars[2,1])
                    if(om > 3) tmp$par.true[4] <- -1 + 2/(1 + exp(-mod$parList$Ecov_process_pars[3,1]))
                    tmp$par.true[5:6] <- exp(mod$parList$mean_rec_pars)
                    if(om > 1) tmp$par.true[7] <- mod$parList$Ecov_beta[1,1]
                    if(om %in% c(3,5)) tmp$par.true[8] <- mod$parList$Ecov_beta[2,1]
                  }
                  if(names(id.list)[ni] == "sel"){ # c("sigR","sigA","sigS","phiY","phiA")
                    tmp$par.est[1:2] <- exp(s1[rownames(s1) == "log_NAA_sigma",1])[1:2]
                    if(em > 1){
                      sel_repars <- s1[rownames(s1) == "sel_repars",1]
                      tmp$par.est[3] <- exp(sel_repars[1]) # sigS
                      if(em == 3){
                        tmp$par.est[4] <- rho_trans(sel_repars[3]) # phiY
                        tmp$par.est[5] <- rho_trans(sel_repars[2]) # phiA
                      }
                    }
                    
                    tmp$par.true[1] <- exp(mod$parList$log_NAA_sigma)[1]
                    tmp$par.true[2] <- exp(mod$parList$log_NAA_sigma)[2]
                    if(om == 1) tmp$par.true[3:5] <- 0
                    if(om > 1){
                      tmp$par.true[3] <- exp(mod$parList$sel_repars[1,1])
                      if(om == 3) tmp$par.true[4:5] <- rev(rho_trans(mod$parList$sel_repars[1,2:3]))
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
      
      # plot 
      for(ty in sim.types){
        bnds <- qbinom(c(0.025,0.975), n.sim, 0.5)/n.sim # 95% CI bounds for median (sims in a given year
        df.boxplots <- df.par.all[df.par.all$em == df.par.all$om & df.par.all$type==levels(df.par.all$type)[ty],] %>%
          mutate(rel.err = par.est/par.true-1) %>%
          group_by(stock, om2, par.lab)
        df.plot <- df.par.all[df.par.all$em == df.par.all$om & df.par.all$type==levels(df.par.all$type)[ty],] %>%
          mutate(rel.err = par.est/par.true) %>%
          group_by(stock, om2, par.lab) %>%
          summarize(med.rel.err = median(rel.err, na.rm=T)-1, # median
                    rel_err_lo = quantile(rel.err, probs=bnds[1], na.rm=T)-1,
                    rel_err_hi = quantile(rel.err, probs=bnds[2], na.rm=T)-1)
          # summarize(med.rel.err = mean(rel.err, na.rm=T)-1, rel.err.se = sd(rel.err)/sqrt(n.sim), .groups='keep') %>% 
          # mutate(rel_err_lo = med.rel.err - 1.96*rel.err.se,
          #        rel_err_hi = med.rel.err + 1.96*rel.err.se)
        
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
        
        if(names(id.list)[ni] == "Ecov2") png(file.path(plots_dir, paste0("estpar_",names(id.list)[ni],"_",types[ty],".png")), width=10, height=3, units='in',res=300)
        if(names(id.list)[ni] == "sel") png(file.path(plots_dir, paste0("estpar_",names(id.list)[ni],"_",types[ty],".png")), width=7, height=3, units='in',res=300)
        if(n.pars == 4 & names(id.list)[ni] == "NAA") png(file.path(plots_dir, paste0("estpar_",names(id.list)[ni],"_",types[ty],".png")), width=7, height=7.5, units='in',res=300)
        if(n.pars == 4 & names(id.list)[ni] == "M") png(file.path(plots_dir, paste0("estpar_",names(id.list)[ni],"_4par_",types[ty],".png")), width=6, height=5, units='in',res=300)
        if(n.pars == 3) png(file.path(plots_dir, paste0("estpar_",names(id.list)[ni],"_3par_",types[ty],".png")), width=4, height=5, units='in',res=300)
        print(p)
        # grid::grid.text(unit(0.98,"npc"),0.5, label = 'Operating model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()
      }
    }
  }
}
