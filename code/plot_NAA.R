# plot NAA (one multipanel)

# bc.type     bias corrected obs only (= 1) or obs + process (= 2)
# sim.types   simulated obs only (= 1) or obs + process (= 2)
# n.mods      4 if all NAA models converged
# multipanel  TRUE makes a 5-panel plot (B, F, relB, relF, recruit), FALSE makes individual plots
# plot.eps    FALSE does not plot the TMB epsilon results, TRUE does
plot_NAA <- function(stock.id="SNEMAYT", bc.type=2, sim.types=1:2, n.mods=4, n.sim=100, multipanel=TRUE, plot.eps=FALSE,
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
  res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
  res.list <- lapply(res.files, readRDS)
  flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
  results <- do.call(rbind, flatten.nested.list(res.list)) %>% as.data.frame
  results <- sapply(results, as.numeric)
  results <- as.data.frame(results)
  types <- c("OE","OEPE")
  mlabs = c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)")
  tylabs = c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)")
  results$om <- factor(results$om, levels=1:n.mods, labels=mlabs)
  results$em <- factor(results$em, levels=1:n.mods, labels=mlabs)
  results$type <- factor(results$type, levels=1:2, labels=tylabs)
  results$em.x <- fct_recode(results$em, m1="m1: SCAA (iid)", m2="m2: SCAA (AR1_y)", m3="m3: NAA (iid)", m4="m4: NAA (2D AR1)")
  
  # calculate relative error
  results$SSB.rel = results$SSB_fit / results$SSB_sim
  results$SSB.rel.bc = results$SSB_fit_bc / results$SSB_sim
  results$F.rel = results$F_fit / results$F_sim
  results$F.rel.bc = results$F_fit_bc / results$F_sim
  results$relB.rel = results$relB_fit / results$relB_sim
  results$relB.rel.bc = results$relB_fit_bc / results$relB_sim
  results$relF.rel = results$relF_fit / results$relF_sim
  results$relF.rel.bc = results$relF_fit_bc / results$relF_sim
  results$catch.rel = results$catch_fit / results$catch_sim
  results$catch.rel.bc = results$catch_fit_bc / results$catch_sim
  
  simdata <- lapply(1:n.mods, function(x) readRDS(file.path(simdata_dir, bc, id, paste0("simdata_om",x,".rds"))))
  results$R.sim = NA
  for(om in 1:n.mods){
    for(em in 1:n.mods){
      for(i in 1:n.sim){
        for(ty in 1:2){
          res.ind <- which(results$om == mlabs[om] & results$em == mlabs[em] & results$sim == i & results$ty == tylabs[ty])
          results$R.sim[res.ind] <- simdata[[om]][[i]][[ty]]$NAA[,1]
        }
      }
    }
  }
  results$R.rel <- results$NAA1 / results$R.sim
  results$R.rel.bc <- results$NAA1_bc / results$R.sim
  
  # ---------------------------------------------------------------------------------
  # Multipanel plots
  if(multipanel){
    for(ty in sim.types){
    # collapse across years, group by om/em
    	df.plot <- filter(results, type==levels(results$type)[ty]) %>%
                  select(om, em, em.x, SSB.rel, F.rel, relB.rel, relF.rel, R.rel) %>%
                  pivot_longer(-c(om,em,em.x), names_to = "variable", values_to = "val") %>%
    	            group_by(om, em)
    	df.plot$val = df.plot$val - 1 # relative error
    	
    	df.plot$variable <- factor(df.plot$variable, levels=c("SSB.rel", "F.rel", "relB.rel", "relF.rel", "R.rel"), 
    	                       labels=c("SSB", "F", expression(B/B[40]["%"]), expression(F/F[40]["%"]), "Recruitment"))
    	df.plot$om2 <- factor(df.plot$om, labels=c(expression(paste("m1:")~paste("SCAA")~paste("(iid)")), 
    	                                     expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
                	                         expression(paste("m3:")~paste("NAA")~paste("(iid)")), 
    	                                     expression(paste("m4:")~paste("NAA")~paste("(2D AR1)"))))
    	df.plot$em2 <- factor(df.plot$em, labels=c(expression(paste("m1:")~paste("SCAA")~paste("(iid)")), 
    	                                     expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
                	                         expression(paste("m3:")~paste("NAA")~paste("(iid)")), 
    	                                     expression(paste("m4:")~paste("NAA")~paste("(2D AR1)"))))	
                                         	 # "m4: NAA (2D~AR1)"))
    
    	 p <- ggplot(df.plot, aes(x=em.x, y=val)) +
                	  geom_boxplot(aes(fill=em2), outlier.shape = NA) +
                	  scale_fill_jco(name="", labels=lapply(levels(df.plot$em2), function(x) parse(text=x))) +
                	  coord_cartesian(ylim=c(-1,1)) +
                    xlab("Estimation model") +
    	              ylab(NULL) +
                	  geom_hline(yintercept = 0, linetype=2, color='black') +
    	              facet_grid(rows=vars(variable), cols=vars(om2), labeller = label_parsed, switch='y') +
                	  theme_bw() +
                	  theme(legend.position="bottom", strip.background.y = element_blank(), strip.placement = "outside", 
                	        strip.text.y = element_text(size = 12), strip.text.x = element_text(size = 8, margin = margin(3,1,1,1, "pt")),
                	        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8),
                	        legend.text = element_text(margin = margin(r = 6, l=1,unit = "pt"), hjust = 0, size=8), 
                	        legend.box.margin = margin(0,0,0,0), legend.margin = margin(0,0,0,0))
      title <- ggdraw() + draw_label("Operating model", hjust = 0.3, vjust=1) + theme(plot.margin = margin(0, 0, 0, 0))
      p1 <- plot_grid(title, p, ncol = 1, rel_heights = c(0.045, 1))

      png(file.path(plots_dir, paste0(stock.id,"_boxplots_",types[ty],".png")), height=8.5, width=6, units='in', res=300)
      print(p1)
      dev.off()
      # return(p1)
      # plot_grid(title, p, ncol = 1, rel_heights = c(0.045, 1))
    }
  } else { # individual plots -----------------------------------------------------------------------------
    for(ty in sim.types){
      df.plot <- filter(results, type==levels(results$type)[ty])

      # 1. SSB
      p <- ggplot(df.plot, aes(x=year, y=SSB.rel)) +
          stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
          stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Year") +
        ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_grid(rows=vars(em), cols=vars(om)) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("1_ssb_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()

      png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=SSB.rel)) +
        geom_boxplot(aes(fill=em), outlier.shape = NA) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Estimation model") +
        ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
        labs(title="Operating model") +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_wrap(vars(om), nrow=1) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()

      # 2. Fishing mortality
      p <- ggplot(df.plot, aes(x=year, y=F.rel)) +
          stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
          stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Year") +
        ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_grid(rows=vars(em), cols=vars(om)) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

      png(file.path(plots_dir,paste0("2_F_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()

      png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=F.rel)) +
        geom_boxplot(aes(fill=em), outlier.shape = NA) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Estimation model") +
        ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
        labs(title="Operating model") +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_wrap(vars(om), nrow=1) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()      

      # 3. SSB / SSB_40
      p <- ggplot(df.plot, aes(x=year, y=relB.rel)) +
          stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
          stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
        # geom_line(aes(group=sim),color='grey', alpha=.2) +
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Year") +
        ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_grid(rows=vars(em), cols=vars(om)) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

      png(file.path(plots_dir,paste0("3_relB_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()

      png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=relB.rel)) +
        geom_boxplot(aes(fill=em), outlier.shape = NA) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Estimation model") +
        ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
        labs(title="Operating model") +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_wrap(vars(om), nrow=1) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()   

      # 4. F / F40
      p <- ggplot(df.plot, aes(x=year, y=relF.rel)) +
          stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
          stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Year") +
        ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_grid(rows=vars(em), cols=vars(om)) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("4_relF_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()

      png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=relF.rel)) +
        geom_boxplot(aes(fill=em), outlier.shape = NA) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Estimation model") +
        ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
        labs(title="Operating model") +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_wrap(vars(om), nrow=1) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()    

      # 5. Catch
      p <- ggplot(df.plot, aes(x=year, y=catch.rel)) +
          stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
          stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Year") +
        ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_grid(rows=vars(em), cols=vars(om)) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("5_catch_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()

      # boxplots (collapse time series)
      png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=catch.rel)) +
        geom_boxplot(aes(fill=em), outlier.shape = NA) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Estimation model") +
        ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
        labs(title="Operating model") +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_wrap(vars(om), nrow=1) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()      

      # 6. Recruitment
      p <- ggplot(df.plot, aes(x=year, y=R.rel)) +
          stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
          stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Year") +
        ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_grid(rows=vars(em), cols=vars(om)) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("6_R_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()

      png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=R.rel)) +
        geom_boxplot(aes(fill=em), outlier.shape = NA) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(0,2)) +
        xlab("Estimation model") +
        ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
        labs(title="Operating model") +
        geom_hline(yintercept = 1, linetype=2, color='black') +
        facet_wrap(vars(om), nrow=1) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()           

# ---------------------------------------------------------------------------
      # now with TMB bias correction
      if(plot.eps){
        # 1. SSB
        p <- ggplot(df.plot, aes(x=year, y=SSB.rel.bc)) +
            stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
            stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Year") +
          ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_grid(rows=vars(em), cols=vars(om)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

        png(file.path(plots_dir,paste0("1_ssb_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=SSB.rel.bc)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Estimation model") +
          ylab(expression(SSB["sim fit"]~"/"~SSB["sim data"])) +
          labs(title="Operating model") +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_wrap(vars(om), nrow=1) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off() 

        # 2. Fishing mortality
        p <- ggplot(df.plot, aes(x=year, y=F.rel.bc)) +
            stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
            stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Year") +
          ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_grid(rows=vars(em), cols=vars(om)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

        png(file.path(plots_dir,paste0("2_F_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=F.rel.bc)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Estimation model") +
          ylab(expression(F["sim fit"]~"/"~F["sim data"])) +
          labs(title="Operating model") +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_wrap(vars(om), nrow=1) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off()     

        # 3. SSB / SSB40
        p <- ggplot(df.plot, aes(x=year, y=relB.rel.bc)) +
            stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
            stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
          # geom_line(aes(group=sim),color='grey', alpha=.2) +
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Year") +
          ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_grid(rows=vars(em), cols=vars(om)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

        png(file.path(plots_dir,paste0("3_relB_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=relB.rel.bc)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Estimation model") +
          ylab(expression(frac(B,B[40]["%"])~"(sim fit)"~"/"~frac(B,B[40]["%"])~"(sim data)")) +
          labs(title="Operating model") +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_wrap(vars(om), nrow=1) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off()   

        # 4. F / F40
        p <- ggplot(df.plot, aes(x=year, y=relF.rel.bc)) +
            stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
            stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Year") +
          ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_grid(rows=vars(em), cols=vars(om)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
        png(file.path(plots_dir,paste0("4_relF_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=relF.rel.bc)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Estimation model") +
          ylab(expression(frac(F,F[40]["%"])~"(sim fit)"~"/"~frac(F,F[40]["%"])~"(sim data)")) +
          labs(title="Operating model") +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_wrap(vars(om), nrow=1) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off()   

        # 5. Catch
        p <- ggplot(df.plot, aes(x=year, y=catch.rel.bc)) +
            stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
            stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Year") +
          ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_grid(rows=vars(em), cols=vars(om)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
        png(file.path(plots_dir,paste0("5_catch_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=catch.rel.bc)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Estimation model") +
          ylab(expression(Catch["sim fit"]~"/"~Catch["sim data"])) +
          labs(title="Operating model") +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_wrap(vars(om), nrow=1) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off() 

        # 6. Recruitment                 
        p <- ggplot(df.plot, aes(x=year, y=R.rel.bc)) +
            stat_flquantiles(probs=c(0.25, 0.75), alpha=0.5, fill="grey", geom="ribbon") + # middle 50%
            stat_flquantiles(probs=c(0.10, 0.90), alpha=0.35, fill="grey", geom="ribbon") + # middle 80%
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Year") +
          ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_grid(rows=vars(em), cols=vars(om)) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
        png(file.path(plots_dir,paste0("6_R_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=R.rel.bc)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(0,2)) +
          xlab("Estimation model") +
          ylab(expression(Recruitment["sim fit"]~"/"~Recruitment["sim data"])) +
          labs(title="Operating model") +
          geom_hline(yintercept = 1, linetype=2, color='black') +
          facet_wrap(vars(om), nrow=1) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off()         
      }
    }
  }
}
