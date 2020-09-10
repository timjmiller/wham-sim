plot_conv <- function(df.conv, plots_dir){
  for(bc.ty in unique(df.conv$bc.type)){
    df <- df.conv %>% filter(bc.type == bc.ty)
    re.lab <- unique(df$re)
    for(re in 1:length(re.lab)){
      df.plot <- df[df$re == re.lab[re],]
      n.mods <- length(unique(df.plot$om))
      
      if(re.lab[re] == "NAA" & n.mods == 4){
        mlabs = c("m1: SCAA (IID)","m2: SCAA (AR1_y)","m3: NAA (IID)","m4: NAA (2D AR1)")
        mlabs_expr = c(expression(paste("m1:")~paste("SCAA")~paste("(IID)")), 
                       expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
                       expression(paste("m3:")~paste("NAA")~paste("(IID)")), 
                       expression(paste("m4:")~paste("NAA")~paste("(2D AR1)")))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1:4)
      }
      if(re.lab[re] == "M" & n.mods == 3){
        mlabs = c("m1: none","m2: IID","m3: 2D AR1")
        mlabs_expr = c(expression(paste("m1:")~paste("none")), 
                       expression(paste("m2:")~paste("IID")), 
                       expression(paste("m3:")~paste("2D")~paste("AR1")))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1:3)
      }
      if(re.lab[re] == "M" & n.mods == 2){ # NScod no m3
        mlabs = c("m1: none","m2: IID")
        mlabs_expr = c(expression(paste("m1:")~paste("none")), 
                       expression(paste("m2:")~paste("IID")))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1:2)
      }
      df.plot$om <- factor(df.plot$om, levels=1:n.mods, labels=mlabs)
      df.plot$em <- factor(df.plot$em, levels=1:n.mods, labels=mlabs)
      df.plot$em.x <- fct_recode(df.plot$em, !!!mlabs_short)
      df.plot$om2 <- factor(df.plot$om, labels=mlabs_expr)
      df.plot$em2 <- factor(df.plot$em, labels=mlabs_expr)
      
      # ----------------------------------------------------------
      p <- ggplot(df.plot, aes(x=id, y=p.conv, fill=em)) +
        geom_point(aes(shape=sim.type), size=2, position=position_dodge(width=1)) +
        xlab("Stock") +
        ylab("Proportion converged") +
        coord_cartesian(ylim=c(.6,1)) +
        # scale_fill_jco(name="Estimation model", labels = scales::parse_format(), guide = guide_legend(override.aes = list(shape=21))) +
        scale_fill_jco(name="Estimation model", guide = guide_legend(override.aes = list(shape=21))) +
        scale_shape_manual(name="Simulation type",values=c(21,24)) +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              strip.text.x = element_text(size = 10))
      title <- ggdraw() + draw_label("Operating model", hjust = 1, vjust=1) + theme(plot.margin = margin(0, 0, 0, 0))
      p1 <- plot_grid(title, p, ncol = 1, rel_heights = c(0.12, 1))
      
      png(file.path(plots_dir, paste0("prop_converged_",re.lab[re],".png")), units='in', res=300, width=9.5, height = 3.5)
      print(p1)
      dev.off()
      
      # plot one sim type at a time
      simtypes <- unique(df.plot$sim.type)
      for(ty in 1:length(simtypes)){
        df.plot2 <- filter(df.plot, sim.type == simtypes[ty])
        
        p <- ggplot(df.plot2, aes(x=id, y=p.conv, fill=em)) +
          geom_point(shape=21, size=2, position=position_dodge(width=1)) +
          coord_cartesian(ylim=c(.6,1)) +
          xlab("Stock") +
          ylab("Proportion converged") +
          # scale_fill_jco(name="Estimation model", labels = scales::parse_format(), guide = guide_legend(override.aes = list(shape=21))) +
          scale_fill_jco(name="Estimation model") +
          facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                strip.text.x = element_text(size = 10))
        title <- ggdraw() + draw_label("Operating model", hjust = 1, vjust=1) + theme(plot.margin = margin(0, 0, 0, 0))
        p1 <- plot_grid(title, p, ncol = 1, rel_heights = c(0.12, 1))
        
        png(file.path(plots_dir, paste0("prop_converged_",re.lab[re],"_",simtypes[ty],".png")), units='in', res=300, width=9.5, height = 3.5)
        print(p1)
        dev.off() 
      }
    }
  }
}
