# setwd("/home/bstock/Documents/ms/wham-sim")
# bc.type=2; sim.types=2
# plots_dir = file.path(getwd(),"plots","bias_correct_oepe")
# res_dir=file.path(getwd(),"results")
# simdata_dir=file.path(getwd(),"data","simdata")
# library(wham)
# library(tidyverse)

plot_aic_cross <- function(df.aic, plots_dir = file.path(getwd(),"plots","bias_correct_oepe"), bystock=TRUE){ 
  df <- df.aic %>% group_by(id, re, bc.type, sim.type, om, sim) %>%
    mutate(minAIC = min(aic, na.rm=T)) %>% ungroup()
  df$daic <- df$aic - df$minAIC
  df$best <- 0
  df$best[which(df$daic == 0)] <- 1
  re.labs = unique(df$re)
  n.re = length(re.labs)
  types <- c("OE","OEPE")
  mlabs = list(c(expression(paste("m1:")~paste("SCAA")~paste("(IID)")),
                      expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")),
                      expression(paste("m3:")~paste("NAA")~paste("(IID)")),
                      expression(paste("m4:")~paste("NAA")~paste("(2D AR1)"))),
               c(expression(paste("m1:")~paste("none")),
                      expression(paste("m2:")~paste("IID")),
                      expression(paste("m3:")~paste("2D")~paste("AR1"))),
               c(expression(paste("m1:")~paste("none")),
                 expression(paste("m2:")~paste("IID")),
                 expression(paste("m3:")~paste("2D")~paste("AR1"))))
  mlabs_m = list(c(expression(paste("SCAA")~paste("(IID)")),
                 expression(paste("SCAA (")*AR1[y]*paste(")")),
                 expression(paste("NAA")~paste("(IID)")),
                 expression(paste("NAA")~paste("(2D AR1)"))),
               c(expression(paste("none")),
                 expression(paste("IID")),
                 expression(paste("2D")~paste("AR1"))),
               c(expression(paste("none")),
                 expression(paste("IID")),
                 expression(paste("2D")~paste("AR1"))))  
  
  if(!bystock) plots <- plots_m <- list()
  for(re in 1:n.re){
    df.plot <- df[df$re==re.labs[re],]
    df.plot2 <- df.plot %>% group_by(id, bc.type, sim.type, om, em) %>% summarize(n.selected = sum(best)) %>% ungroup()
    # for(ty in unique(df.plot2$sim.type)){
    ty=2
      df.plot3 <- df.plot2[df.plot2$sim.type == ty,]
      # df.plot3$om2 <- factor(df.plot3$om, levels=1:max(df.plot3$om), labels=paste0("m",1:max(df.plot3$om)))
      # df.plot3$em2 <- factor(df.plot3$em, levels=max(df.plot3$em):1, labels=paste0("m",max(df.plot3$em):1))
      df.plot3$om2 <- factor(df.plot3$om, levels=1:max(df.plot3$om), labels=mlabs[[re]])
      df.plot3$em2 <- factor(df.plot3$em, levels=max(df.plot3$em):1, labels=rev(mlabs[[re]]))
      df.plot3$om3 <- factor(df.plot3$om, levels=1:max(df.plot3$om), labels=mlabs_m[[re]])
      df.plot3$em3 <- factor(df.plot3$em, levels=max(df.plot3$em):1, labels=rev(mlabs_m[[re]]))
      id.labs <- unique(df.plot3$id)
      n.ids <- length(id.labs)
      if(bystock){
        for(id in 1:n.ids){
          plots_dir_i <- file.path(plots_dir, paste0(id.labs[id],"_",re.labs[re]))
          df.plot4 <- df.plot3[df.plot3$id == id.labs[id],] %>% group_by(om2, em2) %>% summarize(p.selected=sum(n.selected)/100, .groups='keep')
          png(file.path(plots_dir_i, paste0("9_aic_cross_",id.labs[id],"_",re.labs[re],"_",types[ty],".png")), res=300, units='in', height=7, width=7)
          print(ggplot(data = df.plot4,
                 mapping = aes(x = om2,
                               y = em2)) +
            geom_tile(aes(fill = log(p.selected))) +
            # geom_text(aes(label = sprintf("%1.0f", n.selected)), vjust = 1) +
            geom_text(aes(label = sprintf('%#.2f', p.selected)), size=5, vjust = 1) +
            scale_fill_distiller(palette="YlOrRd", direction=1, na.value = "white") +
            # scale_x_discrete(expand = c(0,0), position='top') +
            # scale_y_discrete(expand = c(0,0)) +
            scale_x_discrete(expand = c(0,0), position='top', labels=function(l) parse(text=l)) +
            scale_y_discrete(expand = c(0,0), labels=function(l) parse(text=l), guide = guide_axis(angle = 90)) +
            xlab("Operating model") +
            ylab("Estimation model") +
            theme_bw() +
            theme(axis.text = element_text(size=12), axis.title = element_text(size=14),
                  legend.position='none'))
          dev.off()             
        }
      } else {
        df.plot4 <- df.plot3 %>% group_by(om2, em2) %>% summarize(p.selected=sum(n.selected)/(100*n.ids), .groups='keep')
        df.plot4_m <- df.plot3 %>% group_by(om3, em3) %>% summarize(p.selected=sum(n.selected)/(100*n.ids), .groups='keep')
        plots[[re]] <- ggplot(data = df.plot4,
                              mapping = aes(x = om2,
                                            y = em2)) +
          geom_tile(aes(fill = log(p.selected))) +
          # geom_text(aes(label = sprintf("%1.0f", n.selected)), vjust = 1) +
          geom_text(aes(label = sprintf('%#.3f', p.selected)), size=5, vjust = 1) +
          scale_fill_distiller(palette="YlOrRd", direction=1, na.value = "white") +
          # scale_x_discrete(expand = c(0,0), position='top') +
          # scale_y_discrete(expand = c(0,0)) +
          scale_x_discrete(expand = c(0,0), position='top', labels=function(l) parse(text=l)) +
          scale_y_discrete(expand = c(0,0), labels=function(l) parse(text=l), guide = guide_axis(angle = 90)) +
          xlab("Operating model") +
          ylab("Estimation model") +
          theme_bw() +
          theme(axis.text = element_text(size=12), axis.title = element_text(size=14),
                legend.position='none')
        
        plots_m[[re]] <- ggplot(data = df.plot4_m, # diff text sizes for multipanel
                     mapping = aes(x = om3,
                                   y = em3)) +
                geom_tile(aes(fill = log(p.selected))) +
                # geom_text(aes(label = sprintf("%1.0f", n.selected)), vjust = 1) +
                geom_text(aes(label = sprintf('%#.3f', p.selected)), size=4, vjust = 1) +
                scale_fill_distiller(palette="YlOrRd", direction=1, na.value = "white") +
                # scale_x_discrete(expand = c(0,0), position='top') +
                # scale_y_discrete(expand = c(0,0)) +
                scale_x_discrete(expand = c(0,0), position='top', labels=function(l) parse(text=l)) +
                scale_y_discrete(expand = c(0,0), labels=function(l) parse(text=l), guide = guide_axis(angle = 90)) +
                xlab("Operating model") +
                ylab("Estimation model") +
                theme_bw() +
                theme(axis.text = element_text(size=10), axis.title = element_text(size=14),
                      legend.position='none')
        png(file.path(plots_dir, paste0("aic_cross_",re.labs[re],"_",types[ty],".png")), res=300, units='in', height=7, width=7)
        print(plots[[re]])
        dev.off()             
      }
    # }
  }
  # multipanel aggregated across stocks
  if(!bystock){
    plots_m[[1]] <- plots_m[[1]] + theme(axis.text = element_text(size=8))
    png(file.path(plots_dir,"aic_cross_multipanel.png"), width=4, height=10, units='in', res=300)
    print(plot_grid(plotlist=plots_m, labels = LETTERS[1:n.re], label_size = 18, label_fontface = 'plain', ncol = 1))
    dev.off()
  }
}
  