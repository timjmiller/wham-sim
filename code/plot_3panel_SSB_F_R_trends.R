# 3-panel plot of SSB, F, and R for each stock and process
#   - base/m1
#   - best model by AIC (usually 2D AR1)

# # testing
# library(wham)
# library(tidyverse)
# library(ggplotFL)
# library(ggsci)
# library(cowplot)
# library(data.table)
# ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring", "SNEMAYT","butterfish","NScod","GBhaddock","SNEMAYT")
# re = c(rep("NAA",5), rep("M",3),"sel","Ecov2")
# bc.type = 2 # _oepe
# res_dir=file.path("/home/bstock/Documents/ms/wham-sim/results")
# simdata_dir=file.path("/home/bstock/Documents/ms/wham-sim/data/simdata")
# plots_dir=file.path("/home/bstock/Documents/ms/wham-sim/plots")

plot_3panel_SSB_F_R_trends <- function(ids, re, bc.type=2, plot.eps=FALSE,
                      res_dir=file.path(getwd(),"results"), 
                      simdata_dir=file.path(getwd(),"data","simdata"),
                      plots_dir=file.path(getwd(),"plots")){ 
  ind.naa <- which(re == "NAA")
  if(length(ind.naa) == 0) id.naa <- NULL else id.naa <- paste0(ids[ind.naa],"_NAA")
  ind.m <- which(re == "M")
  if(length(ind.m) == 0) id.m <- NULL else id.m <- paste0(ids[ind.m],"_M")
  ind.ecov <- which(re == "Ecov2")
  if(length(ind.ecov) == 0) id.ecov <- NULL else id.ecov <- paste0(ids[ind.ecov],"_Ecov2")
  ind.sel <- which(re == "sel")
  if(length(ind.sel) == 0) id.sel <- NULL else id.sel <- paste0(ids[ind.sel],"_sel")
  id.list <- list(NAA = id.naa, M = id.m, Ecov2=id.ecov, sel=id.sel)
  
  for(ni in 1:length(id.list)){
    n.stocks <- length(id.list[[ni]])
    if(n.stocks > 0){
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
        mlabs_expr = c(expression(paste("m1:")~paste("RW-none")), 
                       expression(paste("m2:")~paste("RW-linear")),
                       expression(paste("m3:")~paste("RW-poly")),
                       expression(paste("m4:")~paste("AR1-linear")),
                       expression(paste("m5:")~paste("AR1-poly")))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1:n.mods)
      }
      if(names(id.list)[ni] == "sel"){
        n.mods=3
        mlabs = c("m1: none","m2: IID","m3: 2D AR1")
        mlabs_expr = c(expression(paste("m1:")~paste("none")), 
                       expression(paste("m2:")~paste("IID")), 
                       expression(paste("m3:")~paste("2D")~paste("AR1")))
        mlabs_short <- mlabs
        names(mlabs_short) = paste0("m",1:3)
      }
      
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
        plots_dir_i <- file.path(plots_dir, bc)
        df.aic <- readRDS(file.path(plots_dir_i, "daic.rds"))
        st.lab <- strsplit(id.list[[ni]][st],"_")[[1]][1]
        re.lab <- strsplit(id.list[[ni]][st],"_")[[1]][2]
        best <- df.aic %>% filter(Stock == st.lab & re == re.lab & daic == 0) 
        mod.list <- c(file.path(res_dir_i, "m1.rds"), file.path(res_dir_i, paste0("m",best$m,".rds")))
        mods <- lapply(mod.list, readRDS)
        mod.labs <- mlabs[c(1,best$m)]

        # plot
        # years = mods[[1]]$years_full # include projections
        years = mods[[1]]$years # remove projections
        alpha = .05 # 95% CI
        df <- data.frame(matrix(NA, nrow=0, ncol=6))
        colnames(df) <- c("Year","var","val","lo","hi","Model")
        for(i in 1:length(mods)){
          # get SSB + 95% CI
          std = summary(mods[[i]]$sdrep)
          ssb.ind <- which(rownames(std) == "log_SSB")
          log.ssb <- std[ssb.ind,1]
          ssb = exp(log.ssb)/1000
          ssb.cv <- std[ssb.ind,2]
          log.ssb.ci <- log.ssb + cbind(qnorm(1-alpha/2)*ssb.cv, -qnorm(1-alpha/2)*ssb.cv)
          ssb.ci = exp(log.ssb.ci)/1000
          df <- rbind(df, data.frame(Year=years, var="SSB", val=ssb, lo=ssb.ci[,2], hi=ssb.ci[,1], Model=mod.labs[i]))
          
          # get F + 95% CI
          n_ages = mods[[i]]$env$data$n_ages
          faa.ind <- which(rownames(std) == "log_FAA_tot")
          log.faa <- matrix(std[faa.ind,1], length(years), n_ages)
          faa.cv <- matrix(std[faa.ind,2], length(years), n_ages)
          age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
          full.f.ind = cbind(1:length(years), age.full.f)
          log.full.f <- log.faa[full.f.ind]
          full.f.cv <- faa.cv[full.f.ind]
          log.f.ci <- log.full.f + cbind(qnorm(1-alpha/2)*full.f.cv, -qnorm(1-alpha/2)*full.f.cv)
          full.f = exp(log.full.f)
          df <- rbind(df, data.frame(Year=years, var="F", val=full.f, lo=exp(log.f.ci[,2]), hi=exp(log.f.ci[,1]), Model=mod.labs[i]))
          
          # get Recruitment + 95% CI
          scale.R = 10000 # divide by 10^4
          ind = rownames(std) == "log_NAA_rep"
          templo = exp(array(std[ind,1] - qnorm(0.975)*std[ind,2], dim = c(length(years), n_ages)))/scale.R
          temphi = exp(array(std[ind,1] + qnorm(0.975)*std[ind,2], dim = c(length(years), n_ages)))/scale.R
          temp = exp(array(std[ind,1], dim = c(length(years), n_ages)))/scale.R
          df <- rbind(df, data.frame(Year=years, var="Recruitment (x 10000)", val=temp[,1], lo=templo[,1], hi=temphi[,1], Model=mod.labs[i]))
        }
        df$Model <- factor(df$Model, levels=mod.labs, labels=mod.labs)
        df$Year <- as.integer(df$Year)
        dat <- data.table(df)
        
        # force y axes to 0
        dat[,y_min := 0, by = var]
        # trim large CIs to zoom to 120% of max MLE
        dat <- dat %>% group_by(var) %>% mutate(y_max = 1.2*max(val)) %>% as.data.frame
        dat$hi[dat$hi > dat$y_max] = dat$y_max[dat$hi > dat$y_max]
        
        # # exclude years before 1993 in dat bc ggplot doesn't recalculate y-axis limits
        # dat <- subset(dat, Year > 1992)
        
        png(file.path(plots_dir_i, paste0("3panel_SSB_F_R_trends_",id.list[[ni]][st],".png")), units='in', res=300, width=5.5, height=7.5)
        print(ggplot(dat, aes(x=Year, y=val, color=Model, group=Model)) +
          geom_ribbon(aes(ymin=lo, ymax=hi, fill=Model), color=NA, alpha=.15) +
          geom_line(size=.8) +
          # geom_vline(xintercept = tail(mods[[1]]$years,1), linetype=2, size=.4) +
          facet_wrap(vars(var), scales="free_y", ncol=1, strip.position = "left") +
          ylab(NULL) +
          geom_blank(aes(y = y_min)) +
          scale_y_continuous(expand=c(0.01,0.01)) +
          scale_x_continuous(expand=c(0.01,0.01)) +
          # coord_cartesian(xlim=c(1995,2021)) +
          scale_color_jco() +
          scale_fill_jco() +
          theme_bw() +
          theme(strip.background = element_blank(), strip.placement = "outside", 
                legend.position="top", legend.box.margin = margin(0,0,0,0), legend.margin = margin(0,0,0,0)))
        dev.off()
      }
    }
  }
}
