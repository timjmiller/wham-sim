# SNEMAYT NAA test results


id = "SNEMAYT_NAA"
res_dir <- file.path(getwd(),"results",id)
res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
res.list <- lapply(res.files, readRDS)
flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
results <- do.call(rbind, flatten.nested.list(res.list)) %>% as.data.frame
results <- sapply(results, as.numeric)
results <- as.data.frame(results)
types <- c("OE","OEPE")
mlabs = c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)")
tylabs = c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)")
results$om <- factor(results$om, levels=1:4, labels=mlabs)
results$em <- factor(results$em, levels=1:4, labels=mlabs)
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

simdata <- lapply(1:4, function(x) readRDS(file.path(getwd(),"data","simdata",id,paste0("simdata_om",x,".rds"))))
results$R.sim = NA
for(om in 1:4){
  for(em in 1:4){
    for(i in 1:100){
      for(ty in 1:2){
        res.ind <- which(results$om == mlabs[om] & results$em == mlabs[em] & results$sim == i & results$ty == tylabs[ty])
        results$R.sim[res.ind] <- simdata[[om]][[i]][[ty]]$NAA[,1]
      }
    }
  }
}
results$R.rel <- results$NAA1 / results$R.sim
results$R.rel.bc <- results$NAA1_bc / results$R.sim

ty=2
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
            	        legend.text = element_text(margin = margin(r = 6, l=1,unit = "pt"), hjust = 0, size=8), legend.box.margin = margin(0,0,0,0), legend.margin = margin(0,0,0,0))
  title <- ggdraw() + draw_label("Operating model", hjust = 0.3, vjust=1) + theme(plot.margin = margin(0, 0, 0, 0))
  plot_grid(title, p, ncol = 1, rel_heights = c(0.045, 1))