# Brian Stock
# June 15 2020
# Simulation test WHAM
# Step 4: Collect simulation fits
#   NEFSC server version 
#   separate results files for each om/em cross test, for XX/YY in 1-4
#     results_omXX_emYY.rds
#     sdreps_omXX_emYY.rds
#     reps_omXX_emYY.rds
#   assume they are copied to /home/bstock/Documents/ms/wham-sim/results/SNEMAYT/NAA

# source("/home/bstock/Documents/ms/wham-sim/code/SNEMAYT_NAA/4_results_SNEMAYT_NAA.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggplotFL)
library(ggsci)

# get results into data frame
res_dir <- here("results","old-incorrect-bias-correction","SNEMAYT_NAA")
plots_dir <- here("plots","SNEMAYT_NAA")
res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
res.list <- lapply(res.files, readRDS)

sdrep.files <- list.files(path=res_dir, pattern = "sdreps", full.names = TRUE)
sdrep.list <- lapply(sdrep.files, readRDS)
rep.files <- list.files(path=res_dir, pattern = "reps", full.names = TRUE)
rep.list <- lapply(rep.files, readRDS)
simdata <- lapply(1:4, function(x) readRDS(here("data","simdata","SNEMAYT_NAA",paste0("simdata_om",x,".rds"))))

for(m in 1:16){
	for(ty in 1:2){
		for(i in 1:100){
			# print(paste0("m: ",m,", ty: ",ty,", sim:",i))
			if(class(sdrep.list[[m]][[ty]][[i]]) != 'character'){
				ind.naa <- which(rownames(sdrep.list[[m]][[ty]][[i]]) == "log_NAA_rep")
				n.yrs <- dim(res.list[[m]][[ty]][[i]])[1]
				n.ages <- sum(rownames(sdrep.list[[m]][[ty]][[i]])=="log_N1_pars")
	  			naa <- matrix(sdrep.list[[m]][[ty]][[i]][ind.naa,1], n.yrs, n.ages)
	  			logR <- naa[,1]
				logSPR <- rep.list[[m]][[ty]][[i]]$log_SPR_FXSPR
				B40.new <- exp(logR + logSPR)

				ind.ssb <- which(rownames(sdrep.list[[m]][[ty]][[i]]) == "log_SSB")
		  		B <- exp(sdrep.list[[m]][[ty]][[i]][ind.ssb,1])
				relB.fit.new <- B/B40.new

				om <- floor((m-1)/4)+1
				dat <- simdata[[om]][[i]][[ty]]
				logR.sim <- log(dat$NAA[,1])
				logSPR.sim <- dat$log_SPR_FXSPR
				B40.new.sim <- exp(logR.sim + logSPR.sim)
				B.sim <- dat$SSB
				relB.sim.new <- B.sim / B40.new.sim

				B.sim2 <- exp(dat$log_SSB)
				relB.sim.new2 <- B.sim2 / B40.new.sim

				ind.b40 <- which(rownames(sdrep.list[[m]][[ty]][[i]]) == "log_SSB_FXSPR")
		  		B40.fit.old <- exp(sdrep.list[[m]][[ty]][[i]][ind.b40,1])



				relB.rel.new <- relB.fit.new / relB.sim.new
				relB.rel.new2 <- relB.fit.new / relB.sim.new2
				res.list[[m]][[ty]][[i]] <- cbind(res.list[[m]][[ty]][[i]], relB.rel.new, relB.rel.new2)
				rownames(res.list[[m]][[ty]][[i]]) <- NULL
			}
		}
	}
}

flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
results <- do.call(rbind, flatten.nested.list(res.list)) %>% as.data.frame
results <- sapply(results, as.numeric)
results <- as.data.frame(results)

types <- c("OE","OEPE")
results$om <- factor(results$om, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
results$em <- factor(results$em, levels=1:4, labels=c("m1: SCAA (iid)","m2: SCAA (AR1_y)","m3: NAA (iid)","m4: NAA (2D AR1)"))
results$type <- factor(results$type, levels=1:2, labels=c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)"))
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

plots_dir <- "/home/bstock/Documents/ms/wham-sim/plots/SNEMAYT_NAA"
ty=2
	df.plot <- filter(results, type==levels(results$type)[ty]) %>%
              select(om, em, em.x, relB.rel, relB.rel.new) %>%
              pivot_longer(-c(om,em,em.x), names_to = "variable", values_to = "val") %>%
	            group_by(om, em)
	df.plot$val = df.plot$val - 1 # relative error
	df.plot$om2 <- factor(df.plot$om, labels=c(expression(paste("m1:")~paste("SCAA")~paste("(iid)")), 
	                                     expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
            	                         expression(paste("m3:")~paste("NAA")~paste("(iid)")), 
	                                     expression(paste("m4:")~paste("NAA")~paste("(2D AR1)"))))
	df.plot$em2 <- factor(df.plot$em, labels=c(expression(paste("m1:")~paste("SCAA")~paste("(iid)")), 
	                                     expression(paste("m2:")~paste("SCAA (")*AR1[y]*paste(")")), 
            	                         expression(paste("m3:")~paste("NAA")~paste("(iid)")), 
	                                     expression(paste("m4:")~paste("NAA")~paste("(2D AR1)"))))	
                                     	 # "m4: NAA (2D~AR1)"))

	library(cowplot)
  png(file.path(plots_dir,"SNEMAYT_NAA_recalc_relB.png"),height=4, width=6,res=300,units='in')
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
  plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
  dev.off()
