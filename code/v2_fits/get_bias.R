setwd("/home/bstock/Documents/ms/wham-sim")
bc.type = 2 # _oepe

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(tidyverse)
inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
rho_trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
source("/home/bstock/Documents/ms/wham-sim/code/v2_fits/get_results.R")

# summarize bias across stocks and processes (self-tests)
ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring", "SNEMAYT","butterfish","NScod","GBhaddock")
re = c(rep("NAA",5), rep("M",3),"sel")
em.less.complex <- em.more.complex <- self <- c()
for(j in 1:length(ids)){
	results <- get_results(stock.id=ids[j], re=re[j], bc.type=bc.type)
	df.plot <- filter(results, type==levels(results$type)[2]) %>%
		group_by(om2, em.x, em2, sim) %>%
		summarize(SSB.rel=median(SSB.rel), F.rel=median(F.rel), relB.rel=median(relB.rel, na.rm=T),
	        relF.rel=median(relF.rel,na.rm=T), R.rel=median(R.rel), .groups = 'drop_last') %>%
		select(om2, em2, em.x, SSB.rel, F.rel, relB.rel, relF.rel, R.rel) %>% as.data.frame
	df.plot <- do.call(data.frame, lapply(df.plot, function(x) replace(x, is.infinite(x), NA)))
	# df.plot[df.plot == 0] = NA
	df.plot <- df.plot %>%
	pivot_longer(-c(om2,em2,em.x), names_to = "variable", values_to = "val") %>%
	group_by(om2, em.x, em2, variable) %>%
	summarize(rel_err_var = var(val, na.rm=T), rel_err_mean = mean(val, na.rm=T)-1, K=sum(!is.na(val)), .groups = 'keep') %>%
	mutate(rel_err_se = sqrt(rel_err_var/K)) %>%
	mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se)       
	df.plot$variable <- factor(df.plot$variable, levels=c("SSB.rel", "F.rel", "relB.rel", "relF.rel", "R.rel"), 
	                   labels=c("SSB", "F", expression(B/B[40]["%"]), expression(F/F[40]["%"]), "Recruitment"))
	levels(df.plot$om2) = levels(df.plot$em2)

	df.plot$emcat <- "self"
	df.plot$emcat[as.numeric(df.plot$em2) > as.numeric(df.plot$om2)] <- "more"
	df.plot$emcat[as.numeric(df.plot$em2) < as.numeric(df.plot$om2)] <- "less"
	df <- df.plot %>% filter(emcat == "self")
	self <- c(self, df$rel_err_mean)
	df <- df.plot %>% filter(emcat == "more")
	em.more.complex <- c(em.more.complex, df$rel_err_mean)
	df <- df.plot %>% filter(emcat == "less")
	em.less.complex <- c(em.less.complex, df$rel_err_mean)
}

hist(self, col='grey', breaks=20)
abline(v=.02, lty=2, lwd=2)
abline(v=-.02, lty=2, lwd=2)
sum(abs(self) > .02)/length(self)

hist(em.more.complex, col='grey', breaks=20)
abline(v=.02, lty=2, lwd=2)
abline(v=-.02, lty=2, lwd=2)
sum(abs(em.more.complex) > .02)/length(em.more.complex)

hist(em.less.complex, col='grey', breaks=20)
abline(v=.02, lty=2, lwd=2)
abline(v=-.02, lty=2, lwd=2)
sum(abs(em.less.complex) > .02)/length(em.less.complex)

