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

# source("/home/bstock/Documents/ms/wham-sim/code/v3_final_figure_pdfs/fig7.R")

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(here)
library(tidyverse)
library(ggsci)
library(viridis)
library(cowplot)
library(gridGraphics)

# get results into data frame
n.mods = 5
res_dir <- here("results","bias_correct_oepe","SNEMAYT_Ecov2_oepe")
# plots_dir <- here("plots","bias_correct_oepe","SNEMAYT_Ecov2")
plots_dir = here("plots","v3","final_pdfs")
mod.list <- file.path(res_dir,paste0("m",1:n.mods,".rds"))
mods <- lapply(mod.list, readRDS)
# df.mods <- data.frame(Model = 1:n.mods, stringsAsFactors=FALSE)
# df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F, do.print=F)$tab)
# minAIC <- min(df.aic$AIC, na.rm=T)
# df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
# df.mods <- cbind(df.mods, df.aic)
# rownames(df.mods) <- NULL
#   Model dAIC     AIC
# 1     1 33.0 -1757.8
# 2     2 12.7 -1778.1
# 3     3 14.0 -1776.8
# 4     4  0.0 -1790.8
# 5     5  1.3 -1789.5

# ----------------------------------------------------
# Panel 1: CPI
mod = mods[[4]]
dat = mod$env$data
years <- seq(from=dat$year1_Ecov, by=1, length.out=dat$n_years_Ecov)
years_full <- seq(from=dat$year1_Ecov, by=1, length.out=dat$n_years_Ecov+dat$n_years_proj_Ecov)
ecov.pred = mod$rep$Ecov_x
ecov.obs = dat$Ecov_obs
sdrep = summary(mod$sdrep)
if("Ecov_obs_logsigma" %in% names(mod$env$par)){
ecov.obs.sig = matrix(exp(sdrep[rownames(sdrep) %in% "Ecov_obs_logsigma",1]), ncol=dat$n_Ecov) # all the same bc obs_sig_var --> 0
} else {
ecov.obs.sig = exp(mod$input$par$Ecov_obs_logsigma)
}
ecov.use = dat$Ecov_use_obs
ecov.obs.sig[ecov.use == 0] <- NA
ecov.pred.se = matrix(sdrep[rownames(sdrep) %in% "Ecov_x",2], ncol=dat$n_Ecov)
ecov.obs[ecov.use == 0] <- NA
ecov.res = (ecov.obs - ecov.pred[1:dat$n_years_Ecov,]) / ecov.obs.sig # standard residual (obs - pred)
df.cpi <- data.frame(Year=years_full, ecov.pred=ecov.pred[,1], 
						ecov.pred.low = ecov.pred[,1] - 1.96 * ecov.pred.se[,1],
						ecov.pred.high = ecov.pred[,1] + 1.96*ecov.pred.se[,1],
						ecov.low = ecov.obs[,1] - 1.96 * ecov.obs.sig[,1],
						ecov.high = ecov.obs[,1] + 1.96 * ecov.obs.sig[,1],
						ecov.obs = ecov.obs[,1])
# p1 <- ggplot(df.cpi, aes(x=Year, y=ecov.obs)) +
# 	geom_ribbon(aes(ymin=ecov.pred.low, ymax=ecov.pred.high), fill="dodgerblue", alpha=.35) +
# 	geom_line(aes(y=ecov.pred), color="dodgerblue", size=1) +
# 	geom_pointrange(aes(ymin=ecov.low, ymax=ecov.high)) +
# 	scale_x_continuous(breaks=seq(1975, 2015, by=10)) +
# 	ylab("Cold Pool Index (CPI)") +
# 	theme_bw()
vals <- approx(x = df.cpi$Year, y = df.cpi$ecov.pred, n = 5000)
df1 <- data.frame(x = vals$x, y = vals$y)
p1 <- ggplot(df.cpi, aes(x=Year, y=ecov.obs)) +
  geom_ribbon(aes(ymin=ecov.pred.low, ymax=ecov.pred.high), alpha=.2) +
  geom_point(data=df1, aes(x=x, y=y, color=y), size=.9) +
  geom_pointrange(aes(ymin=ecov.low, ymax=ecov.high)) +
  scale_x_continuous(breaks=seq(1975, 2015, by=10)) +
  scale_color_viridis_c(name="CPI", option = "plasma") +
  ylab("Cold Pool Index (CPI)") +
  theme_bw()

# ------------------------------------------------------
# panel 2: m1 SSB-Rec
# p2 <- plot.SR.pred.line(mods[[1]])
ssb.units = "mt"; recruits.units = "thousands"; scale.ssb = 1; scale.recruits = 1; age.recruit = 1
mod = mods[[1]]
std = summary(mod$sdrep)
ssb.ind = which(rownames(std) == "log_SSB")
years = mod$years
nyrs = length(years)
log.ssb <- std[ssb.ind,1]
R <- mod$rep$NAA[,1]
SR <- matrix(NA, (nyrs-age.recruit), 3)
SR[,1] <- years[1:(nyrs-age.recruit)]/scale.ssb
SR[,2] <- exp(log.ssb[1:(nyrs-age.recruit)])
SR[,3] <- R[age.recruit +1:(nyrs-age.recruit)]/scale.recruits
log_a <- mod$parList$mean_rec_pars[1]
log_b <- mod$parList$mean_rec_pars[2]
SR.par.year = nyrs
a.b.ind = which(rownames(std) == "mean_rec_pars")
l.ab = std[a.b.ind,1]
a.b.cov = mod$sdrep$cov.fixed[a.b.ind,a.b.ind]
lR.fn = function(la, lb, S) la  + log(S) - log(1 + exp(lb)*S)
dlR.dp = Deriv::Deriv(lR.fn, x = c("la","lb"))
seq.ssb <- seq(0, max(SR[,2]), length.out=300)
sd.pred.lR = sapply(seq.ssb, function(x)
{
  d = dlR.dp(l.ab[1], l.ab[2], S= x)
  return(c(sqrt(t(d) %*% a.b.cov %*% d)))
})
pred.lR <- lR.fn(l.ab[1], l.ab[2], seq.ssb)
ci.pred.lR = pred.lR + qnorm(0.975)*cbind(-sd.pred.lR, sd.pred.lR) - log(scale.recruits)

scale.ssb = 1000
scale.recruits = 1000
df2 <- data.frame(SSB=SR[,2]/scale.ssb, Recruitment=SR[,3]/scale.recruits)
df3 <- data.frame(pred.ssb=seq.ssb/scale.ssb, pred.R=exp(pred.lR)/scale.recruits, pred.R.low=exp(ci.pred.lR[,1])/scale.recruits, pred.R.high=exp(ci.pred.lR[,2])/scale.recruits)
# p2 <- ggplot(df3, aes(x=pred.ssb, y=pred.R)) +
# 	geom_ribbon(data=df3, aes(x=pred.ssb, ymin=pred.R.low, ymax=pred.R.high), alpha=.2) +
# 	geom_line(data=df3, aes(x=pred.ssb, y=pred.R), size=1) +
# 	geom_point(data=df2, aes(x=SSB, y=Recruitment)) +
# 	scale_x_continuous(expand=c(0.01,0)) +
# 	scale_y_continuous(expand=c(0.01,0)) +
#   annotate("text", x=1, y=160, label=expression(hat(R)~"="~italic(f)*phantom()*"(SSB)"), hjust=0, size = 4, parse = TRUE) +
# 	# geom_label(aes(x=pred.ssb, y=pred.R, label=lab), size=4, alpha=1, label.r=unit(0, "lines"), label.size=NA,
# 				# data=data.frame(pred.ssb=4, pred.R=160, lab="R ~ f(SSB)")) +
#   ylab(expression("Age-1"~"Recruits"~"(x"~10^3*phantom()*")")) +
# 	xlab(expression("SSB"~"(x"~10^3~"mt)")) +
# 	theme_bw()

# ---------------------------------------------------
# Panel 3
mod = mods[[4]]
std = summary(mod$sdrep)
ssb.ind = which(rownames(std) == "log_SSB")
years = mod$years
nyrs = length(years)
R <- mod$rep$NAA[age.recruit +1:(nyrs-age.recruit),1]

# predR_t is function of SSB_t, a, b_t
ssb <- exp(std[ssb.ind,1])[1:(nyrs-age.recruit)] # need to lag ssb
cpi <- mod$rep$Ecov_out[1:(nyrs-age.recruit),1] # modeled CPI in year t-1 affects R in year t (Ecov_out is lagged)
log_a <- mod$rep$log_SR_a # constant
log_b <- mod$rep$log_SR_b # varies by year according to CPI

seq.ssb <- seq(0, max(ssb), length.out=300)
lR.fn = function(la, lb, S) la  + log(S) - log(1 + exp(lb)*S)
cnames <- c("Year","CPI","SSB.seq","R.seq")
predR <- as.data.frame(matrix(NA, ncol = length(cnames), nrow = 0))
colnames(predR) <- cnames
for(y in 1:(nyrs-age.recruit)){
  tmp <- as.data.frame(matrix(NA, ncol = length(cnames), nrow = length(seq.ssb)))
  colnames(tmp) <- cnames
  tmp$SSB.seq = seq.ssb
  tmp$Year = years[y]
  tmp$CPI = cpi[y]
  tmp$R.seq = exp(lR.fn(log_a[y], log_b[y], seq.ssb))
  predR <- rbind(predR, tmp)
}

scale.ssb = 1000
scale.recruits = 1000
df4 <- data.frame(SSB=ssb/scale.ssb, Recruitment=R/scale.recruits, CPI=cpi)
df4 <- df4[seq(dim(df4)[1],1),]

# ymax = max(df3$pred.R, predR$R.seq/scale.recruits)
ymax = max(df3$pred.R.high)*1.02
xmax = max(df3$pred.ssb, predR$SSB.seq/scale.ssb)
p2 <- ggplot(df3, aes(x=pred.ssb, y=pred.R)) +
  geom_ribbon(data=df3, aes(x=pred.ssb, ymin=pred.R.low, ymax=pred.R.high), alpha=.2) +
  geom_line(data=df3, aes(x=pred.ssb, y=pred.R), size=1) +
  geom_point(data=df2, aes(x=SSB, y=Recruitment)) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,xmax)) +
  scale_y_continuous(expand=c(0.01,0), limits=c(0,ymax)) +
  annotate("text", x=1, y=ymax*.95, label=expression(hat(R)~"="~italic(f)*phantom()*"(SSB)"), hjust=0, size = 4, parse = TRUE) +
  ylab(expression("Age-1"~"Recruits"~"(x"~10^3*phantom()*")")) +
  xlab(expression("SSB"~"(x"~10^3~"mt)")) +
  theme_bw()
p3 <- ggplot(predR, aes(x=SSB.seq/scale.ssb, y=R.seq/scale.recruits, color=CPI, group=Year)) +
  geom_line(size=0.7) +
  # geom_point(data=df4, aes(x=SSB, y=Recruitment), inherit.aes = F) +
  geom_point(data=df4, aes(x=SSB, y=Recruitment, fill=CPI), inherit.aes = F, pch=21, size=2) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0,xmax)) +
  scale_y_continuous(expand=c(0.01,0), limits=c(0,ymax)) +
  scale_color_viridis_c(name="CPI", option = "plasma") +
  scale_fill_viridis_c(name="CPI", option = "plasma") +
  annotate("text", x=1, y=ymax*.95, label=expression(hat(R)~"="~italic(f)*phantom()*"(SSB,CPI)"), hjust=0, size = 4, parse = TRUE) +
  ylab(expression("Age-1"~"Recruits"~"(x"~10^3*phantom()*")")) +
  xlab(expression("SSB"~"(x"~10^3~"mt)")) +
  theme_bw() +
  theme(legend.position="none")

# first align the top-row plot (p1) with the left-most plot of the bottom row (p2)
plots <- align_plots(p1, p2, align = 'v', axis = 'l')
bottom_row <- plot_grid(plots[[2]], p3, labels = c('B', 'C'), label_size = 16, label_fontface = 'plain')

# png(file.path(plots_dir,"CPI_Recruit_SNEMAYT.png"), width=8, height=7.5, units='in', res=300)
cairo_pdf(file.path(plots_dir, "fig7.pdf"), width=8, height=7.5)
print(plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 16, label_fontface = 'plain', ncol = 1))
dev.off()




