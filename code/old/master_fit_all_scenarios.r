library(wham)
library(parallel)
library(foreach)
library(doParallel)
# convert Lowestoft input files to vanilla ASAP
myCluster <- makeCluster(50, # number of cores to use
  type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

rvals = paste0("r",1:36)

source("convert_ICES_to_ASAP.r")
for(r in rvals[17:36])
{
  print(r)
  write.dir = paste0("~/work/ICES/MGWG/stock_recruit/vpa/",r,"/")
  system.time(foreach(i = 1:300) %dopar% {
    #for(i in 1:1){
    user.od <- user.wd <- paste0("~/work/ICES/MGWG/stock_recruit/vpa/",r,"/iter", i, "/")
    model.id <- paste0("iter",i)
    print(user.wd)
    ICES2ASAP(user.wd, user.od, model.id = model.id, ices.id= '')
  })
  system.time(foreach(i = 1:300) %dopar% {
    x = wham::read_asap3_dat(paste0("~/work/ICES/MGWG/stock_recruit/vpa/",r,"/iter", i, "/ASAP_iter",i,".dat"))
    x = wham::prepare_wham_input(x)
    #x$data$agg_index_catch[] = 0.1  
    #x$data$agg_index_sigma[] = 0.2  
    x$data$index_Neff[] = 100
    x$data$catch_Neff[] = 200
    x$par$log_R_sigma[] = log(0.5) 
    x$data$fracyr_indices[] = 0
    x$data$Fbar_ages = 4:7
    x$random = "log_R"
    x$map = x$map[!(names(x$map) %in% c("mean_rec_pars","log_R_sigma"))]
    x$par$mean_rec_pars = 0
    x$data$random_recruitment = 1
    x$data$recruit_model = 2 #random about mean
    saveRDS(x, file = paste0("~/work/ICES/MGWG/stock_recruit/vpa/",r,"/iter", i, "/WHAM_input_iter",i,".RDS"))
  })
}
for(r in rvals[-1])
{
  print(r)
  source("fit_all_models.r")
}

for(r in rvals[7:9])
{
  print(r)
  source("make_RDS_smaller.r")
}

system.time(res <- foreach(i = 1:50, .combine = "rbind") %dopar% {
  fit = readRDS(file = paste0("~/work/ICES/MGWG/stock_recruit/vpa/",r,"/iter", i, "/WHAM_fit0_iter",i,".RDS"))
  x = summary(fit$sdrep)[rownames(summary(fit$sdrep)) == "logit_SR_h",][1,]
  #0.2 + 0.8/(1 + exp(-(x[1])))
  out = 0.2 + 0.8/(1 + exp(-(x[1])))
  out = c(out, 0.2 + 0.8/(1 + exp(-(x[1] + c(-1,1)*qnorm(0.975)*x[2]))))
  x = summary(fit$sdrep)[rownames(summary(fit$sdrep)) == "log_SR_R0",][1,]
  x[1] = x[1] + fit$rep$log_SPR0[1]
  out = c(out,  exp(x[1]), exp(x[1] + c(-1,1)*qnorm(0.975)*x[2]))
  x = summary(fit$sdrep)[rownames(summary(fit$sdrep)) == "log_SR_a",][1,]
  out = c(out,  exp(x[1]), exp(x[1] + c(-1,1)*qnorm(0.975)*x[2]))
  x = summary(fit$sdrep)[rownames(summary(fit$sdrep)) == "log_SR_b",][1,]
  out = c(out,  exp(x[1]), exp(x[1] + c(-1,1)*qnorm(0.975)*x[2]))
  x = summary(fit$sdrep)[rownames(summary(fit$sdrep)) == "log_R_sigma",]
  out = c(out,  exp(x[1]), exp(x[1] + c(-1,1)*qnorm(0.975)*x[2]))
  out
})
library(plotrix)
h = 0.63
S0 = 925.24
sigR = 0.3
cairo_pdf(paste0(r,"_res.pdf"), family = "Times", height = 7, width = 10)
par(mfrow = c(1,3), mar = c(2,5,2,1), oma = c(4,1,2,1))
plotCI(x = 1:50, y = res[,1], ui = res[,3], li = res[,2], xlab = "", ylab = "Steepness", cex.lab = 1.5)
abline(h = h, lwd = 2)
tmp = c(mean(res[,1])+c(-1,1)*qnorm(0.975)*sd(res[,1])/sqrt(NROW(res))) - h
title(bquote(paste("95% Bias CI =" ~ .(paste(round(tmp,2),collapse = " - ")))))

plotCI(x = 1:50, y = res[,4], ui = res[,6], li = res[,5], xlab = "", ylab = "S0", cex.lab = 1.5)
abline(h = S0, lwd = 2)
tmp = c(mean(res[,4])+c(-1,1)*qnorm(0.975)*sd(res[,4])/sqrt(NROW(res))) - S0
title(bquote(paste("95% Bias CI =" ~ .(paste(round(tmp,2),collapse = " - ")))))

plotCI(x = 1:50, y = res[,13], ui = res[,15], li = res[,14], xlab = "", ylab = "Sigma_R", cex.lab = 1.5)
abline(h = sigR, lwd = 2)
tmp = c(mean(res[,13])+c(-1,1)*qnorm(0.975)*sd(res[,13])/sqrt(NROW(res))) - sigR
title(bquote(paste("95% Bias CI =" ~ .(paste(round(tmp,2),collapse = " - ")))))
mtext(side = 1, "Iteration", outer = TRUE, line = 2.5)
title(r, outer = TRUE, cex = 1.5)
dev.off()


stopCluster(myCluster)
