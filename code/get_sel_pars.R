library(wham)
mod = readRDS("/home/bstock/Documents/ms/wham-sim/results/v2_fits/GBhaddock_sel/Sel-3.rds")

get_sel_pars <- function(mod){
	pars <- as.data.frame(matrix(summary(mod$sdrep)[rownames(summary(mod$sdrep))=="sel_repars",1:2], ncol=2))
	pars$low <- pars[,1] - 1.96*pars[,2]
	pars$hi <- pars[,1] + 1.96*pars[,2]
	rho_trans <- function(x){
		y <- 2/(1 + exp(-2*x)) - 1
		return(y)
	}
	pars[1,] = exp(pars[1,])
    pars[2,] <- rho_trans(pars[2,])
    pars[3,] <- rho_trans(pars[3,])
    rownames(pars) <- c("sigmaSel","phiPar","phiYear")
    return(round(pars,2))
}

get_sel_pars(mod)

