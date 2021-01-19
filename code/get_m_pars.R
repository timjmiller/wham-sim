library(wham)
mods3 <- list(butterfish = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/butterfish_M_oepe/m3.rds"),
			SNEMAYT = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/SNEMAYT_M_oepe/m3.rds"))

mods2 <- list(butterfish = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/butterfish_M_oepe/m2.rds"),
			NScod = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/NScod_M_oepe/m2.rds"),
			SNEMAYT = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/SNEMAYT_M_oepe/m2.rds"))

get_M_pars <- function(mod){
	sigM <- matrix(summary(mod$sdrep)[rownames(summary(mod$sdrep))=="sigma_M",1:2], ncol=2)
	Mrho_y <- matrix(summary(mod$sdrep)[rownames(summary(mod$sdrep))=="rho_M_y",1:2], ncol=2)
	Mrho_a <- matrix(summary(mod$sdrep)[rownames(summary(mod$sdrep))=="rho_M_a",1:2], ncol=2)
	pars <- as.data.frame(rbind(sigM, Mrho_y, Mrho_a))
	pars$low <- pars[,1] - 1.96*pars[,2]
	pars$hi <- pars[,1] + 1.96*pars[,2]
	return(round(pars,2))
}

allpars2 <- lapply(mods2, get_M_pars)
allpars2
# $butterfish
#     V1   V2  low   hi
# 1 0.33 0.06 0.21 0.44
# 2 0.00 0.00 0.00 0.00
# 3 0.00 0.00 0.00 0.00

# $NScod
#     V1   V2  low   hi
# 1 0.52 0.05 0.43 0.62
# 2 0.00 0.00 0.00 0.00
# 3 0.00 0.00 0.00 0.00

# $SNEMAYT
#     V1  V2  low  hi
# 1 1.21 0.1 1.02 1.4
# 2 0.00 0.0 0.00 0.0
# 3 0.00 0.0 0.00 0.0

allpars3 <- lapply(mods3, get_M_pars)
allpars3
# $butterfish
#      V1   V2   low   hi
# 1  0.16 0.03  0.10 0.22
# 2  0.97 0.02  0.92 1.01
# 3 -0.38 0.30 -0.97 0.20

# $SNEMAYT
#     V1   V2  low   hi
# 1 0.79 0.14 0.51 1.06
# 2 0.63 0.16 0.31 0.94
# 3 0.40 0.16 0.09 0.70

