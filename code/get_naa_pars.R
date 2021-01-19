library(wham)
# mods <- list(ICEherring = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/ICEherring_NAA_oepe/m4.rds"),
# 			GBhaddock = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/GBhaddock_NAA_oepe/m4.rds"),
# 			butterfish = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/butterfish_NAA_oepe/m4.rds"),
# 			NScod = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/NScod_NAA_oepe/m4.rds"),
# 			SNEMAYT = readRDS("/home/bstock/Documents/ms/wham-sim/results/bias_correct_oepe/SNEMAYT_NAA_oepe/m4.rds"))
mods <- list(ICEherring = readRDS("/home/bstock/Documents/ms/wham-sim/results/v2_fits/ICEherring_NAA/NAA-4.rds"),
			GBhaddock = readRDS("/home/bstock/Documents/ms/wham-sim/results/v2_fits/GBhaddock_NAA/NAA-4.rds"),
			butterfish = readRDS("/home/bstock/Documents/ms/wham-sim/results/v2_fits/butterfish_NAA/NAA-4.rds"),
			NScod = readRDS("/home/bstock/Documents/ms/wham-sim/results/v2_fits/NScod_NAA/NAA-4.rds"),
			SNEMAYT = readRDS("/home/bstock/Documents/ms/wham-sim/results/v2_fits/SNEMAYT_NAA/NAA-4.rds"))

get_NAA_pars <- function(mod){
	naasig <- summary(mod$sdrep)[rownames(summary(mod$sdrep))=="NAA_sigma",1:2]
	naarho_y <- summary(mod$sdrep)[rownames(summary(mod$sdrep))=="NAA_rho_y",1:2]
	naarho_a <- summary(mod$sdrep)[rownames(summary(mod$sdrep))=="NAA_rho_a",1:2]

	pars <- as.data.frame(rbind(naasig, naarho_y, naarho_a))
	pars$low <- pars[,1] - 1.96*pars[,2]
	pars$hi <- pars[,1] + 1.96*pars[,2]
	return(round(pars,2))
}

allpars <- lapply(mods, get_NAA_pars)
allpars

# $ICEherring
#             Estimate Std. Error   low   hi
# NAA_sigma       0.61       0.09  0.42 0.79
# NAA_sigma.1     0.32       0.03  0.27 0.37
# naarho_y        0.17       0.11 -0.04 0.39
# naarho_a        0.48       0.10  0.28 0.68

# $GBhaddock
#             Estimate Std. Error   low   hi
# NAA_sigma       1.64       0.14  1.38 1.91
# NAA_sigma.1     0.32       0.03  0.26 0.37
# naarho_y        0.73       0.07  0.58 0.87
# naarho_a       -0.12       0.10 -0.32 0.08

# $butterfish
#             Estimate Std. Error   low   hi
# NAA_sigma       0.13       0.06  0.01 0.25
# NAA_sigma.1     0.25       0.05  0.15 0.36
# naarho_y        0.89       0.05  0.79 1.00
# naarho_a       -0.13       0.26 -0.63 0.37

# $NScod
#             Estimate Std. Error   low   hi
# NAA_sigma       0.91       0.10  0.72 1.10
# NAA_sigma.1     0.23       0.02  0.20 0.27
# naarho_y        0.40       0.09  0.23 0.57
# naarho_a        0.12       0.11 -0.09 0.33

# $SNEMAYT
#             Estimate Std. Error  low   hi
# NAA_sigma       0.73       0.11 0.53 0.94
# NAA_sigma.1     0.47       0.05 0.37 0.57
# naarho_y        0.53       0.09 0.35 0.71
# naarho_a        0.33       0.12 0.10 0.56

