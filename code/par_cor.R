# Look at correlation between fixed effect pars that constrain random effects
#   e.g. sigma_R, sigma_a, rho_year, and rho_age for NAA

get_cor <- function(mod, parnames){
	std <- summary(mod$sdrep)
 	cov <- mod$sdrep$cov.fixed
	ind <- which(rownames(cov) %in% parnames)
	naa.cov <- cov[ind, ind]
	naa.cor <- cov2cor(naa.cov)
	return(naa.cor)
}

library(wham)
library(here)
library(tidyverse)

ids = c("SNEMAYT","butterfish","NScod","GBhaddock","ICEherring")
cor <- vector("list",length(ids))
names(cor) <- ids
for(j in 1:length(ids)){
  res_dir <- here("results","v2_fits",paste0(ids[j],"_NAA"))
  mod.list <- list.files(res_dir, full.names = TRUE)
  mods <- lapply(mod.list, readRDS)
  cor[[j]] <- lapply(mods[3:5], get_cor, parnames=c("log_NAA_sigma","trans_NAA_rho"))
}

avg.cor <- vector("list",3)
names(avg.cor) <- c("NAA-2","NAA-3","NAA-4")
for(i in 1:3){
	tmp <- sapply(ids, function(x) cor[[x]][i]) 
	avg.cor[[i]] <- Reduce(`+`, tmp)/length(tmp)
	# avg <- Reduce(`+`, tmp)/length(tmp)
	# avg.cor[[i]] <- avg[upper.tri(avg)]
}
avg.cor
# $`NAA-2`
#                    sigma_R      rho_year
# log_NAA_sigma     1.0000000    -0.1451555
# trans_NAA_rho    -0.1451555     1.0000000

# $`NAA-3`
#                   sigma_R       sigma_a
# log_NAA_sigma    1.00000000    0.02482425
# log_NAA_sigma    0.02482425    1.00000000

# $`NAA-4`
#                  sigma_R       sigma_a   		rho_age        rho_year
# log_NAA_sigma    1.00000000    0.19356092   -0.07486228    -0.1809648
# log_NAA_sigma    0.19356092    1.00000000    0.06867479    -0.4261449
# trans_NAA_rho   -0.07486228    0.06867479    1.00000000    -0.1811612
# trans_NAA_rho   -0.18096481   -0.42614493   -0.18116122     1.0000000

# only correlation > 0.2 is sigma_a and rho_year

ids = c("SNEMAYT","butterfish")
corM <- vector("list",length(ids))
names(corM) <- ids
for(j in 1:length(ids)){
  res_dir <- here("results","v2_fits",paste0(ids[j],"_M"))
  mod.list <- list.files(res_dir, full.names = TRUE)
  mods <- lapply(mod.list, readRDS)
  # cor[[j]] <- lapply(mods[2:3], get_cor, parnames=c("M_repars"))
  corM[[j]] <- get_cor(mods[[3]], parnames=c("M_repars"))
}

# only 1 model, M-3, and 2 stocks
avg.cor <- Reduce(`+`, corM)/length(corM)
avg.cor
#            sigma_M,  rho_M_a,    rho_M_y
# M_repars  1.0000000  0.12091471 -0.56110770
# M_repars  0.1209147  1.00000000 -0.08331715
# M_repars -0.5611077 -0.08331715  1.00000000

# only correlation > 0.2 is sigma_M and rho_M_year

# Selectivity (only Sel-3)
ids = c("GBhaddock")
corS <- vector("list",length(ids))
names(corS) <- ids
for(j in 1:length(ids)){
  res_dir <- here("results","v2_fits",paste0(ids[j],"_sel"))
  mod.list <- list.files(res_dir, full.names = TRUE)
  mods <- lapply(mod.list, readRDS)
  # cor[[j]] <- lapply(mods[2:3], get_cor, parnames=c("M_repars"))
  corS[[j]] <- get_cor(mods[[3]], parnames=c("sel_repars"))
}

corS
# $GBhaddock
#             sigma,    rho,        rho_y
# sel_repars  1.0000000 -0.2664685 -0.2202882
# sel_repars -0.2664685  1.0000000 -0.3877260
# sel_repars -0.2202882 -0.3877260  1.0000000

# highest correlation is rho and rho_year

# average correlation between sigma and rho_year across all stocks and processes
mean(c(sapply(cor, function(x) x[[3]][2,4]), sapply(corM, function(x) x[1,3]), corS[[1]][1,3]))
# -0.4341535

# all NAA correlations
# $SNEMAYT
# $SNEMAYT[[1]]
#               log_NAA_sigma trans_NAA_rho
# log_NAA_sigma       1.00000      -0.15183
# trans_NAA_rho      -0.15183       1.00000

# $SNEMAYT[[2]]
#               log_NAA_sigma log_NAA_sigma
# log_NAA_sigma     1.0000000    -0.1491755
# log_NAA_sigma    -0.1491755     1.0000000

# $SNEMAYT[[3]]
#               log_NAA_sigma log_NAA_sigma trans_NAA_rho trans_NAA_rho
# log_NAA_sigma     1.0000000     0.2978950    -0.2630345    -0.3467301
# log_NAA_sigma     0.2978950     1.0000000    -0.1349108    -0.4336773
# trans_NAA_rho    -0.2630345    -0.1349108     1.0000000    -0.2332206
# trans_NAA_rho    -0.3467301    -0.4336773    -0.2332206     1.0000000


# $butterfish
# $butterfish[[1]]
#               log_NAA_sigma trans_NAA_rho
# log_NAA_sigma     1.0000000    -0.4026259
# trans_NAA_rho    -0.4026259     1.0000000

# $butterfish[[2]]
#               log_NAA_sigma log_NAA_sigma
# log_NAA_sigma     1.0000000     0.1323002
# log_NAA_sigma     0.1323002     1.0000000

# $butterfish[[3]]
#               log_NAA_sigma log_NAA_sigma trans_NAA_rho trans_NAA_rho
# log_NAA_sigma     1.0000000     0.3370353     0.2781503    -0.3314475
# log_NAA_sigma     0.3370353     1.0000000     0.2883951    -0.6671762
# trans_NAA_rho     0.2781503     0.2883951     1.0000000    -0.1455724
# trans_NAA_rho    -0.3314475    -0.6671762    -0.1455724     1.0000000


# $NScod
# $NScod[[1]]
#               log_NAA_sigma trans_NAA_rho
# log_NAA_sigma     1.0000000    -0.1010072
# trans_NAA_rho    -0.1010072     1.0000000

# $NScod[[2]]
#               log_NAA_sigma log_NAA_sigma
# log_NAA_sigma    1.00000000    0.07460788
# log_NAA_sigma    0.07460788    1.00000000

# $NScod[[3]]
#               log_NAA_sigma log_NAA_sigma trans_NAA_rho trans_NAA_rho
# log_NAA_sigma     1.0000000     0.1548010  -0.209988598  -0.215070200
# log_NAA_sigma     0.1548010     1.0000000   0.104742621  -0.301786971
# trans_NAA_rho    -0.2099886     0.1047426   1.000000000  -0.008194331
# trans_NAA_rho    -0.2150702    -0.3017870  -0.008194331   1.000000000


# $GBhaddock
# $GBhaddock[[1]]
#               log_NAA_sigma trans_NAA_rho
# log_NAA_sigma    1.00000000   -0.03944064
# trans_NAA_rho   -0.03944064    1.00000000

# $GBhaddock[[2]]
#               log_NAA_sigma log_NAA_sigma
# log_NAA_sigma   1.000000000  -0.001443196
# log_NAA_sigma  -0.001443196   1.000000000

# $GBhaddock[[3]]
#               log_NAA_sigma log_NAA_sigma trans_NAA_rho trans_NAA_rho
# log_NAA_sigma    1.00000000    0.02139332    0.09340251    0.08005809
# log_NAA_sigma    0.02139332    1.00000000    0.33936393   -0.59602826
# trans_NAA_rho    0.09340251    0.33936393    1.00000000   -0.41119793
# trans_NAA_rho    0.08005809   -0.59602826   -0.41119793    1.00000000


# $ICEherring
# $ICEherring[[1]]
#               log_NAA_sigma trans_NAA_rho
# log_NAA_sigma    1.00000000   -0.03087353
# trans_NAA_rho   -0.03087353    1.00000000

# $ICEherring[[2]]
#               log_NAA_sigma log_NAA_sigma
# log_NAA_sigma    1.00000000    0.06783182
# log_NAA_sigma    0.06783182    1.00000000

# $ICEherring[[3]]
#               log_NAA_sigma log_NAA_sigma trans_NAA_rho trans_NAA_rho
# log_NAA_sigma    1.00000000     0.1566799    -0.2728412   -0.09163432
# log_NAA_sigma    0.15667994     1.0000000    -0.2542170   -0.13205584
# trans_NAA_rho   -0.27284121    -0.2542170     1.0000000   -0.10762083
# trans_NAA_rho   -0.09163432    -0.1320558    -0.1076208    1.00000000
