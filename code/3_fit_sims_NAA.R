# Brian Stock
# June 15 2020
# Simulation test WHAM

# Assumes you open R in project directory
# library(here); source(here("code","3_fit_sims_NAA.R"))

# Step 1: Fit 4 NAA models to 2019 SNE-MA yellowtail flounder data

# Step 2: Simulate 4 NAA operating models:
#  1. rec 	iid
#  2. rec 	ar1_y
#  3. rec+1 iid
#  4. rec+1 2dar1
# For each OM, simulate using
#  a. obs error only (fix NAA at estimated values)
#  b. obs + process error (new NAA)

# Step 3: Fit OMs to simulated datasets
# Step 4: Collect + plot results
# -----------------------------------------------------------------------
# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "om_mode")
library(wham)
library(here)
library(tidyverse)

# functions to calculate F/F40 and B/B40
calc_relF <- function(mod, type="fit"){
	if(type == "fit"){
		F.t <- mod$rep$log_FXSPR
		faa <- log(mod$rep$FAA_tot)		
	}
	if(type == "sim"){
		F.t <- mod$env$data$log_FXSPR
		faa <- log(mod$env$data$FAA_tot)			
	}
  	age.full.f <- apply(faa,1, function(x) max(which(x == max(x))))
  	full.f <- faa[cbind(seq_along(age.full.f),age.full.f)]
	rel.f <- exp(full.f - F.t)
	return(rel.f)
}
calc_relB <- function(mod, type="fit"){
	if(type == "fit"){
		SSB.t <- exp(mod$rep$log_SSB_FXSPR)
		ssb <- mod$rep$SSB
	}
	if(type == "sim"){
		SSB.t <- exp(mod$env$data$log_SSB_FXSPR)
		ssb <- mod$env$data$SSB
	}
	rel.ssb <- ssb / SSB.t	
	return(rel.ssb)
}

# Assumes you open R in project directory
# Load fit models from step 1
# mod.list <- here("results","NAA",paste0("m",1:4,".rds"))
# mods <- lapply(mod.list, readRDS)
# input.list <- here("results","NAA",paste0("m",1:4,"_input.rds"))
# inputs <- lapply(input.list, readRDS)

n.mods <- 4
n.sim <- 100
n.types <- 2
mod <- readRDS(here("results","NAA","m1.rds"))
n.years <- length(mod$years)
rm("mod")
sim.seeds <- readRDS(here("data","sim_seeds.rds"))

options(warn=-1) # suppress warning messages

nested.list <- function(len) if(length(len) == 1) vector("list", len) else lapply(1:len[1], function(...) nested.list(len[-1]))
sdreps <- nested.list(c(n.types, n.mods, n.sim))
res.colnames <- c("om","em","type","year","sim","F_fit","F_sim","relF_fit","relF_sim","SSB_fit","SSB_sim","relB_fit","relB_sim","catch_fit","catch_sim",paste0("NAA",1:6))
results <- rep(list(rep(list(rep(list(matrix(NA, ncol = length(res.colnames), nrow = n.years)),n.sim)),n.mods)),n.types) # nested lists with preallocated matrices
# colnames(results) <- res.colnames
for(m in 1:n.mods){
	simdata <- readRDS(here("data","simdata","SNEMAYT","NAA",paste0("simdata_om",m,".rds")))
	for(i in 1:n.sim){
		print(paste0("Model: ",m," Sim: ", i))
		set.seed(sim.seeds[i])

		# a) obs error, keep all parameters (incl NAA) as in fit model, simulate catch + index data
		input1 <- readRDS(here("results","NAA",paste0("m",m,"_input.rds")))
		input1$data <- simdata[[i]][[1]]
		fit1 <- fit_wham(input1, do.sdrep=T, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)
		if(exists("err")) rm("err") # need to clean this up
		df <- as.matrix(data.frame(om=m, em=m, type=1, year=fit1$years, sim=i, 
			F_fit=fit1$rep$F[,1], F_sim=fit1$env$data$F[,1], relF_fit=calc_relF(fit1, type="fit"), relF_sim=calc_relF(fit1, type="sim"), 
			SSB_fit=fit1$rep$SSB, SSB_sim=fit1$env$data$SSB, relB_fit=calc_relB(fit1, type="fit"), relB_sim=calc_relB(fit1, type="sim"),  
			catch_fit=fit1$rep$pred_catch[,1], catch_sim=fit1$env$data$pred_catch[,1]))
		dfnaa <- fit1$rep$NAA; colnames(dfnaa) <- paste0("NAA",1:6)
		# results[seq(((i-1)*n.types*n.years)) + 1:n.years,] <- cbind(df, dfnaa)
		# results <- rbind(results, cbind(df, dfnaa))
		results[[1]][[m]][[i]] <- cbind(df, dfnaa)
		# results[[1]] <- cbind(df, dfnaa)
		sdreps[[1]][[m]][[i]] <- summary(fit1$sdrep)

		# b) process + obs error, keep parameters except NAA as in fit model, simulate NAA (process error) and catch + index data (obs error)
		input2 <- readRDS(here("results","NAA",paste0("m",m,"_input.rds")))
		input2$data <- simdata[[i]][[2]]
		fit2 <- fit_wham(input2, do.sdrep=T, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)		
		if(exists("err")) rm("err") # need to clean this up
		df <- as.matrix(data.frame(om=m, em=m, type=2, year=fit2$years, sim=i, 
			F_fit=fit2$rep$F[,1], F_sim=fit2$env$data$F[,1], relF_fit=calc_relF(fit2, type="fit"), relF_sim=calc_relF(fit2, type="sim"), 
			SSB_fit=fit2$rep$SSB, SSB_sim=fit2$env$data$SSB, relB_fit=calc_relB(fit2, type="fit"), relB_sim=calc_relB(fit2, type="sim"),  
			catch_fit=fit2$rep$pred_catch[,1], catch_sim=fit2$env$data$pred_catch[,1]))
		dfnaa <- fit2$rep$NAA; colnames(dfnaa) <- paste0("NAA",1:6)
		# results <- rbind(results, cbind(df, dfnaa))
		results[[2]][[m]][[i]] <- cbind(df, dfnaa)
		# results[[2]] <- cbind(df, dfnaa)
		sdreps[[2]][[m]][[i]] <- summary(fit2$sdrep)

		# save simulated data for cross-tests
		# saveRDS(results, file=here("results","NAA",paste0("res_om",m,"_sim_",i,".rds")))
		# saveRDS(simdata, file=here("data","simdata","SNEMAYT","NAA",paste0("om",m,"_sim_",i,".rds")))
		# rm(list=c("simdata","fit_ai","fit_bi","df","dfnaa"))
		rm(list=c("input1","input2","fit1","fit2","df","dfnaa"))
	}
	rm(list=c("simdata"))
}

# Collect results into one matrix
flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
res2 <- do.call(rbind, flatten.nested.list(results)) %>% as.data.frame
saveRDS(res2, file=here("results","NAA","results.rds"))
saveRDS(sdreps, file=here("results","NAA","sdreps.rds"))

# un-suppress warnings
options(warn=0)

# 1. files to delete (test first)
# find | grep res_om | xargs ls -lh
# 2. delete the files
# find | grep res_om | xargs rm -f


# # improved list of objects
# .ls.objects <- function (pos = 1, pattern, order.by,
#                         decreasing=FALSE, head=FALSE, n=5) {
#     napply <- function(names, fn) sapply(names, function(x)
#                                          fn(get(x, pos = pos)))
#     names <- ls(pos = pos, pattern = pattern)
#     obj.class <- napply(names, function(x) as.character(class(x))[1])
#     obj.mode <- napply(names, mode)
#     obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
#     obj.prettysize <- napply(names, function(x) {
#                            format(utils::object.size(x), units = "auto") })
#     obj.size <- napply(names, object.size)
#     obj.dim <- t(napply(names, function(x)
#                         as.numeric(dim(x))[1:2]))
#     vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
#     obj.dim[vec, 1] <- napply(names, length)[vec]
#     out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
#     names(out) <- c("Type", "Size", "PrettySize", "Length/Rows", "Columns")
#     if (!missing(order.by))
#         out <- out[order(out[[order.by]], decreasing=decreasing), ]
#     if (head)
#         out <- head(out, n)
#     out
# }

# # shorthand
# lsos <- function(..., n=10) {
#     .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
# }

# lsos()
