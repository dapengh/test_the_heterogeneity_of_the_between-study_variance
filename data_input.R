BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
BRD_long <- wide2long(BRD_r)
BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
BRD_pair %>% set_baseline("Florfenicol") %>% set_baseline("Enrofloxacin") %>% set_baseline("No active control") -> BRD_pair
BRD_f <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed =T,comb.random = F)
BRD_pair$index <- ifelse(BRD_pair$treat1 == "No active control", 1, 2)
BRD_pair %>% group_by(studlab) %>% mutate(narms = (sqrt(8*n() + 1) + 1)/2) %>% slice( 1 : ((sqrt(8*n() + 1) - 1)/2) ) %>% ungroup() -> BRD_pair

BRD_pair$studlab <- as.numeric(BRD_pair$studlab)
BRD_pair <- BRD_pair[order(BRD_pair$studlab),]

A <- design_matrix(BRD_pair$treat1, BRD_pair$treat2, baseline = "No active control")
C1 <- C1_func(BRD_pair)

dat <- list(studlab = BRD_pair$studlab, Y = BRD_pair$TE, C1 = C1, A = A, index = BRD_pair$index)

MLE_single <- optim(c(rep(1,ncol(dat$A)), BRD_f$tau^2), wrap_single, dat=dat,
                    method = "L-BFGS-B", lower = c(rep(-Inf, ncol(dat$A)),0), 
                    upper = rep(Inf, ncol(dat$A)+1),
                    control = list(fnscale = -1), hessian = T)

MLE_multiple <- optim(c(rep(2,ncol(dat$A)), c(BRD_f$tau^2, 0.1)), wrap_multiple, dat=dat,
                      method = "L-BFGS-B", lower = c(rep(-Inf, ncol(dat$A)),0, 0), 
                      upper = rep(Inf, ncol(dat$A)+2),
                      control = list(fnscale = -1), hessian = T)

deviance_single <- abs(MLE_single$value)
deviance_multiple <- abs(MLE_multiple$value)
LRT <- - 2 * (MLE_single$value - MLE_multiple$value)
p_multiple <- 1 - pchisq(LRT, 1)
MTCdata <- BRD_pair

n_basic <- length(unique(c(MTCdata$treat1, MTCdata$treat2))) - 1
fisher_info<-solve(-MLE_single$hessian)
prop_sigma<-sqrt(diag(fisher_info))
#C2_single_MLE <- C2_1tau2(dat, MLE_single$par[ncol(A)+1])
#prop_sigma<-sqrt(diag(solve(t(A) %*% solve(C1 + C2_single_MLE) %*% A)))
#prop_sigma[13] <- 0
upper<-MLE_single$par+1.96*prop_sigma
lower<-MLE_single$par-1.96*prop_sigma
MLE_single_res <- data.frame(MLE_single_value=MLE_single$par, MLE_single_se = prop_sigma, MLE_single_lower=lower, MLE_single_upper=upper)
rownames(MLE_single_res) <- c(colnames(A),"tau2")
##MLE_single CI all
MLE_single_res_all <- list()
MLE_single_res_all[[1]] <- t(outer(MLE_single$par[1:n_basic], MLE_single$par[1:n_basic], FUN = "-"))
diag(MLE_single_res_all[[1]]) <- MLE_single$par[1:n_basic]
MLE_single_res_all[[2]] <- outer(prop_sigma[1:n_basic]^2, prop_sigma[1:n_basic]^2, FUN = "+") - 2*fisher_info[1:n_basic, 1:n_basic]
diag(MLE_single_res_all[[2]]) <- prop_sigma[1:n_basic]^2
MLE_single_res_all[[3]] <- MLE_single_res_all[[1]] - qnorm(0.975) * MLE_single_res_all[[2]]
MLE_single_res_all[[4]] <- MLE_single_res_all[[1]] + qnorm(0.975) * MLE_single_res_all[[2]]
names(MLE_single_res_all) <- c("value", "se", "lower", "upper")

fisher_info<-solve(-MLE_multiple$hessian)
prop_sigma<-sqrt(diag(fisher_info))
#C2_multiple_MLE <- C2_2tau2(dat, MLE_multiple$par[(ncol(A)+1) : (ncol(A)+2)])
#prop_sigma<-sqrt(diag(solve(t(A) %*% solve(C1 + C2_multiple_MLE) %*% A)))
#prop_sigma[13:14] <- 0
upper<-MLE_multiple$par+1.96*prop_sigma
lower<-MLE_multiple$par-1.96*prop_sigma
MLE_multiple_res <- data.frame(MLE_multiple_value=MLE_multiple$par, MLE_multiple_se = prop_sigma, MLE_multiple_lower=lower, MLE_multiple_upper=upper)
rownames(MLE_multiple_res) <- c(colnames(A),"tau2_N2A", "tau2_A2A")

MLE_multiple_res_all <- list()
MLE_multiple_res_all[[1]] <- t(outer(MLE_multiple$par[1:n_basic], MLE_multiple$par[1:n_basic], FUN = "-"))
diag(MLE_multiple_res_all[[1]]) <- MLE_multiple$par[1:n_basic]
MLE_multiple_res_all[[2]] <- outer(prop_sigma[1:n_basic]^2, prop_sigma[1:n_basic]^2, FUN = "+") - 2*fisher_info[1:n_basic, 1:n_basic]
diag(MLE_multiple_res_all[[2]]) <- prop_sigma[1:n_basic]^2
MLE_multiple_res_all[[3]] <- MLE_multiple_res_all[[1]] - qnorm(0.975) * MLE_multiple_res_all[[2]]
MLE_multiple_res_all[[4]] <- MLE_multiple_res_all[[1]] + qnorm(0.975) * MLE_multiple_res_all[[2]]
names(MLE_multiple_res_all) <- c("value", "se", "lower", "upper")

###all comparisons
trt_name <- gsub(".*:","",colnames(A))
sig_comp <- NULL
counter <- 1
sig_index <- matrix(0, ncol = 2, nrow = 78)
for(i in 1:ncol(A)){
  for(j in 1:i){
    if(i == j){
      sig_comp[counter] <- paste0("No active control:", trt_name[i])
      sig_index[counter,] <- c(i,i)
    }else{
      sig_comp[counter] <- paste0(trt_name[j],":",trt_name[i])
      sig_index[counter,] <- c(j,i)
    }
    counter <- counter + 1
  }
}
sig_comp_all_df <- data.frame(MLE_single_lower = MLE_single_res_all$lower[sig_index],
                              MLE_single_upper = MLE_single_res_all$upper[sig_index],
                              MLE_multiple_lower = MLE_multiple_res_all$lower[sig_index],
                              MLE_multiple_upper = MLE_multiple_res_all$upper[sig_index],
                              row.names = sig_comp)
row.names(sig_comp_all_df)[11] <- "Enrofloxacin:Florfenicol"
sig_comp_all_df[11, c(1,3,4)] <- -sig_comp_all_df[11, c(1,3,4)]
data_trt <- unique(apply(MTCdata[ , c("treat1", "treat2")], 1, FUN = function(x){paste0(x[1],":",x[2])}))
