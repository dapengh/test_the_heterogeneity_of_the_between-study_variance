packages <- c("dplyr", "netmeta", "MASS", "Matrix", "parallel", "foreach", "doParallel", "doRNG")
lapply(packages, library, character.only = TRUE)
source("./BRD_func_3arm.R")
source("./data_input.R")

### model constrain: 0
### study size for each comparison: 1*(66, 40)
### True one tau2: MLE_single$par[ncol(A) + 1]
### True barcesic parameter: MLE_single estimates
### Within-study variance dist: same as real data
set.seed(20190915)
dat <- real_setting(MTCdata)
cl <- makeCluster(16)
registerDoParallel(cl)
set.seed(201991234)
tmp <- foreach(i = 1:100) %dorng% {
  simu_MLE_boot(dat = dat, d_true = MLE_multiple$par[1:ncol(dat$A)], tau2 = MLE_multiple$par[(ncol(dat$A) + 1):(ncol(dat$A) + 2)], n = 1, n_boot = 1000)
}

stopCluster(cl)
res_real_ns1_1tau2_1_zero <- mapply(unlist, tmp) %>% t() %>% as.data.frame()
write.csv(res_real_ns1_1tau2_1_zero, file = "./res/res_real_4.csv", row.names = F)