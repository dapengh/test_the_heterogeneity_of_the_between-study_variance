
wide2long <- function(MTCdata){
  N <- max(MTCdata$Number.of.arms)
  study <- rep(MTCdata$Study.number,N)
  t <- NULL
  n <- NULL
  r <- NULL
  for(i in c(2,1,3)){
    r <- c(r, eval(parse(text = paste0("MTCdata$Number.of.Event.in.arm.",i, sep = ""))))
    n <- c(n, eval(parse(text = paste0("MTCdata$Total.number.in.arm.",i, sep = ""))))
    t <- c(t, eval(parse(text = paste0("MTCdata$Arm.",i, sep = ""))))
  }
  res <- data.frame(id = study, t = t, r = r, n = n)
  res <- res %>% dplyr::filter(!is.na(n)) %>% arrange(id)
  res
}

set_baseline <- function(data, baseline){
  for(i in 1:dim(data)[1]){
    if(data$treat2[i] == baseline){
      tmp <- data$treat1[i]
      data$treat1[i] <- data$treat2[i]
      data$treat2[i] <- tmp
      data$TE[i] <- -data$TE[i]
    }
  }
  data
}


design_matrix <- function(t1,t2,baseline){
  treat_all <- unique(c(t1,t2))
  treat_others <- treat_all[!(baseline == treat_all)]
  namecol <- paste0(baseline, ":", treat_others, sep = "")
  namecol_2nd <- sub(".*\\:", "", namecol)
  namerow <- paste0(t1, ":", t2, sep = "")
  namerow_1st <- sub("\\:.*", "", namerow)
  namerow_2nd <- sub(".*\\:", "", namerow)
  design_mat <- matrix(0, nrow = length(namerow), ncol = length(namecol))
  for(i in 1:length(namerow)){
    if(namerow[i] %in% namecol){
      design_mat[i,which(namerow[i] == namecol)] <- 1
    }else{
      design_mat[i,which(namerow_1st[i] == namecol_2nd)] <- -1
      design_mat[i,which(namerow_2nd[i] == namecol_2nd)] <- 1
    }
  }
  colnames(design_mat) <- namecol
  design_mat
}

C1_func <- function(dat){
  study <- unique(dat$studlab)
  mat_list <- list()
  for(i in study){
    subgraph <- dat[dat$studlab == i,]
    number = which(study == i)
    if(nrow(subgraph)==1){
      mat_list[[number]] <- subgraph$seTE^2
    }else{
      mat_list[[number]] <- matrix(subgraph$seTE^2, ncol = nrow(subgraph), nrow = nrow(subgraph))
      for(j in 1:nrow(mat_list[[number]])){
        for(k in 1:nrow(mat_list[[number]])){
          if( j != k ){
            mat_list[[number]][j,k] <- 1/subgraph$event1[1] + 1/(subgraph$n1[1] - subgraph$event1[1])
          }
        }
      }
    }
  }
  return(as.matrix(Matrix::bdiag(mat_list)))
}

C2_1tau2 <- function(dat, tau2){
  study <- unique(dat$studlab)
  mat_list <- list()
  for(i in study){
    subgraph <- sum(dat$studlab == i)
    number = which(study == i)
    if(subgraph==1){
      mat_list[[number]] <- tau2
    }else{
      mat_list[[number]] <- matrix(tau2, ncol = subgraph, nrow = subgraph)
      for(j in 1:nrow(mat_list[[number]])){
        for(k in 1:nrow(mat_list[[number]])){
          if( j != k ){
            mat_list[[number]][j,k] <- 1/2 * tau2
          }
        }
      }
    }
  }
  return(as.matrix(Matrix::bdiag(mat_list)))
}

C2_2tau2 <- function(dat, tau2){
  if(tau2[2] > 4 * tau2[1]){tau2[2] = 4*tau2[1]}
  study <- unique(dat$studlab)
  mat_list <- list()
  for(i in study){
    idx <- unique(dat$index[dat$studlab == i])
    subgraph <- sum(dat$studlab == i)
    number = which(study == i)
    if(subgraph==1){
      mat_list[[number]] <- tau2[idx]
    }else{
      mat_list[[number]] <- matrix(tau2[idx], ncol = subgraph, nrow = subgraph)
      for(j in 1:nrow(mat_list[[number]])){
        for(k in 1:nrow(mat_list[[number]])){
          if( j != k ){
            if(idx == 1){mat_list[[number]][j,k] <- tau2[1] - 1/2 *tau2[2]}
            if(idx == 2){mat_list[[number]][j,k] <- 1/2 *tau2[2]}
          }
        }
      }
    }
  }
  return(as.matrix(Matrix::bdiag(mat_list)))
}


wrap_single <- function(parms, dat){
  tau2 <- parms[ncol(dat$A)+1]
  ll_single <- function(mu, tau2, dat){
    C2 <- C2_1tau2(dat, tau2)
    U <- diag(Matrix::chol(dat$C1+C2))
    ll <- -sum(log(U)) - t(dat$Y - dat$A %*% mu) %*% solve(dat$C1 + C2) %*% (dat$Y - dat$A %*% mu) /2 - 1/2 * nrow(dat$A) * log(2*pi)
    return(ll)
  }
  ll_single(mu = parms[1:ncol(dat$A)], tau2 = tau2, dat = dat)
}




wrap_multiple <- function(parms, dat){
  tau2 <- parms[(ncol(dat$A)+1):(ncol(dat$A)+2)]
  ll_multiple <- function(mu, tau2, dat){
    C2 <- C2_2tau2(dat, tau2)
    U <- diag(Matrix::chol(dat$C1+C2))
    ll <- -sum(log(U)) - t(dat$Y - dat$A %*% mu) %*% solve(dat$C1+C2) %*% (dat$Y - dat$A %*% mu) /2 - 1/2 * nrow(dat$A) * log(2*pi)
    return(ll)
  }
  ll_multiple(mu = parms[1:ncol(dat$A)], tau2 = tau2, dat = dat)
}

wrap_single_MLE <- function(parms, dat){
  tau2 <- parms
  ll_single <- function(tau2, dat){
    C1 <- dat$C1 
    C2 <- C2_1tau2(dat, tau2)
    mu <- solve(t(dat$A) %*% solve(dat$C1 + C2) %*% dat$A, t(dat$A) %*% solve(dat$C1 + C2) %*% dat$Y)
    U <- diag(Matrix::chol(dat$C1+C2))
    ll <- -sum(log(U)) - t(dat$Y - dat$A %*% mu) %*% solve(dat$C1 + C2) %*% (dat$Y - dat$A %*% mu) /2 - 1/2 * nrow(dat$A) * log(2*pi)
    return(ll)
  }
  ll_single(tau2 = tau2, dat = dat)
}

wrap_multiple_MLE <- function(parms, dat){
  tau2 <- parms[1:2]
  ll_multiple <- function(tau2, dat){
    C2 <- C2_2tau2(dat, tau2)
    mu <- solve(t(dat$A) %*% solve(dat$C1 + C2) %*% dat$A, t(dat$A) %*% solve(dat$C1 + C2) %*% dat$Y)
    U <- diag(Matrix::chol(dat$C1+C2))
    ll <- -sum(log(U)) - t(dat$Y - dat$A %*% mu) %*% solve(dat$C1 + C2) %*% (dat$Y - dat$A %*% mu) /2 - 1/2 * nrow(dat$A) * log(2*pi)
    return(ll)
  }
  ll_multiple(tau2 = tau2,  dat = dat)
}



smry_func <- function(tmp){
  LRT <- -2*(tmp$ll_single - tmp$ll_multiple)
  p_value <- 1 - pchisq(LRT, 1)
  if(tmp$n_tau2[1] == 1){
    return(list(Type_1_error = mean(p_value <= 0.05), zero_prop_single = mean(tmp$tau2_single == 0),
                neg_prop_single = mean(tmp$tau2_single < 0),
                zero_prop_multiple = (mean(tmp$tau2_multiple1 == 0) + mean(tmp$tau2_multiple2 == 0))/2,
                neg_prop_multiple = (mean(tmp$tau2_multiple1 < 0) + mean(tmp$tau2_multiple2 < 0))/2))
  }
  if(tmp$n_tau2[1] == 2){
    return(list(Power = mean(p_value <= 0.05), zero_prop_single = mean(tmp$tau2_single == 0),
                neg_prop_single = mean(tmp$tau2_single < 0),
                zero_prop_multiple = (mean(tmp$tau2_multiple1 == 0) + mean(tmp$tau2_multiple2 == 0))/2,
                neg_prop_multiple = (mean(tmp$tau2_multiple1 < 0) + mean(tmp$tau2_multiple2 < 0))/2))
  }
}


simu_MLE_data <- function(dat, d_true, tau2, n_boot){
  dat_tmp <- dat
  d <- dat$A %*% d_true
  LRT <- numeric(n_boot)
  if(length(tau2) == 1){
    C2 <- C2_1tau2(dat, tau2)
    for(i in 1:n_boot){
      dat_tmp$Y <- MASS::mvrnorm(1, d, dat$C1 + C2)
      tmp_single_MLE <- optim(c(tau2), wrap_single_MLE, dat=dat_tmp,
                              method = "L-BFGS-B", lower = c( 0), 
                              upper = rep(Inf, 1),
                              control = list(fnscale = -1, maxit = 200))
      tmp_multiple_MLE <- optim(c(tau2, tau2), wrap_multiple_MLE, dat=dat_tmp,
                                method = "L-BFGS-B", lower = c( 0), 
                                upper = rep(Inf, 1),
                                control = list(fnscale = -1, maxit = 200))
      LRT[i] <- -2*(tmp_single_MLE$value - tmp_multiple_MLE$value)
    }
  }
  return(LRT)
}

simu_MLE_boot <- function(dat, d_true, tau2, n, n_boot){
  dat_tmp <- dat
  d <- dat$A %*% d_true
  d_all <- t(outer(d_true, d_true, FUN = "-"))
  diag(d_all) <- d_true
  single_cover_cumm <- numeric(choose(ncol(dat$A)+1,2))
  multiple_cover_cumm <- numeric(choose(ncol(dat$A)+1,2))
  n_basic <- ncol(dat$A)
  if(length(tau2) == 2){
    C2 <- C2_2tau2(dat, tau2)
    for(i in 1:n){
      dat_tmp$Y <- MASS::mvrnorm(1, d, dat$C1 + C2)
      tmp_single_MLE <- optim(c(tau2[1]), wrap_single_MLE, dat=dat_tmp,
                              method = "L-BFGS-B", lower = c( 0), 
                              upper = rep(Inf, 1),
                              control = list(fnscale = -1, maxit = 200))
      C2_single_MLE <- C2_1tau2(dat_tmp, tmp_single_MLE$par)
      mu <- solve(t(dat_tmp$A) %*% solve(dat_tmp$C1 + C2_single_MLE) %*% dat_tmp$A, t(dat_tmp$A) %*% solve(dat_tmp$C1 + C2_single_MLE) %*% dat_tmp$Y)
      fisher_info <- solve(t(dat_tmp$A) %*% solve(dat_tmp$C1 + C2_single_MLE) %*% dat_tmp$A)
      prop_sigma<-sqrt(diag(fisher_info))
      upper<-mu+qnorm(0.975)*prop_sigma
      lower<-mu-qnorm(0.975)*prop_sigma
      MLE_single_res <- data.frame(MLE_single_value=mu, MLE_single_se = prop_sigma, MLE_single_lower=lower, MLE_single_upper=upper)
      rownames(MLE_single_res) <- c(colnames(dat_tmp$A))
      ##MLE_single CI all
      MLE_single_res_all <- list()
      MLE_single_res_all[[1]] <- t(outer(MLE_single$par[1:n_basic], MLE_single$par[1:n_basic], FUN = "-"))
      diag(MLE_single_res_all[[1]]) <- MLE_single$par[1:n_basic]
      MLE_single_res_all[[2]] <- outer(prop_sigma[1:n_basic]^2, prop_sigma[1:n_basic]^2, FUN = "+") - 2*fisher_info[1:n_basic, 1:n_basic]
      diag(MLE_single_res_all[[2]]) <- prop_sigma[1:n_basic]^2
      MLE_single_res_all[[3]] <- MLE_single_res_all[[1]] - qnorm(0.975) * MLE_single_res_all[[2]]
      MLE_single_res_all[[4]] <- MLE_single_res_all[[1]] + qnorm(0.975) * MLE_single_res_all[[2]]
      names(MLE_single_res_all) <- c("value", "se", "lower", "upper")
      
      
      tmp_multiple_MLE <- optim(c(tau2), wrap_multiple_MLE, dat=dat_tmp,
                                method = "L-BFGS-B", lower = c( 0,0), 
                                upper = rep(Inf, 2),
                                control = list(fnscale = -1, maxit = 200))
      C2_multiple_MLE <- C2_2tau2(dat = dat_tmp, tmp_multiple_MLE$par)
      mu <- solve(t(dat_tmp$A) %*% solve(dat_tmp$C1 + C2_multiple_MLE) %*% dat_tmp$A, t(dat_tmp$A) %*% solve(dat_tmp$C1 + C2_multiple_MLE) %*% dat_tmp$Y)
      fisher_info <- solve(t(dat_tmp$A) %*% solve(dat_tmp$C1 + C2_multiple_MLE) %*% dat_tmp$A)
      prop_sigma<-sqrt(diag(fisher_info))
      upper<-mu+qnorm(0.975)*prop_sigma
      lower<-mu-qnorm(0.975)*prop_sigma
      MLE_multiple_res <- data.frame(MLE_multiple_value=mu, MLE_multiple_se = prop_sigma, MLE_multiple_lower=lower, MLE_multiple_upper=upper)
      rownames(MLE_multiple_res) <- c(colnames(dat_tmp$A))
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
            sig_comp[counter] <- paste0("NAC:", trt_name[i])
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
                                    eff_true = d_all[sig_index],
                                    row.names = sig_comp)
      single_cover <- (sig_comp_all_df$eff_true >= sig_comp_all_df$MLE_single_lower) &  
        (sig_comp_all_df$eff_true <= sig_comp_all_df$MLE_single_upper)
      multiple_cover <- (sig_comp_all_df$eff_true >= sig_comp_all_df$MLE_multiple_lower) &  
        (sig_comp_all_df$eff_true <= sig_comp_all_df$MLE_multiple_upper)
      single_cover_cumm <- single_cover_cumm + single_cover
      multiple_cover_cumm <- multiple_cover_cumm + multiple_cover
    }
  }
  return(list(single_cover_rate = single_cover_cumm/n, multiple_cover_rate = multiple_cover_cumm/n))
}

