
rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, gamma_sd = 1, 
                         gamma_mn = 1, C = 1){
  tot_res <- NULL;
  fea_res <- NULL;
  vars_t  <- NULL;
  for(i in 2:max_sp){
    iter           <- iters;
    tot_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 7);
    fea_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 7);
    while(iter > 0){
      r_vec    <- rnorm(n = i, mean = 0, sd = rmx);
      A0_dat   <- rnorm(n = i * i, mean = 0, sd = 0.4);
      A0       <- matrix(data = A0_dat, nrow = i, ncol = i);
      A0       <- species_interactions(mat = A0, type = int_type);
      C_dat    <- rbinom(n = i * i, size = 1, prob = C);
      C_mat    <- matrix(data = C_dat, nrow = i, ncol = i);
      A0       <- A0 * C_mat;
      diag(A0) <- -1;
      gam0     <- make_gammas(nn = i, distribution = 0, sdd = gamma_sd);
      gam1     <- make_gammas(nn = i, distribution = 1, sdd = gamma_sd);
      gam2     <- make_gammas(nn = i, distribution = 2, sdd = gamma_sd);
      gam3     <- make_gammas(nn = i, distribution = 3, sdd = gamma_sd);
      gam4     <- make_gammas(nn = i, distribution = 4, sdd = gamma_sd);
      A1       <- A0 * gam1;
      A2       <- A0 * gam2;
      A3       <- A0 * gam3;
      A4       <- A0 * gam4;
      A0       <- A0 * gam0;
      A0_stb   <- max(Re(eigen(A0)$values)) < 0;
      A1_stb   <- max(Re(eigen(A1)$values)) < 0;
      A2_stb   <- max(Re(eigen(A2)$values)) < 0;
      A3_stb   <- max(Re(eigen(A3)$values)) < 0;
      A4_stb   <- max(Re(eigen(A4)$values)) < 0;
      A0_fea   <- min(-1*solve(A0) %*% r_vec) > 0;
      A1_fea   <- min(-1*solve(A1) %*% r_vec) > 0;
      A2_fea   <- min(-1*solve(A2) %*% r_vec) > 0;
      A3_fea   <- min(-1*solve(A3) %*% r_vec) > 0;
      A4_fea   <- min(-1*solve(A4) %*% r_vec) > 0;
      if(A0_stb == TRUE){
        tot_res[[i-1]][iter, 1] <- 1;
      }
      if(A1_stb == TRUE){
        tot_res[[i-1]][iter, 2] <- 1;
      }
      if(A2_stb == TRUE){
        tot_res[[i-1]][iter, 3] <- 1;
      }
      if(A3_stb == TRUE){
        tot_res[[i-1]][iter, 4] <- 1;
      }
      if(A4_stb == TRUE){
        tot_res[[i-1]][iter, 5] <- 1;
      }
      if(A0_fea == TRUE){
        fea_res[[i-1]][iter, 1] <- 1;
      }
      if(A1_fea == TRUE){
        fea_res[[i-1]][iter, 2] <- 1;
      }
      if(A2_fea == TRUE){
        fea_res[[i-1]][iter, 3] <- 1;
      }
      if(A3_fea == TRUE){
        fea_res[[i-1]][iter, 4] <- 1;
      }
      if(A4_fea == TRUE){
        fea_res[[i-1]][iter, 5] <- 1;
      }
      new_v   <- get_var_A(A1);
      vars_t  <- rbind(vars_t, c(i, new_v));
      iter    <- iter - 1;
    }
    print(i);
  }
  all_res <- summarise_randmat(tot_res = tot_res, fea_res = fea_res);
  return(list(all_res = all_res, vars_t = vars_t));
}

get_var_A <- function(A){
  off_vec  <- NULL;
  diag_vec <- NULL;
  for(i in 1:dim(A)[1]){
    for(j in 1:dim(A)[2]){
      if(i == j){
        diag_vec <- c(diag_vec, A[i,j]);
      }else{
        off_vec  <- c(off_vec, A[i, j])
      }
    }
  }
  v_diag <- var(diag_vec);
  v_off  <- var(off_vec);
  vals   <- c(v_diag, v_off);
  return(vals);
}


summ <- NULL;
for(i in 2:30){
  dat <- sim$vars_t[sim$vars_t[,1] == i,];
  row <- apply(X = dat, MARGIN = 2, FUN = mean);
  summ <- rbind(summ, row);
}

