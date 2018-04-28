
rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, gamma_sd = 1, 
                         gamma_mn = 1, C = 1){
  tot_res  <- NULL;
  fea_res  <- NULL;
  vars_t   <- NULL;
  off_cors <- NULL;
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
      
      gam1     <- runif(n = i, min = 0, max = 999) #make_gammas(nn = i, distribution = 1, sdd = gamma_sd);
      gam0     <- rep(mean(gam1), i);
      
      gam2     <- make_gammas(nn = i, distribution = 2, sdd = gamma_sd);
      gam3     <- make_gammas(nn = i, distribution = 3, sdd = gamma_sd);
      gam4     <- make_gammas(nn = i, distribution = 4, sdd = gamma_sd);
      
      A1       <- A0 * gam1;
      A2       <- A0 * gam1;
      A3       <- A0 * gam0;
      diag(A2) <- -1 * gam0;
      diag(A3) <- -1 * gam1;
      

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
      
      stbs     <- c(A0_stb, A1_stb, A2_stb, A3_stb, A4_stb);
      
      v0       <- get_var_A(A0);
      v1       <- get_var_A(A1);
      v2       <- get_var_A(A2);
      v3       <- get_var_A(A3);
      v4       <- get_var_A(A4);
      vs       <- c(v0, v1, v2, v3, v4);
      vars_t   <- rbind(vars_t, c(i, vs, stbs));
      
      o0       <- get_off_corr(A0);
      o1       <- get_off_corr(A1);
      o2       <- get_off_corr(A2);
      o3       <- get_off_corr(A3);
      o4       <- get_off_corr(A4);
      os       <- c(o0, o1, o2, o3, o4);
      off_cors <- rbind(off_cors, c(i, os, stbs));      
      
      iter    <- iter - 1;
    }
    print(i);
  }
  all_res <- summarise_randmat(tot_res = tot_res, fea_res = fea_res);
  return(list(all_res = all_res, vars_t = vars_t, off_cors = off_cors));
}

get_var_A <- function(A){
  off_vec  <- NULL;
  for(i in 1:dim(A)[1]){
    for(j in 1:dim(A)[2]){
      if(i != j){
        off_vec  <- c(off_vec, A[i, j])
      }
    }
  }
  v_off  <- var(off_vec);
  return(v_off);
}

get_off_corr <- function(A){
  pairs <- NULL;
  for(i in 1:dim(A)[1]){
    for(j in 1:dim(A)[2]){
      if(i > j){
        pairs  <- rbind(pairs, c(A[j, i], A[i, j]));
      }
    }
  }
  rho <- cor(pairs[,1], pairs[,2]);
  return(rho);
}

summ <- NULL;
for(i in 2:30){
  dat <- sim$vars_t[sim$vars_t[,1] == i,];
  row <- apply(X = dat, MARGIN = 2, FUN = mean);
  summ <- rbind(summ, row);
}















################################################################################




get_var_A <- function(A){
    off_vec  <- NULL;
    for(i in 1:dim(A)[1]){
        for(j in 1:dim(A)[2]){
            if(i != j){
                off_vec  <- c(off_vec, A[i, j])
            }
        }
    }
    v_off  <- var(off_vec);
    return(v_off);
}

get_mn_A <- function(A){
    off_vec  <- NULL;
    for(i in 1:dim(A)[1]){
        for(j in 1:dim(A)[2]){
            if(i != j){
                off_vec  <- c(off_vec, A[i, j])
            }
        }
    }
    v_off  <- mean(off_vec);
    return(v_off);
}

get_mn_abs_A <- function(A){
    off_vec  <- NULL;
    for(i in 1:dim(A)[1]){
        for(j in 1:dim(A)[2]){
            if(i != j){
                off_vec  <- c(off_vec, A[i, j])
            }
        }
    }
    v_off  <- mean(abs(off_vec));
    return(v_off);
}



summ_stats <- NULL;
summ_eigs0 <- NULL;
summ_eigs1 <- NULL;
for(i in 1:10000){
    A_comp <- NULL;
    A_dat  <- rnorm(n = 256, mean = 0, sd = 0.4);
    A_mat  <- matrix(data = A_dat, nrow = 16);
    gammas <- runif(n = 16, min = 0, max = 1);
    mu_gam <- mean(gammas);
    diag(A_mat) <- -1;
    A1     <- gammas * A_mat;
    A0     <- mu_gam * A_mat;
    vA0    <- get_var_A(A0);
    vA1    <- get_var_A(A1);
    mA0    <- get_mn_A(A0);
    mA1    <- get_mn_A(A1);
    aA0    <- get_mn_abs_A(A0);
    aA1    <- get_mn_abs_A(A1);
    A0_e   <- eigen(A0)$values;
    A0_r   <- Re(A0_e);
    A0_i   <- Im(A0_e);
    A1_e   <- eigen(A1)$values;
    A1_r   <- Re(A1_e);
    A1_i   <- Im(A1_e);

    A0_s   <- max(A0_r) < 0;
    A1_s   <- max(A1_r) < 0;
    eigV0  <- cbind(A0_r, A0_i);
    eigV1  <- cbind(A1_r, A1_i); 
    sum_v  <- c(vA0, vA1, mA0, mA1, aA0, aA1, A0_s, A1_s);
    summ_eigs0 <- rbind(summ_eigs0, eigV0);
    summ_eigs1 <- rbind(summ_eigs1, eigV1);
    summ_stats <- rbind(summ_stats, sum_v);
    if(i %% 1000 == 0){
        print(i);    
    }
}





t0 <- A0;
diag(t0) <- NA;
t1 <- A1;
diag(t1) <- NA;



summ_stats <- NULL;
summ_eigs0 <- NULL;
summ_eigs1 <- NULL;
for(i in 1:10000){
    A_comp <- NULL;
    A_dat  <- rnorm(n = 256, mean = 0, sd = 0.4);
    A_mat  <- matrix(data = A_dat, nrow = 16);
    gammas <- runif(n = 16, min = 0, max = 2);
    mu_gam <- mean(gammas);
    diag(A_mat) <- -1;
    s_gams <- sort(gammas, decreasing = FALSE);
    rsums  <- apply(abs(A_mat), 1, sum);
    osums  <- (length(rsums) + 1) - rank(rsums);
    ogams  <- s_gams[osums];
    A1     <- ogams  * A_mat;
    A0     <- mu_gam * A_mat;
    vA0    <- get_var_A(A0);
    vA1    <- get_var_A(A1);
    mA0    <- get_mn_A(A0);
    mA1    <- get_mn_A(A1);
    aA0    <- get_mn_abs_A(A0);
    aA1    <- get_mn_abs_A(A1);
    A0_e   <- eigen(A0)$values;
    A0_r   <- Re(A0_e);
    A0_i   <- Im(A0_e);
    A1_e   <- eigen(A1)$values;
    A1_r   <- Re(A1_e);
    A1_i   <- Im(A1_e);
    A0_s   <- max(A0_r) < 0;
    A1_s   <- max(A1_r) < 0;
    eigV0  <- cbind(A0_r, A0_i);
    eigV1  <- cbind(A1_r, A1_i);
    sum_v  <- c(vA0, vA1, mA0, mA1, aA0, aA1, A0_s, A1_s);
    
    summ_eigs0 <- rbind(summ_eigs0, eigV0);
    summ_eigs1 <- rbind(summ_eigs1, eigV1);
    summ_stats <- rbind(summ_stats, sum_v);
    if(i %% 1000 == 0){
        print(i);    
    }
}
















summ_stats <- NULL;
summ_eigs0 <- NULL;
summ_eigs1 <- NULL;
for(i in 1:10000){
    A_comp <- NULL;
    A_dat  <- rnorm(n = 256, mean = 0, sd = 0.4);
    A_mat  <- matrix(data = A_dat, nrow = 16);
    gammas <- runif(n = 16, min = 0, max = 2);
    mu_gam <- mean(gammas);
    A1     <- A_mat;
    A0     <- A_mat;
    diag(A1) <- runif(n = 16, min = -2, max = 0);
    diag(A0) <- mean(diag(A1));
    vA0    <- get_var_A(A0);
    vA1    <- get_var_A(A1);
    mA0    <- get_mn_A(A0);
    mA1    <- get_mn_A(A1);
    aA0    <- get_mn_abs_A(A0);
    aA1    <- get_mn_abs_A(A1);
    A0_e   <- eigen(A0)$values;
    A0_r   <- Re(A0_e);
    A0_i   <- Im(A0_e);
    A1_e   <- eigen(A1)$values;
    A1_r   <- Re(A1_e);
    A1_i   <- Im(A1_e);
    A0_s   <- max(A0_r) < 0;
    A1_s   <- max(A1_r) < 0;
    eigV0  <- cbind(A0_r, A0_i);
    eigV1  <- cbind(A1_r, A1_i);
    sum_v  <- c(vA0, vA1, mA0, mA1, aA0, aA1, A0_s, A1_s);
    
    summ_eigs0 <- rbind(summ_eigs0, eigV0);
    summ_eigs1 <- rbind(summ_eigs1, eigV1);
    summ_stats <- rbind(summ_stats, sum_v);
    if(i %% 1000 == 0){
        print(i);    
    }
}



























A_comp <- NULL;
A_dat  <- rnorm(n = 1000000, mean = 0, sd = 0.4);
A_mat  <- matrix(data = A_dat, nrow = 1000);
gammas <- runif(n = 1000, min = 0, max = 1);
mu_gam <- mean(gammas);
diag(A_mat) <- -1;
A1     <- gammas * A_mat;
A0     <- mu_gam * A_mat;
A0_e   <- eigen(A0)$values;
A0_r   <- Re(A0_e);
A0_i   <- Im(A0_e);
A1_e   <- eigen(A1)$values;
A1_r   <- Re(A1_e);
A1_i   <- Im(A1_e);

A0_vm       <- A0;
diag(A0_vm) <- NA;
A0vec       <- as.vector(A0_vm);
A0vec       <- A0vec[is.na(A0vec) == FALSE];
A1_vm       <- A1;
diag(A1_vm) <- NA;
A1vec       <- as.vector(A1_vm);
A1vec       <- A1vec[is.na(A1vec) == FALSE];

par(mar = c(5, 5, 1, 1));
plot(A0_r, A0_i, xlim = c(-8.75, 8.25), ylim = c(-8.5,8.5), pch = 4, cex = 0.7,
     xlab = "Real", ylab = "Imaginary", cex.lab = 1.5, cex.axis = 1.5, asp = 1);
points(A1_r, A1_i, col = "red", pch = 4, cex = 0.7)

vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(1000) * sd(A0vec) * cos(vl) + mean(diag(A0));
y0 <- sqrt(1000) * sd(A0vec) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3)


vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(1000) * sd(A1vec) * cos(vl) + mean(diag(A1));
y0 <- sqrt(1000) * sd(A1vec) * sin(vl);
points(x = x0, y = y0, type = "l", col = "red", lwd = 3)
















A_comp <- NULL;
A_dat  <- rnorm(n = 400, mean = 0, sd = 0.4);
A_mat  <- matrix(data = A_dat, nrow = 20);
gammas <- runif(n = 20, min = 0, max = 1);
mu_gam <- mean(gammas);
diag(A_mat) <- -1;
A1     <- gammas * A_mat;
A0     <- mu_gam * A_mat;
A0_e   <- eigen(A0)$values;
A0_r   <- Re(A0_e);
A0_i   <- Im(A0_e);
A1_e   <- eigen(A1)$values;
A1_r   <- Re(A1_e);
A1_i   <- Im(A1_e);

par(mar = c(5, 5, 1, 1));
plot(A0_r, A0_i, xlim = c(-2.75, 2.25), ylim = c(-2.5,2.5), pch = 4, cex = 0.7,
     xlab = "Real", ylab = "Imaginary", cex.lab = 1.5, cex.axis = 1.5)
points(A1_r, A1_i, col = "red", pch = 4, cex = 0.7)

vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(20) * sd(A0vec) * cos(vl) + mean(diag(A0));
y0 <- sqrt(20) * sd(A0vec) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3)


vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(20) * sd(A1vec) * cos(vl) + mean(diag(A1));
y0 <- sqrt(20) * sd(A1vec) * sin(vl);
points(x = x0, y = y0, type = "l", col = "red", lwd = 3)












































rvls0  <- NULL;
ivls0  <- NULL;
rvls1  <- NULL;
ivls1  <- NULL;
iters  <- 1000;

while(iters > 0){
    A_comp <- NULL;
    A_dat  <- rnorm(n = 400, mean = 0, sd = 0.4);
    A_mat  <- matrix(data = A_dat, nrow = 20);
    gammas <- runif(n = 20, min = 0, max = 100);
    mu_gam <- mean(gammas);
    diag(A_mat) <- -1;
    A1     <- gammas * A_mat;
    A0     <- mu_gam * A_mat;
    A0_e   <- eigen(A0)$values;
    A0_r   <- Re(A0_e);
    A0_i   <- Im(A0_e);
    A1_e   <- eigen(A1)$values;
    A1_r   <- Re(A1_e);
    A1_i   <- Im(A1_e);
    rvls0  <- c(rvls0, A0_r);
    ivls0  <- c(ivls0, A0_i);
    rvls1  <- c(rvls1, A1_r);
    ivls1  <- c(ivls1, A1_i);
    iters  <- iters - 1;
}  

sum(rvls0 < 0);
sum(rvls1 < 0);

sum(rvls0 < 0) / sum(rvls1 < 0);

par(mar = c(5, 5, 1, 1));
plot(rvls0, ivls0, xlim = c(-2, 1), ylim = c(-1.5,1.5), pch = 4, cex = 0.7,
     xlab = "Real", ylab = "Imaginary", cex.lab = 1.5, cex.axis = 1.5)
points(rvls1, ivls1, col = "red", pch = 4, cex = 0.7)

vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(20) * sd(A0vec) * cos(vl) + mean(diag(A0));
y0 <- sqrt(20) * sd(A0vec) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3)


vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(20) * sd(A1vec) * cos(vl) + mean(diag(A1));
y0 <- sqrt(20) * sd(A1vec) * sin(vl);
points(x = x0, y = y0, type = "l", col = "red", lwd = 3)











###########################################################################

# Pr stability is higher given variation in gamma (S = 16 as an example)

summ_stats <- NULL;
summ_eigs0 <- NULL;
summ_eigs1 <- NULL;
for(i in 1:10000){
    A_comp <- NULL;
    A_dat  <- rnorm(n = 256, mean = 0, sd = 0.4);
    A_mat  <- matrix(data = A_dat, nrow = 16);
    gammas <- runif(n = 16, min = 0, max = 1);
    mu_gam <- mean(gammas);
    diag(A_mat) <- -1;
    A1     <- gammas * A_mat;
    A0     <- mu_gam * A_mat;
    vA0    <- get_var_A(A0);
    vA1    <- get_var_A(A1);
    mA0    <- get_mn_A(A0);
    mA1    <- get_mn_A(A1);
    aA0    <- get_mn_abs_A(A0);
    aA1    <- get_mn_abs_A(A1);
    A0_e   <- eigen(A0)$values;
    A0_r   <- Re(A0_e);
    A0_i   <- Im(A0_e);
    A1_e   <- eigen(A1)$values;
    A1_r   <- Re(A1_e);
    A1_i   <- Im(A1_e);
    
    A0_s   <- max(A0_r) < 0;
    A1_s   <- max(A1_r) < 0;
    eigV0  <- cbind(A0_r, A0_i);
    eigV1  <- cbind(A1_r, A1_i); 
    sum_v  <- c(vA0, vA1, mA0, mA1, aA0, aA1, A0_s, A1_s);
    summ_eigs0 <- rbind(summ_eigs0, eigV0);
    summ_eigs1 <- rbind(summ_eigs1, eigV1);
    summ_stats <- rbind(summ_stats, sum_v);
    if(i %% 1000 == 0){
        print(i);    
    }
}

# But the radius of the eigen value distribution is bigger?



A_comp <- NULL;
A_dat  <- rnorm(n = 256, mean = 0, sd = 0.4);
A_mat  <- matrix(data = A_dat, nrow = 16);
gammas <- runif(n = 16, min = 0, max = 1);
mu_gam <- mean(gammas);
diag(A_mat) <- -1;
A1     <- gammas * A_mat;
A0     <- mu_gam * A_mat;
A0_e   <- eigen(A0)$values;
A0_r   <- Re(A0_e);
A0_i   <- Im(A0_e);
A1_e   <- eigen(A1)$values;
A1_r   <- Re(A1_e);
A1_i   <- Im(A1_e);

par(mar = c(5, 5, 1, 1));
plot(A0_r, A0_i, xlim = c(-2.75, 2.25), ylim = c(-2.5,2.5), pch = 4, cex = 0.7,
     xlab = "Real", ylab = "Imaginary", cex.lab = 1.5, cex.axis = 1.5)
points(A1_r, A1_i, col = "red", pch = 4, cex = 0.7)

vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(16) * sd(A0vec) * cos(vl) + mean(diag(A0));
y0 <- sqrt(16) * sd(A0vec) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3)


vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(16) * sd(A1vec) * cos(vl) + mean(diag(A1));
y0 <- sqrt(16) * sd(A1vec) * sin(vl);
points(x = x0, y = y0, type = "l", col = "red", lwd = 3)

# Bigger example









A_comp <- NULL;
A_dat  <- rnorm(n = 1000000, mean = 0, sd = 0.4);
A_mat  <- matrix(data = A_dat, nrow = 1000);
C_dat  <- rbinom(n = 1000 * 1000, size = 1, prob = 1);
C_mat  <- matrix(data = C_dat, nrow = 1000, ncol = 1000);
A_mat  <- A_mat * C_mat;
#gammas <- c(20, rep(1, 999));
#gammas <- c(rep(1.5, 500), rep(0.5, 500))
gammas <- runif(n = 1000, min = 0, max = 2);
mu_gam <- mean(gammas);
diag(A_mat) <- -1;
A1     <- gammas * A_mat;
A0     <- mu_gam * A_mat;
A0_e   <- eigen(A0)$values;
A0_r   <- Re(A0_e);
A0_i   <- Im(A0_e);
A1_e   <- eigen(A1)$values;
A1_r   <- Re(A1_e);
A1_i   <- Im(A1_e);

A0_vm       <- A0;
diag(A0_vm) <- NA;
A0vec       <- as.vector(A0_vm);
A0vec       <- A0vec[is.na(A0vec) == FALSE];
A1_vm       <- A1;
diag(A1_vm) <- NA;
A1vec       <- as.vector(A1_vm);
A1vec       <- A1vec[is.na(A1vec) == FALSE];

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 0, 0));
plot(A0_r, A0_i, xlim = c(-16.75, 16.25), ylim = c(-16.5,16.5), pch = 4, cex = 0.7,
     xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, asp = 1);
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(1000) * sd(A0vec) * cos(vl) + mean(diag(A0));
y0 <- sqrt(1000) * sd(A0vec) * sin(vl);
x1 <- sqrt(1000) * sd(A1vec) * cos(vl) + mean(diag(A1));
y1 <- sqrt(1000) * sd(A1vec) * sin(vl);
text(x = -19, y = 14, labels = "a", cex = 2);
points(x = x0, y = y0, type = "l", lwd = 3);
points(x = x1, y = y1, type = "l", col = "red", lwd = 3, lty = "dashed");
plot(A1_r, A1_i, xlim = c(-16.75, 16.25), ylim = c(-16.5,16.5), pch = 4, cex = 0.7,
     xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, asp = 1, col = "red",
     yaxt = "n");
text(x = -19, y = 14, labels = "b", cex = 2);
points(x = x1, y = y1, type = "l", col = "red", lwd = 3)
points(x = x0, y = y0, type = "l", lwd = 3, lty = "dashed");
mtext(side = 1, "Real", outer = TRUE, line = 3, cex = 2);
mtext(side = 2, "Imaginary", outer = TRUE, line = 2.5, cex = 2);




