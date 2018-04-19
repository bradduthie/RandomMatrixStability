species_interactions <- function(mat, type = 0){
    if(type == 1){
        mat[mat > 0] <- -1*mat[mat > 0];
    }
    if(type == 2){
        mat[mat < 0] <- -1*mat[mat < 0];
    }
    if(type == 3){
        for(i in 1:dim(mat)[1]){
            for(j in 1:dim(mat)[2]){
                if(mat[i, j] * mat[j, i] > 0){
                    mat[j, i] <- -1 * mat[j, i];
                }
            }
        }
    }
    return(mat);
}

make_gammas <- function(nn = 10, distribution = 1, mn = 1, sdd = 1){
    if(distribution == 0){
        dat          <- rep(x = mn, times = nn);
    }
    if(distribution == 1){
        mval         <- sdd * sqrt(12);
        dat          <- runif(n = nn, min = 0, max = mval);
    }
    if(distribution == 2){
        dat          <- rexp(n = nn);
        dat          <- sdd * ((dat - mean(dat)) / sd(dat));
        dat          <- dat - min(dat) + min(abs(dat));
    }
    if(distribution == 3){
        dat          <- rbeta(n = nn, shape1 = 0.5, shape2 = 0.5);
        dat          <- sdd * ((dat - mean(dat)) / sd(dat));
        dat          <- dat - min(dat) + min(abs(dat));
    }
    if(distribution == 4){
        dat          <- rgamma(n = nn, shape = 2, scale = 2);
        dat          <- sdd * ((dat - mean(dat)) / sd(dat));
        dat          <- dat - min(dat) + min(abs(dat));
    }
    return(dat);
}

summarise_randmat <- function(tot_res, fea_res){
    sims    <- length(tot_res);
    all_res <- matrix(data = 0, nrow = sims, ncol = 37);
    for(i in 1:sims){
        all_res[i, 1]  <- i + 1;
        # Stable and unstable
        all_res[i, 2]  <- sum(tot_res[[i]][,1] == FALSE);
        all_res[i, 3]  <- sum(tot_res[[i]][,1] == TRUE);
        all_res[i, 4]  <- sum(tot_res[[i]][,2] == FALSE);
        all_res[i, 5]  <- sum(tot_res[[i]][,2] == TRUE);
        all_res[i, 6]  <- sum(tot_res[[i]][,3] == FALSE);
        all_res[i, 7]  <- sum(tot_res[[i]][,3] == TRUE);
        all_res[i, 8]  <- sum(tot_res[[i]][,4] == FALSE);
        all_res[i, 9]  <- sum(tot_res[[i]][,4] == TRUE);
        all_res[i, 10] <- sum(tot_res[[i]][,5] == FALSE);
        all_res[i, 11] <- sum(tot_res[[i]][,5] == TRUE);
        # Stabilised and destabilised
        all_res[i, 12] <- sum(tot_res[[i]][,1] == FALSE & 
                                  tot_res[[i]][,2] == TRUE);
        all_res[i, 13] <- sum(tot_res[[i]][,1] == FALSE & 
                                  tot_res[[i]][,3] == TRUE);
        all_res[i, 14] <- sum(tot_res[[i]][,1] == FALSE & 
                                  tot_res[[i]][,4] == TRUE);
        all_res[i, 15] <- sum(tot_res[[i]][,1] == FALSE & 
                                  tot_res[[i]][,5] == TRUE);
        all_res[i, 16] <- sum(tot_res[[i]][,1] == TRUE & 
                                  tot_res[[i]][,2] == FALSE);
        all_res[i, 17] <- sum(tot_res[[i]][,1] == TRUE & 
                                  tot_res[[i]][,3] == FALSE);
        all_res[i, 18] <- sum(tot_res[[i]][,1] == TRUE & 
                                  tot_res[[i]][,4] == FALSE);
        all_res[i, 19] <- sum(tot_res[[i]][,1] == TRUE & 
                                  tot_res[[i]][,5] == FALSE);
        # Feasible and infeasible
        all_res[i, 20]  <- sum(fea_res[[i]][,1] == FALSE);
        all_res[i, 21]  <- sum(fea_res[[i]][,1] == TRUE);
        all_res[i, 22]  <- sum(fea_res[[i]][,2] == FALSE);
        all_res[i, 23]  <- sum(fea_res[[i]][,2] == TRUE);
        all_res[i, 24]  <- sum(fea_res[[i]][,3] == FALSE);
        all_res[i, 25]  <- sum(fea_res[[i]][,3] == TRUE);
        all_res[i, 26]  <- sum(fea_res[[i]][,4] == FALSE);
        all_res[i, 27]  <- sum(fea_res[[i]][,4] == TRUE);
        all_res[i, 28]  <- sum(fea_res[[i]][,4] == FALSE);
        all_res[i, 29]  <- sum(fea_res[[i]][,4] == TRUE);
        # Feased and defeased
        all_res[i, 30] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,2] == TRUE);
        all_res[i, 31] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,3] == TRUE);
        all_res[i, 32] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,4] == TRUE);
        all_res[i, 33] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,5] == TRUE);
        all_res[i, 34] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,2] == FALSE);
        all_res[i, 35] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,3] == FALSE);
        all_res[i, 36] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,4] == FALSE);
        all_res[i, 37] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,5] == FALSE);  
    }
    cnames <- c("N", "A0_unstable", "A0_stable", "A1_unstable", "A1_stable", 
                "A2_unstable", "A2_stable", "A3_unstable", "A3_stable", 
                "A4_unstable", "A4_stable", "A1_stabilised", "A2_stabilised",
                "A3_stabilised", "A4_stabilised", "A1_destabilised",
                "A2_destabilised", "A3_destabilised", "A4_destabilised", 
                "A0_infeasible", "A0_feasible", "A1_infeasible", "A1_feasible", 
                "A2_infeasible", "A2_feasible", "A3_infeasible", "A3_feasible", 
                "A4_infeasible", "A4_feasible", "A1_made_feasible",
                "A2_made_feasible", "A3_made_feasible", "A4_made_feasible",
                "A1_made_infeasible",  "A2_made_infeasible", 
                "A3_made_infeasible", "A4_made_infeasible");
    colnames(all_res) <- cnames;
    return(all_res);
}

rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, gamma_sd = 1, 
                         gamma_mn = 1, C = 1){
    tot_res <- NULL;
    fea_res <- NULL;
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
            iter    <- iter - 1;
        }
        print(i);
    }
    all_res <- summarise_randmat(tot_res = tot_res, fea_res = fea_res);
    return(all_res);
}
