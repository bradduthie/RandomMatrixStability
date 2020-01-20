rand_gen_swn <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1, by = 4,
                         sigma = 0.4, Kdiv = 2, beta = 0.2){
    tot_res <- NULL;
    fea_res <- NULL;
    real_Cs <- NULL;
    rho_res <- NULL;
    cmplxty <- NULL;
    sp_try  <- seq(from = by, to = max_sp, by = by);
    for(i in 1:length(sp_try)){
        iter           <- iters;
        tot_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        fea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        rho_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 4);
        cmplxty[[i]]   <- matrix(data = 0, nrow = iter, ncol = 2);
        real_Cs[[i]]   <- matrix(data = 0, nrow = iter, ncol = 3);
        while(iter > 0){
            r_vec    <- rnorm(n = sp_try[i], mean = 0, sd = rmx);
            A0_dat   <- rnorm(n = sp_try[i] * sp_try[i], mean = 0, sd = sigma);
            A0       <- matrix(data = A0_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- species_interactions(mat = A0, type = int_type);
            swn      <- create_swn(N = sp_try[i], beta = beta,
                                   K = (sp_try[i] / Kdiv));
            real_Cs[[i]][iter, 1] <- sp_try[i];
            real_Cs[[i]][iter, 2] <- iter;
            real_Cs[[i]][iter, 3] <- get_C(swn);
            C_dat    <- rbinom(n = sp_try[i] * sp_try[i], size = 1, prob = C);
            C_mat    <- matrix(data = C_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- A0 * C_mat;
            A0       <- A0 * swn;
            diag(A0) <- -1;
            gam1     <- runif(n = sp_try[i], min = 0, max = 2);
            A1       <- A0 * gam1;
            A0       <- A0 * mean(gam1);
            A0_stb   <- max(Re(eigen(A0)$values)) < 0;
            A1_stb   <- max(Re(eigen(A1)$values)) < 0;
            A0_fea   <- min(-1*solve(A0) %*% r_vec) > 0;
            A1_fea   <- min(-1*solve(A1) %*% r_vec) > 0;
            A0_rho   <- mat_rho(A0);
            A1_rho   <- mat_rho(A1);
            rho_diff <- A1_rho - A0_rho;
            rho_abs  <- abs(A1_rho) - abs(A0_rho);
            if(A0_stb == TRUE){
                tot_res[[i]][iter, 1] <- 1;
            }
            if(A1_stb == TRUE){
                tot_res[[i]][iter, 2] <- 1;
            }
            if(A0_fea == TRUE){
                fea_res[[i]][iter, 1] <- 1;
            }
            if(A1_fea == TRUE){
                fea_res[[i]][iter, 2] <- 1;
            }
            rho_res[[i]][iter, 1] <- A0_rho;
            rho_res[[i]][iter, 2] <- A1_rho;
            rho_res[[i]][iter, 3] <- rho_diff;
            rho_res[[i]][iter, 4] <- rho_abs;
            cmplxty[[i]][iter, 1] <- get_complexity(A0);
            cmplxty[[i]][iter, 2] <- get_complexity(A1);
            iter                  <- iter - 1;
        }
        print(sp_try[i]);
    }
    all_res     <- summarise_randmat(tot_res = tot_res, fea_res = fea_res,
                                     rho_res = rho_res, cmplxty = cmplxty);
    all_res[,1] <- sp_try;
    full_res    <- list(all_res = all_res, real_Cs = real_Cs);
    res_table   <- add_C_stats(sim = full_res);
    return(res_table);
}

add_C_stats <- function(sim){
    real_Cs <- sim$real_Cs;
    all_res <- sim$all_res;
    res_mns <- matrix(data = 0, nrow = length(sim$real_Cs), ncol = 2);
    for(i in 1:dim(res_mns)[1]){
        mn_vals       <- apply(X = real_Cs[[i]], MARGIN = 2, FUN = mean);
        res_mns[i, 1] <- mn_vals[1];
        res_mns[i, 2] <- mn_vals[3];
    }
    new_all_res <- cbind(all_res, res_mns[,2]);
    colnames(new_all_res)[dim(new_all_res)[2]] <- "C";
    return(new_all_res);
}

get_complexity <- function(mat){
    S         <- length(diag(mat));
    mat_gz    <- sum(mat != 0);
    trace_gz  <- sum(diag(mat) != 0);
    offdiagC  <- mat_gz - trace_gz;
    offdiags  <- (dim(mat)[1] * dim(mat)[2]) - S;
    calc_C    <- offdiagC / offdiags;
    diag(mat) <- NA;
    sigma     <- sd(mat, na.rm = TRUE);
    complx    <- sigma * sqrt(S * calc_C);
    return(complx);
}

create_swn <- function(N = 100, K = 20, beta = 0.05){
    mat <- matrix(data = 0, nrow = N, ncol = N);
    diag(mat) <- 1;
    if(K %% 2 != 0){
        warning("K needs to be an even integer");
    }
    if(K > dim(mat)[1]){
        warning("K is too high");
    }
    mat <- add_sw_edges(mat, K);
    mat <- rewire_sw_edges(mat, K, beta);
    return(mat);
}

add_sw_edges <- function(mat, K){
    N <- dim(mat)[1];
    for(i in 1:N){
        for(j in 1:N){
            if(i > j){
                neig1 <- abs(i - j);
                neig2 <- abs(i - (N+j));
                neig  <- min(c(neig1, neig2));
                if(0 < neig & neig <= (K/2)){
                    mat[i, j] <- 1;
                    mat[j, i] <- 1;
                }
            }
       }
    }
    return(mat);
}

rewire_sw_edges <- function(mat, K, beta){
    N <- dim(mat)[1];
    for(i in 1:N){
        for(j in 1:N){
            if(i != j & mat[i, j] > 0){
                samp_beta <- runif(n = 1) < beta;
                if(samp_beta){ # Rewire node
                    mat[i, j] <- 0;
                    mat[j, i] <- 0;
                    k         <- sample_k(i, N, K);
                    mat[i, k] <- 1;
                    mat[k, i] <- 1;
                }
            }
        }
    }
    return(mat);
}

sample_k <- function(i, N, K){
    tsamp <- 1:N;
    lo    <- i - (K/2);
    hi    <- i + (K/2);
    span1 <- seq(from = lo, to = hi, by = 1);
    span2 <- span1 + N;
    span  <- c(span1, span2);
    nouse <- span[span %in% tsamp];
    tsamp <- tsamp[-nouse];
    getit <- sample(tsamp, size = 1);
    return(getit);
}

get_C <- function(mat){
    c1  <- sum(mat != 0);
    tot <- (dim(mat)[1] * dim(mat)[2]) - length(diag(mat));
    return(c1 / tot);
}

visualise_network <- function(mat){
    N  <- dim(mat)[1];
    rd <- seq(from = 0, to = 2*pi, length = N);
    yy <- sin(rd);
    xx <- cos(rd);
    par(bty = "n");
    plot(x = xx, y = yy, pch = 20, cex = 1.5, 
         xaxt = "n", yaxt = "n", xlab = "", ylab = "");
    for(i in 1:N){
        for(j in 1:N){
            if(i > j & mat[i,j] > 0){
                lines(x = c(xx[i], xx[j]), y = c(yy[i], yy[j]), 
                      lwd = 0.5);
            }
        }
    }
}




