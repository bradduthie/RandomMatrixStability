rand_gen_swn <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1, by = 4,
                         sigma = 0.4, Kdiv = 2, beta = 0.2){
    tot_res <- NULL;
    fea_res <- NULL;
    real_Cs <- NULL;
    sp_try  <- seq(from = by, to = max_sp, by = by);
    for(i in 1:length(sp_try)){
        iter           <- iters;
        tot_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        fea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
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
            iter    <- iter - 1;
        }
        print(sp_try[i]);
    }
    all_res     <- summarise_randmat(tot_res = tot_res, fea_res = fea_res);
    all_res[,1] <- sp_try;
    return(list(all_res = all_res, real_Cs = real_Cs));
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

sim <- rand_gen_swn(max_sp = 56, iters = 100000, rmx = 0.4, C = 1, by = 4,
                    sigma = 0.4, Kdiv = 2, beta = 0.0);




mat       <- matrix(data = 0, nrow = 2, ncol = 2);
mat[1, 2] <- 1;
mat[2, 1] <- 1;
for(N in 3:dim(mat)[1]){
    sum_Kj  <- sum(mat[1:N, 1:N]);
    sum_Kis <- apply(X = mat[1:N, 1:N], MARGIN = 1, FUN = sum);
    # sample bernoulli with pr of sum_Kis/sumkj in each new element
    # remember to flip.
}



