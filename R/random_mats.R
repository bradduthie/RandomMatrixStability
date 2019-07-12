#' Find a stabilised system
#' 
#' Compares random matrices in which variation in component response rate does
#' not vary to random matrices in which this variation for 2 to max_sp 
#' components.
#'
#'@return A table of stability results, where rows summarise for each component
#'number (S) the number of stable or unstable (also, feasible and infeasible)
#'random matrices produced.
#'@param max_sp Maximum number of components to randomise
#'@param iters Number of iterations (i.e., random matrices) per component
#'@param int_type Type of interaction between components (0 is random)
#'@param rmx Standard deviation of non-zero matrix element components
#'@param C Connectedness of matrices (i.e., probability of non-zero matrix 
#'element components.
#'@examples
#'rand_gen_var(max_sp = 2, iters = 4);
#'@export
rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1, by = 1,
                         sigma = 0.4){
    tot_res <- NULL;
    fea_res <- NULL;
    sp_try  <- seq(from = by, to = max_sp, by = by);
    for(i in 1:length(sp_try)){
        iter           <- iters;
        tot_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        fea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        while(iter > 0){
            r_vec    <- rnorm(n = sp_try[i], mean = 0, sd = rmx);
            A0_dat   <- rnorm(n = sp_try[i] * sp_try[i], mean = 0, sd = sigma);
            A0       <- matrix(data = A0_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- species_interactions(mat = A0, type = int_type);
            C_dat    <- rbinom(n = sp_try[i] * sp_try[i], size = 1, prob = C);
            C_mat    <- matrix(data = C_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- A0 * C_mat;
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
    return(all_res);
}

#' Restrict component interaction types
#' 
#' Restricts the interaction types of a matrix to be competitive (1), mutualist
#' (2), or predator-prey (3).
#'
#'@return The matrix is returned with appropriately restricted interaction types
#'@param mat The matrix to be restricted
#'@param type The type of restriction being made
#'@examples
#'eg_mat  <- matrix(data = rnorm(n = 16), nrow = 4);
#'new_mat <- species_interactions(mat = eg_mat, type = 3);
#'@export
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

#' Summarise random matrix results
#' 
#' Takes the list output of the rand_gen_var function and summarises it into a
#' useable format.
#'
#'@return A table of stability and feasibility results.
#'@param tot_res The tot_res list output from the rand_gen_var function
#'@param fea_res The fea_res list output from the rand_gen_var function
#'@examples
#'eg_rand  <- rand_gen_var(max_sp = 2, iters = 4);
#'sum_rand <- summarise_randmat(eg_rand$tot_res, eg_rand$fea_res);
#'@export
summarise_randmat <- function(tot_res, fea_res){
    sims    <- length(tot_res);
    all_res <- matrix(data = 0, nrow = sims, ncol = 13);
    for(i in 1:sims){
        all_res[i, 1]  <- i + 1;
        # Stable and unstable
        all_res[i, 2]  <- sum(tot_res[[i]][,1] == FALSE);
        all_res[i, 3]  <- sum(tot_res[[i]][,1] == TRUE);
        all_res[i, 4]  <- sum(tot_res[[i]][,2] == FALSE);
        all_res[i, 5]  <- sum(tot_res[[i]][,2] == TRUE);
        # Stabilised and destabilised
        all_res[i, 6] <- sum(tot_res[[i]][,1] == FALSE & 
                                 tot_res[[i]][,2] == TRUE);
        all_res[i, 7] <- sum(tot_res[[i]][,1] == TRUE & 
                                 tot_res[[i]][,2] == FALSE);
        # Feasible and infeasible
        all_res[i, 8]  <- sum(fea_res[[i]][,1] == FALSE);
        all_res[i, 9]  <- sum(fea_res[[i]][,1] == TRUE);
        all_res[i, 10]  <- sum(fea_res[[i]][,2] == FALSE);
        all_res[i, 11]  <- sum(fea_res[[i]][,2] == TRUE);
        # Feased and defeased
        all_res[i, 12] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,2] == TRUE);
        all_res[i, 13] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,2] == FALSE);
    }
    cnames <- c("N", "A0_unstable", "A0_stable", "A1_unstable", "A1_stable", 
                "A1_stabilised", "A1_destabilised", "A0_infeasible", 
                "A0_feasible", "A1_infeasible", "A1_feasible", 
                "A1_made_feasible", "A1_made_infeasible");
    colnames(all_res) <- cnames;
    return(all_res);
}


mat_rho <- function(mat){
    if(dim(mat)[1] != dim(mat)[2]){
        stop("Error: Not a square matrix");
    }
    S  <- dim(mat)[1];
    od <- 0.5 * (S * S - S);
    ut <- rep(x = NA, times = od);
    lt <- rep(x = NA, times = od);
    ct <- 1;
    for(i in 1:dim(mat)[1]){
        for(j in 1:dim(mat)[2]){
            if(i < j){
                ut[ct] <- mat[i, j];
                lt[ct] <- mat[j, i]; 
                ct     <- ct + 1;
            }
        }
    }
    rho <- cor(ut, lt);
    return(rho);
}



rand_gen_rho <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1, by = 1,
                         sigma = 0.4){
    sp_try  <- seq(from = by, to = max_sp, by = by);
    ret_mat <- matrix(data = NA, nrow = length(sp_try) * iters, ncol = 8);
    ret_rho <- 1;
    for(i in 1:length(sp_try)){
        iter           <- iters;
        while(iter > 0){
            r_vec    <- rnorm(n = sp_try[i], mean = 0, sd = rmx);
            A0_dat   <- rnorm(n = sp_try[i] * sp_try[i], mean = 0, sd = sigma);
            A0       <- matrix(data = A0_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- species_interactions(mat = A0, type = int_type);
            C_dat    <- rbinom(n = sp_try[i] * sp_try[i], size = 1, prob = C);
            C_mat    <- matrix(data = C_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- A0 * C_mat;
            diag(A0) <- -1;
            gam1     <- runif(n = sp_try[i], min = 0, max = 2);
            A1       <- A0 * gam1;
            A0       <- A0 * mean(gam1);
            A0_stb   <- max(Re(eigen(A0)$values));
            A1_stb   <- max(Re(eigen(A1)$values));
            A0_rho   <- mat_rho(A0);
            A1_rho   <- mat_rho(A1);
            ret_mat[ret_rho, 1] <- sp_try[i];
            ret_mat[ret_rho, 2] <- iter;
            ret_mat[ret_rho, 3] <- A0_stb;
            ret_mat[ret_rho, 4] <- A1_stb;
            ret_mat[ret_rho, 5] <- as.numeric(A0_stb < 0);
            ret_mat[ret_rho, 6] <- as.numeric(A1_stb < 0);
            ret_mat[ret_rho, 7] <- A0_rho;
            ret_mat[ret_rho, 8] <- A1_rho;
            ret_rho             <- ret_rho + 1;
            iter                <- iter - 1;
        }
        print(sp_try[i]);
    }
    ret_mat <- ret_mat[ret_mat[,1] > 2,];
    return(ret_mat);
}
