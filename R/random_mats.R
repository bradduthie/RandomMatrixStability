#' Find stabilised random systems
#' 
#' Compares random matrices in which variation in component response rate does
#' not vary to random matrices in which this variation for 2 to max_sp 
#' components. This function is to be used for random networks, in which
#' off-diagonal elements are randomly sampled from a normal distribution with
#' a mean of 'mn' and a standard deviation of 'sigma'.
#'
#'@return A table of stability results, where rows summarise for each component
#'number (S) the number of stable or unstable (also, feasible and infeasible)
#'random matrices produced.
#'@param max_sp Maximum number of components to randomise
#'@param by Sequence between component number to randomise (e.g., '2': 2, 4, 6)
#'@param iters Number of iterations (i.e., random matrices) per component
#'@param int_type Type of interaction between components including random (0),
#'competitor (1), mutualist (2), predator-prey (3), and cascade model (4)
#'@param rmx Standard deviation of population growth rates (for feasibility)
#'@param C Connectedness of matrices (i.e., probability of non-zero matrix 
#'element components.
#'@param sigma Standard deviation of interaction strength among network elements
#'@param mn Mean interaction strength among network elements
#'@param dval Self-regulation of network elements (1 by default)
#'@examples
#'rand_gen_var(max_sp = 2, iters = 4);
#'@export
rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1, by = 1,
                         sigma = 0.4, mn = 0, dval = 1){
    tot_res <- NULL;
    fea_res <- NULL;
    rho_res <- NULL;
    cmplxty <- NULL;
    sp_try  <- seq(from = by, to = max_sp, by = by);
    for(i in 1:length(sp_try)){
        iter           <- iters;
        tot_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        fea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        rho_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 6);
        cmplxty[[i]]   <- matrix(data = 0, nrow = iter, ncol = 2);
        while(iter > 0){
            r_vec    <- rnorm(n = sp_try[i], mean = mn, sd = rmx);
            A0_dat   <- rnorm(n = sp_try[i] * sp_try[i], mean = mn, sd = sigma);
            A0       <- matrix(data = A0_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- species_interactions(mat = A0, type = int_type);
            C_dat    <- rbinom(n = sp_try[i] * sp_try[i], size = 1, prob = C);
            C_mat    <- matrix(data = C_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- A0 * C_mat;
            diag(A0) <- -1 * dval;
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
            rho_res[[i]][iter, 5] <- max(Re(eigen(A0)$values));
            rho_res[[i]][iter, 6] <- max(Re(eigen(A1)$values));
            cmplxty[[i]][iter, 1] <- get_complexity(A0);
            cmplxty[[i]][iter, 2] <- get_complexity(A1);
            iter                  <- iter - 1;
        }
        print(sp_try[i]);
    }
    all_res     <- summarise_randmat(tot_res = tot_res, fea_res = fea_res,
                                     rho_res = rho_res, cmplxty = cmplxty);
    all_res[,1] <- sp_try;
    return(all_res);
}

#' Restrict component interaction types
#' 
#' Restricts the interaction types of a matrix to be competitive (1), mutualist
#' (2), predator-prey (3), or cascade model (4).
#'
#'@return The matrix is returned with appropriately restricted interaction types
#'@param mat The matrix to be restricted
#'@param type The type of restriction being made. This includes random (0),
#'competitor (1), mutualist (2), predator-prey (3), or cascade model food web 
#'(4)
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
    if(type == 4){
        for(i in 1:dim(mat)[1]){
            for(j in 1:dim(mat)[2]){
                if(j > i & mat[i, j] < 0){
                    mat[i, j] <- -1 * mat[i, j];
                }
                if(j < i & mat[i, j] > 0){
                    mat[i, j] <- -1 * mat[i, j];
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
#'@param rho_res The rho_res list output from the rand_gen_var function
#'@param cmplexty The cmplxty list output from the rand_gen_var function
#'@examples
#'eg_rand  <- rand_gen_var(max_sp = 2, iters = 4);
#'sum_rand <- summarise_randmat(eg_rand$tot_res, eg_rand$fea_res);
#'@export
summarise_randmat <- function(tot_res, fea_res, rho_res, cmplxty){
    sims    <- length(tot_res);
    all_res <- matrix(data = 0, nrow = sims, ncol = 25);
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
        # Mean correlations between A[i,j] and A[j,i] off-diag elements
        all_res[i, 14] <- mean(rho_res[[i]][,1]);
        all_res[i, 15] <- mean(rho_res[[i]][,2]);
        all_res[i, 16] <- mean(rho_res[[i]][,3]);
        all_res[i, 17] <- mean(rho_res[[i]][,4]);
        all_res[i, 18] <- mean(cmplxty[[i]][,1]);
        all_res[i, 19] <- mean(cmplxty[[i]][,2]);
        all_res[i, 20] <- mean(rho_res[[i]][,5]);
        all_res[i, 21] <- mean(rho_res[[i]][,6]);
        all_res[i, 22:23] <- range(rho_res[[i]][,5]);
        all_res[i, 24:25] <- range(rho_res[[i]][,6]);
    }
    cnames <- c("N", "A0_unstable", "A0_stable", "A1_unstable", "A1_stable", 
                "A1_stabilised", "A1_destabilised", "A0_infeasible", 
                "A0_feasible", "A1_infeasible", "A1_feasible", 
                "A1_made_feasible", "A1_made_infeasible", "A0_rho", "A1_rho",
                "rho_diff", "rho_abs", "complex_A0", "complex_A1", "A0_eig",
                "A1_eig", "LCI_A0", "UCI_A0", "LCI_A1", "UCI_A1");
    colnames(all_res) <- cnames;
    return(all_res);
}

#' Find a stabilised system for correlated systems
#' 
#' Compares random matrices in which variation in component response rate does
#' not vary to random matrices in which this variation for a fixed number of S
#' components. This function is to be used for random networks, in which
#' off-diagonal elements are randomly sampled from a normal distribution with
#' a mean of 'mn', a standard deviation of 'sigma', and a pre-specified
#' set of correlation structures.
#'
#'@return A table of stability results, where rows summarise for each 
#'correlation (rho) the number of stable or unstable (also, feasible and 
#'infeasible) random matrices produced.
#'@param S Number of components to randomise
#'@param rhos Vector of correlations to simulate
#'@param iters Number of iterations (i.e., random matrices) per component
#'@param rmx Standard deviation of population growth rates (for feasibility)
#'@param C Connectedness of matrices (i.e., probability of non-zero matrix 
#'element components.
#'@param sigma Standard deviation of interaction strength among network elements
#'@param mn Mean interaction strength among network elements
#'@param dval Self-regulation of network elements (1 by default)
#'@examples
#'rand_rho_var(S = 10, rhos = c(-0.2, 0, 0.2), iters = 4);
#'@export
rand_rho_var <- function(S, rhos, iters, int_type = 0, rmx = 0.4, C = 1, 
                         sigma = 0.4, mn = 0, dval = 1){
    tot_res <- NULL;
    fea_res <- NULL;
    rho_res <- NULL;
    cmplxty <- NULL;
    for(i in 1:length(rhos)){
        iter           <- iters;
        tot_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        fea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        rho_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 6);
        cmplxty[[i]]   <- matrix(data = 0, nrow = iter, ncol = 2);
        while(iter > 0){
            r_vec    <- rnorm(n = S, mean = 0, sd = rmx);
            A0       <- build_rho_mat(S = S, sigma = sigma, rho = rhos[i], 
                                      mn = mn, dval = dval);
            gam1     <- runif(n = S, min = 0, max = 2);
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
            rho_res[[i]][iter, 5] <- max(Re(eigen(A0)$values));
            rho_res[[i]][iter, 6] <- max(Re(eigen(A1)$values));
            cmplxty[[i]][iter, 1] <- get_complexity(A0);
            cmplxty[[i]][iter, 2] <- get_complexity(A1);
            iter                  <- iter - 1;
        }
        print(rhos[i]);
    }
    all_res     <- summarise_randmat(tot_res = tot_res, fea_res = fea_res,
                                     rho_res = rho_res, cmplxty = cmplxty);
    all_res[,1]          <- rhos;
    colnames(all_res)[1] <- "rho_set";
    return(all_res);
}


