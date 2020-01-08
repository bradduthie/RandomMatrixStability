# THINGS TO DO FOR THE PAPER:
#
# 1. Do the cascade model, as suggested by reviewers, but to keep a broad focus
#    (i.e., not just on ecological networks) and do not do the niche or nested-
#    hierarchy models. The latter two are food webs, these two being refinements
#    of the cascade model. 
# 2. Do a small-world network, and note that many types of real (non-ecological)
#    networks are arranged this way. Will need to sample across wide range of
#    relevant parameters.
# 3. Do a scale-free network, and again note that many types of real networks
#    are arranged as such. Sample across a wide range of relevant parameters
# 
# The point of the above is that this is not an ecological paper per se, but is
# more broadly about the stability of complex networks given component response
# rate variation. So I do not want to get into precise ecological networks. I
# have already done predator-prey, mutualist, competitor, and now cascade 
# networks. And it is worth noting that looking for stability in these networks
# is less interesting that it would originally appear since, as noted in the
# paper, feasibility remains unaffected by varying component response rate.
#
# Lastly, I want to do the following:
#
# Vary the diagonal matrix (at least a bit). This was requested. I can show that
# the main point still works given a completely random matrix (show this), but
# might continue with random values selected for the diagonal around some sort
# of mean centre at a negative value (just to ensure some chance of stability).
#
# Analytically, I wonder if I can decompose M into the original circular matrix
# plus the contribution of \gamma. This could somehow show that \gamma causes
# variation in the real parts of the eigenvalues, which could sometimes lead to
# stability? Try this direction, and see where it leads at least
#



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
    rho_res <- NULL;
    cmplxty <- NULL;
    sp_try  <- seq(from = by, to = max_sp, by = by);
    for(i in 1:length(sp_try)){
        iter           <- iters;
        tot_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        fea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        rho_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 3);
        cmplxty[[i]]   <- matrix(data = 0, nrow = iter, ncol = 2);
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
            A0_rho   <- mat_rho(A0);
            A1_rho   <- mat_rho(A1);
            rho_diff <- A1_rho - A0_rho;
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
#'@param fea_res The rho_res list output from the rand_gen_var function
#'@examples
#'eg_rand  <- rand_gen_var(max_sp = 2, iters = 4);
#'sum_rand <- summarise_randmat(eg_rand$tot_res, eg_rand$fea_res);
#'@export
summarise_randmat <- function(tot_res, fea_res, rho_res, cmplxty){
    sims    <- length(tot_res);
    all_res <- matrix(data = 0, nrow = sims, ncol = 18);
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
        all_res[i, 17] <- mean(cmplxty[[i]][,1]);
        all_res[i, 18] <- mean(cmplxty[[i]][,2]);
    }
    cnames <- c("N", "A0_unstable", "A0_stable", "A1_unstable", "A1_stable", 
                "A1_stabilised", "A1_destabilised", "A0_infeasible", 
                "A0_feasible", "A1_infeasible", "A1_feasible", 
                "A1_made_feasible", "A1_made_infeasible", "A0_rho", "A1_rho",
                "rho_diff", "complex_A0", "complex_A1");
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

# Relate to zero law of biology? Variance in the correlation will inherently 
# increase, meaning that there is a higher probability that a negative
# correlation will occur and lead to stability. There is a bound at 1, but 
# multiplying by a vector makes more situations with a -1, assuming the
# initial correlation of M_{i,j} and M_{j,i} is not uniform?



## A1 multiplies by element -- still decreases stability.
rand_A1_test <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1, by = 1,
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
            G1_vals  <- runif(n = sp_try[i] * sp_try[i], min = 0, max = 2);
            G1       <- matrix(data = G1_vals, nrow = sp_try[i]);
            A1       <- A0 * G1;
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

# Build a matrix with fixed correlation between A[i, j] and A[j, i]
build_rho_mat <- function(S, sigma, rho){
    mat  <- matrix(data = 0, nrow = S, ncol = S);
    ia   <- 0.5 * ((S * S) - S);
    ut   <- rnorm(n = ia, mean = 0, sd = sigma);
    lt   <- rnorm(n = ia, mean = 0, sd = sigma);
    od   <- cbind(ut, lt);
    co   <- var(od);
    ch   <- solve(chol(co));
    nx   <- od %*% ch;
    ms   <- matrix(data = c(1, rho, rho, 1), ncol = 2);
    c2   <- chol(ms);
    el   <- nx %*% c2 * sd(ut) + mean(ut);
    fc   <- 1;
    for(i in 1:S){
        for(j in 1:S){
            if(i < j){
                mat[i, j] <- el[fc, 1];
                mat[j, i] <- el[fc, 2];
                fc        <- fc + 1;
            }
        }
    }
    diag(mat)           <- -1;
    return(mat);
}
