#' Find stabilised small world networks
#' 
#' Compares small world networks in which variation in component response rate 
#' does not vary to random matrices in which this variation for 2 to max_sp 
#' components. This function is to be used for small world networks, which
#' are constructed using the method of Watts and Strogatz (Nature vol. 393,
#' 1998). Off-diagonal elements have a mean of 'mn' and a standard deviation of
#' 'sigma'.
#'
#'@return A table of stability results, where rows summarise for each component
#'number (S) the number of stable or unstable (also, feasible and infeasible)
#'random matrices produced.
#'@param max_sp Maximum number of components to randomise
#'@param iters Number of iterations (i.e., random matrices) per component
#'@param int_type Type of interaction between components including random (0),
#'competitor (1), mutualist (2), predator-prey (3), and cascade model (4)
#'@param rmx Standard deviation of population growth rates (for feasibility)
#'@param C Connectedness of matrices (i.e., probability of non-zero matrix 
#'element components.
#'@param sigma Standard deviation of interaction strength among network elements
#'@param beta Probability that a random interaction in a regular network is
#'rewired (parameter p in Watts and Strogatz 1998)
#'@param Kdiv Value to divide the component number by to produce the parameter
#'K for creating the small world network. For example, if S = 32 and K = 8, then
#'the small world network will be created from a regular network in which each
#'component is connected to 32/8 = 4 other components. This needs to be used
#'cautiously to avoid generating non-even values of K.
#'@param mn Mean interaction strength among network elements
#'@param dval Self-regulation of network elements (1 by default)
#'@examples
#'rand_gen_swn(max_sp = 16, iters = 4);
#'@export
rand_gen_swn <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1, by = 4,
                         sigma = 0.4, Kdiv = 2, beta = 0.2, mn = 0, dval = 1){
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
            A0_dat   <- rnorm(n = sp_try[i] * sp_try[i], mean = mn, sd = sigma);
            A0       <- matrix(data = A0_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- species_interactions(mat = A0, type = int_type);
            swn      <- create_swn(S = sp_try[i], beta = beta,
                                   K = (sp_try[i] / Kdiv));
            real_Cs[[i]][iter, 1] <- sp_try[i];
            real_Cs[[i]][iter, 2] <- iter;
            real_Cs[[i]][iter, 3] <- get_C(swn);
            C_dat    <- rbinom(n = sp_try[i] * sp_try[i], size = 1, prob = C);
            C_mat    <- matrix(data = C_dat, nrow = sp_try[i], 
                               ncol = sp_try[i]);
            A0       <- A0 * C_mat;
            A0       <- A0 * swn;
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

#' Build a small world network
#' 
#' Builds a small world network using the method of Watts and Strogatz (Nature 
#' vol. 393, 1998)
#'
#'@return A small world network represented by a square matrix
#'@param S The size of the network (number of components)
#'@param beta Probability that a random interaction in a regular network is
#'rewired (parameter p in Watts and Strogatz 1998)
#'@param K Number of edges that each vertice in the network contains
#'@examples
#'eg_swn <- create_swn(S = 32, K = 4);
#'@export
create_swn <- function(S = 100, K = 20, beta = 0.05){
    mat <- matrix(data = 0, nrow = S, ncol = S);
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