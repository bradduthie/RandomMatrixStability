#' Find stabilised scale-free world networks
#' 
#' Compares scale-free networks in which variation in component response rate 
#' does not vary to random matrices in which this variation for 'm' to 'max_sp' 
#' components. This function is to be used for scale-free networks, which
#' are constructed using the method of Albert and Barab\'{a}si (Reviews of 
#' Modern Physics, Vo. 74, 2002). Off-diagonal elements have a mean of 'mn' and 
#' a standard deviation of 'sigma'.
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
#'@param m The number of vertices that newly added components will have in 
#'scale-free network construction
#'@param mn Mean interaction strength among network elements
#'@param dval Self-regulation of network elements (1 by default)
#'@examples
#'rand_gen_sfn(max_sp = 16, iters = 4);
#'@export
rand_gen_sfn <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1, by = 4,
                         sigma = 0.4, m = 2, mn = 0, dval = 1){
    tot_res <- NULL;
    fea_res <- NULL;
    real_Cs <- NULL;
    rho_res <- NULL;
    cmplxty <- NULL;
    start_S <- max(c(by, m));
    sp_try  <- seq(from = start_S, to = max_sp, by = by);
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
            swn      <- create_sfn(S = sp_try[i], m = m);
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

#' Build a scale-free network
#' 
#' Builds a scale-free network using the method of Albert and Barab\'{a}si
#' (Reviews of Modern Physics, Vo. 74, 2002). A saturated network of size 'm' is
#' first seeded, with new components being added with 'm' edges sequentially
#' until a network of size S is produced.
#' 
#'@return A scale-free network represented by a square matrix
#'@param S The size of the network (number of components)
#'@param m The number of vertices that newly added components will have
#'@examples
#'eg_swn <- create_sfn(S = 8, m = 2);
#'@export
create_sfn <- function(S, m){
    mat <- matrix(data = 0, nrow = S, ncol = S);
    mat <- seed_mat(mat, m);
    for(i in m:S){
        new_p <- get_p(mat = mat[1:i, 1:i], m = m);
        mat   <- add_p(mat = mat, new_p = new_p);
    }
    return(mat);
}

seed_mat <- function(mat, m){
    mat[1:m, 1:m]  <- 1;
    diag(mat)[1:m] <- 0;
    return(mat);
}

get_p <- function(mat, m){ # This function gets links for the new node
    msize    <- dim(mat)[1];
    edges    <- apply(X = mat, MARGIN = 1, FUN = sum);
    prs      <- edges / sum(edges);
    ch       <- sample(x = 1:msize, size = m, prob = prs);
    newe     <- rep(x = 0, times = msize);
    newe[ch] <- 1;
    return(newe);
}

add_p <- function(mat, new_p){
    row             <- length(new_p);
    mat[row, 1:row] <- new_p;
    mat[1:row, row] <- new_p;
    return(mat);
}