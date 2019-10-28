# NOTE FOR RESEARCH 27 OCT 2019
# SEE EQN 2.4-2.9 OF AHMADIAN ET AL 2015
# PARTICULARLY 2.8, IN WHICH THE SPECTRAL DENSITY IS DEPENDENT UP ON THE MATRIX
# L. THE DISTRIBUTION APPEARS TO DEPEND ON L^-1 IN A WAY (TRACE) THAT MIGHT 
# CAUSE THE STABILITY. NEED TO FIGURE OUT HOW THIS WORKS FOR N < 10 AND N > 10?
N  <- 4000;
z  <- 10;
M  <- matrix(data = rnorm(N * N, sd = 1/N), nrow = N, ncol = N);
Mz <- z - M;
ml <- Mz %*% t(Mz);
ms <- solve(ml);
tr <- (1/N) * sum(diag(ms));
tr;


Rg <- Re(eigen(M)$values);
Ig <- Im(eigen(M)$values);
plot(x = Rg, y = Ig);

# Maybe instead, define L as {N, N, N, N}, so the SD is preserved? Then
# use the tools above?

# Another check: Show unusually high values can be shunted into the very
# negative group.
N <- 800;
L <- matrix(data = 0, nrow = N, ncol = N);
J <- matrix(data = rnorm(N*N, sd = sqrt(1/N)), nrow = N, ncol = N);
diag(J)  <- 1;
diag(L)  <- -1 * c(rep(x = 0.01, times = N/2), rep(x = 1.98, times = N/2));
M1       <- L %*% J
eM1      <- eigen(M1)$values;
rM1      <- Re(eM1);
iM1      <- Im(eM1);
plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", xlim = c(-0.015, 0.015),
     ylim = c(-0.015, 0.015));
abline(v = 0, col = "red");
plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-1.5, -1.5), ylim = c(-1.5, 1.5));
rm1 <- tail(sort(rM1));

L        <- matrix(data = 0, nrow = N, ncol = N);
diag(L)  <- -0.995;
M2       <- L %*% J;
eM2      <- eigen(M2)$values;
rM2      <- Re(eM2);
iM2      <- Im(eM2);
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", xlim = c(-1.5, 1.5),
     ylim = c(-1.5, 1.5));
abline(v = 0, col = "red");
rm2 <- tail(sort(rM2));

rm1;
rm2;

plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "x", 
     xlim = c(-1.5, -1.5), ylim = c(-1.5, 1.5));
points(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "x", xlim = c(-1.5, 1.5),
     ylim = c(-1.5, 1.5), col = "blue");
abline(v = 0, col = "red");


plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "x", xlim = c(-0.15, 0.15),
     ylim = c(-0.15, 0.15));
points(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "x", col = "blue");
abline(v = 0, col = "red");

################################################################################
################################################################################
################################################################################


# Let's compare the eigenvalue distributions of the upper-left of M1, versus
# the whole eigenvalues to see what the rest adds.
M_ul  <- M1[1:400, 1:400];
em_ul <- eigen(M_ul)$values;
r_ul  <- Re(em_ul);
i_ul  <- Im(em_ul);

# Zoom in
plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02));
abline(v = 0, col = "red");
points(x = r_ul, y = i_ul, asp = 1, cex = 0.8, pch = "x", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02), col = "orange");
# The cloud is clearly bigger for the eigenvalues of the full matrix.

# Now add in the rest of it.
M_ll  <- M1[401:800, 401:800];
em_ll <- eigen(M_ll)$values;
r_ll  <- Re(em_ll);
i_ll  <- Im(em_ll);

plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-3.5, 0.5), ylim = c(-2, 2));
abline(v = 0, col = "red");
points(x = r_ul, y = i_ul, asp = 1, cex = 0.8, pch = "x", 
       xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02), col = "orange");
points(x = r_ll, y = i_ll, asp = 1, cex = 0.8, pch = "x", 
       xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02), col = "orange");
# The size of the bigger cloud does not appear to have changed much.
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(1/N) * 0.01 * sqrt(N/2) * cos(vl) + mean(diag(M_ul));
y0 <- sqrt(1/N) * 0.01 * sqrt(N/2) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3, col = "orange");

################################################################################
################################################################################
################################################################################
# SOME PROGRESS
################################################################################
################################################################################
################################################################################
# The following code will make a random complex system with the following
# parameter values:
# S       = 800
# sigma^2 = 0.1870829
# sigma   = 0.035
# C       = 1
# Make the random matrix with no variation in gamma
S        <- 800;
sigma    <- 0.036;
M        <- matrix(data = rnorm(S*S, sd = sigma), nrow = S, ncol = S);
diag(M)  <- -1;
M1       <- M;
L        <- matrix(data = 0, nrow = S, ncol = S);
diag(L)  <- rep(x = 1, times = S);
eM1      <- eigen(M1)$values;
rM1      <- Re(eM1);
iM1      <- Im(eM1);
plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", xlim = c(-2.5, 0.5),
     ylim = c(-1.5, 1.5), xlab = "Real", ylab = "Imaginary", cex.lab = 1.25);
abline(v = 0, col = "red");

# Now make the one with half slow components and half fast components
M2       <- M;
L        <- matrix(data = 0, nrow = S, ncol = S);
diag(L)  <- c(rep(x = 0.01, times = S/2), rep(x = 1.99, times = S/2));
diag(M2) <- -1;
M2       <- L %*% M2;
eM2      <- eigen(M2)$values;
rM2      <- Re(eM2);
iM2      <- Im(eM2);
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", xlim = c(-3.5, 0.5),
     ylim = c(-2, 2));
abline(v = 0, col = "red");
# Zoom in on that small circle
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02));
abline(v = 0, col = "red");

# See what the eigenvalue distributions would be if we only had two matrices
# Upper left, and lower right.
M3       <- M2[1:(S/2), 1:(S/2)];
eM3      <- eigen(M3)$values;
rM3      <- Re(eM3);
iM3      <- Im(eM3);
# Look at that circle
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02),
     xlab = "Real", ylab = "Imaginary", cex.lab = 1.25);
abline(v = 0, col = "red");
points(x = rM3, y = iM3, asp = 1, cex = 0.8, pch = "x", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02), col = "orange");
# This corresponds to the expected diameter
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sigma * 0.01 * sqrt(N/2) * cos(vl) + mean(diag(M3));
y0 <- sigma * 0.01 * sqrt(N/2) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3, col = "orange");
# Check out what happens when you change N
x0 <- sigma * 0.01 * sqrt(N) * cos(vl) + mean(diag(M3));
y0 <- sigma * 0.01 * sqrt(N) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 2, col = "orange", lty = "dashed");


# See what both look like if separate matrices, upper left and lower right.
M4       <- M2[(S/2 + 1):S, (S/2 + 1):S];
eM4      <- eigen(M4)$values;
rM4      <- Re(eM4);
iM4      <- Im(eM4);
# Look at that circle
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", xlim = c(-3.5, 0.5),
     ylim = c(-2, 2), xlab = "Real", ylab = "Imaginary", cex.lab = 1.25);
abline(v = 0, col = "red");
# Add the old M3 ones
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sigma * 0.01 * sqrt(N/2) * cos(vl) + mean(diag(M3));
y0 <- sigma * 0.01 * sqrt(N/2) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3, col = "orange");
# Add the new M4 ones
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x1 <- sigma * 1.99 * sqrt(N/2) * cos(vl) + mean(diag(M4));
y1 <- sigma * 1.99 * sqrt(N/2) * sin(vl);
points(x = x1, y = y1, type = "l", lwd = 3, col = "orange");

# Need to figure out where that variation is going.



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
