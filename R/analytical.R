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
x0 <- sigma * 0.01 * sqrt(S/2) * cos(vl) + mean(diag(M3));
y0 <- sigma * 0.01 * sqrt(S/2) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3, col = "orange");
# Check out what happens when you change S
# Something is increasing the size of the circle.
# Solve this analytically.

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
x0 <- sigma * 0.01 * sqrt(S/2) * cos(vl) + mean(diag(M3));
y0 <- sigma * 0.01 * sqrt(S/2) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3, col = "orange");
# Add the new M4 ones
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x1 <- sigma * 1.99 * sqrt(S/2) * cos(vl) + mean(diag(M4));
y1 <- sigma * 1.99 * sqrt(S/2) * sin(vl);
points(x = x1, y = y1, type = "l", lwd = 3, col = "orange");

# Need to figure out where that variation is going.


################################################################################
# See the relationship between rho and the leading eigenvalue real component
################################################################################
rand_rho_var <- function(S, rhos, iters, int_type = 0, rmx = 0.4, C = 1, 
                         by = 1, sigma = 0.4){
    tot_res <- NULL;
    fea_res <- NULL;
    rho_res <- NULL;
    cmplxty <- NULL;
    for(i in 1:length(rhos)){
        iter           <- iters;
        tot_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 8);
        fea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 7);
        rho_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 6);
        cmplxty[[i]]   <- matrix(data = 0, nrow = iter, ncol = 2);
        while(iter > 0){
            r_vec    <- rnorm(n = S, mean = 0, sd = rmx);
            A0       <- build_rho_mat(S = S, sigma = sigma, rho = rhos[i]);
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

simpleboot <- function(freqs,repli=1000,alpha=0.05){    
    vals  <- NULL;                                                          
    i     <- 0;                                         
    while(i < repli){                                                   
        boot  <- sample(x=freqs,size=length(freqs),replace=TRUE);           
        strap <- mean(boot);                                                
        vals  <- c(vals,strap);                                             
        i     <- i + 1;                                                     
    }                                                                       
    vals   <- sort(x=vals,decreasing=FALSE);                            
    lowCI  <- vals[round((alpha*0.5)*repli)];                               
    highCI <- vals[round((1-(alpha*0.5))*repli)];                           
    CIs    <- c(lowCI,highCI);                                          
    return(CIs);                                                            
}  

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

sim10 <- rand_rho_var(S = 10, rhos = seq(from = -0.5, to = 0.5, by = 0.05), 
                      iters = 10000);

sim15 <- rand_rho_var(S = 15, rhos = seq(from = -0.5, to = 0.5, by = 0.05), 
                      iters = 10000);

sim20 <- rand_rho_var(S = 20, rhos = seq(from = -0.5, to = 0.5, by = 0.05), 
                      iters = 10000);

sim25 <- rand_rho_var(S = 25, rhos = seq(from = -0.5, to = 0.5, by = 0.05), 
                      iters = 10000);

sim30 <- rand_rho_var(S = 30, rhos = seq(from = -0.5, to = 0.5, by = 0.05), 
                      iters = 10000);

sim35 <- rand_rho_var(S = 35, rhos = seq(from = -0.5, to = 0.5, by = 0.05), 
                      iters = 10000);

# Completely linear
par(mar = c(1, 1, 1, 1), mfrow = c(3, 2), oma = c(5, 5, 1, 1));

plot(x = sim10[,1], y = sim10[,20], type = "b", pch = 20, lwd = 2, xaxt = "n");
points(x = sim10[,15], y = sim10[,21], type = "b", lwd = 2);

plot(x = sim15[,1], y = sim15[,20], type = "b", pch = 20, lwd = 2, xaxt = "n",
     yaxt = "n");
points(x = sim15[,15], y = sim15[,21], type = "b", lwd = 2);

plot(x = sim20[,1], y = sim20[,20], type = "b", pch = 20, lwd = 2, xaxt = "n");
points(x = sim20[,15], y = sim20[,21], type = "b", lwd = 2);

plot(x = sim25[,1], y = sim25[,20], type = "b", pch = 20, lwd = 2, xaxt = "n",
     yaxt = "n");
points(x = sim25[,15], y = sim25[,21], type = "b", lwd = 2);

plot(x = sim30[,1], y = sim30[,20], type = "b", pch = 20, lwd = 2);
points(x = sim30[,15], y = sim30[,21], type = "b", lwd = 2);

plot(x = sim35[,1], y = sim35[,20], type = "b", pch = 20, lwd = 2, yaxt = "n");
points(x = sim35[,15], y = sim35[,21], type = "b", lwd = 2);

mtext(text = "E correlation between A_ij and A_ji", side = 1,
      line = 0, outer = TRUE, cex = 1.5);
mtext(text = "E leading real part eigenvalue", side = 2,
      line = 1, outer = TRUE, cex = 1.5);






par(mar = c(1, 1, 1, 1), mfrow = c(3, 2), oma = c(5, 5, 1, 1));

plot(x = sim10[,1], y = sim10[,20], type = "b", pch = 20, lwd = 2, xaxt = "n");
points(x = sim10[,15], y = sim10[,21], type = "b", lwd = 2);
abline(h = 0, col = "red");

plot(x = sim15[,1], y = sim15[,20], type = "b", pch = 20, lwd = 2, xaxt = "n");
points(x = sim15[,15], y = sim15[,21], type = "b", lwd = 2);
abline(h = 0, col = "red");

plot(x = sim20[,1], y = sim20[,20], type = "b", pch = 20, lwd = 2, xaxt = "n");
points(x = sim20[,15], y = sim20[,21], type = "b", lwd = 2);
abline(h = 0, col = "red");

plot(x = sim25[,1], y = sim25[,20], type = "b", pch = 20, lwd = 2, xaxt = "n");
points(x = sim25[,15], y = sim25[,21], type = "b", lwd = 2);
abline(h = 0, col = "red");

plot(x = sim30[,1], y = sim30[,20], type = "b", pch = 20, lwd = 2);
points(x = sim30[,15], y = sim30[,21], type = "b", lwd = 2);
abline(h = 0, col = "red");

plot(x = sim35[,1], y = sim35[,20], type = "b", pch = 20, lwd = 2);
points(x = sim35[,15], y = sim35[,21], type = "b", lwd = 2);
abline(h = 0, col = "red");

mtext(text = "E correlation between A_ij and A_ji", side = 1,
      line = 0, outer = TRUE, cex = 1.5);
mtext(text = "E leading real part eigenvalue", side = 2,
      line = 1, outer = TRUE, cex = 1.5);










rand_fishing <- function(S, iters, int_type = 0, rmx = 0.4, C = 1, 
                         by = 1, sigma = 0.4){
    tot_res <- NULL;
    for(i in 1:length(rhos)){
        iter           <- iters;
        while(iter > 0){
            A0_dat   <- rnorm(n = S * S, mean = 0, sd = sigma);
            A0       <- matrix(data = A0_dat, nrow = S, ncol = S);
            gam1     <- runif(n = S, min = 0, max = 2);
            A1       <- A0 * gam1;
            A0       <- A0 * mean(gam1);
            A0_eig   <- max(Re(eigen(A0)$values));
            A1_eig   <- max(Re(eigen(A1)$values));
            if(A0_eig > 0 & A1_eig < 0){
                A0_rho   <- mat_rho(A0);
                A1_rho   <- mat_rho(A1);
                A0_comp  <- get_complexity(A0);
                A1_comp  <- get_complexity(A1);
                it_res   <- c(A0_eig, A1_eig, A0_rho, A1_rho, A0_comp, A1_comp);
                tot_res  <- rbind(tot_res, it_res);
            }
            iter <- iter - 1;
        }
    }
    colnames(tot_res) <- c("A0_eig", "A1_eig", "A0_rho", "A1_rho", "A0_comp",
                           "A1_comp");
    return(tot_res);
}

