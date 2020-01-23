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

sim10 <- rand_rho_var(S = 10, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);

sim15 <- rand_rho_var(S = 15, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);

sim20 <- rand_rho_var(S = 20, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);

sim25 <- rand_rho_var(S = 25, rhos = seq(from = -0.5, to = 0.5, by = 0.05), 
                      iters = 10000, sigma = 0.2);

sim30 <- rand_rho_var(S = 30, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);

sim35 <- rand_rho_var(S = 35, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);

sim40 <- rand_rho_var(S = 40, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);

sim45 <- rand_rho_var(S = 45, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);

write.csv(sim10, "../notebook/sim_results/rhos/sigma0pt2/sim10.csv");
write.csv(sim15, "../notebook/sim_results/rhos/sigma0pt2/sim15.csv");
write.csv(sim20, "../notebook/sim_results/rhos/sigma0pt2/sim20.csv");
write.csv(sim25, "../notebook/sim_results/rhos/sigma0pt2/sim25.csv");


write.csv(sim30, "notebook/sim_results/rhos/sigma0pt2/sim30.csv");
write.csv(sim35, "notebook/sim_results/rhos/sigma0pt2/sim35.csv");
write.csv(sim40, "notebook/sim_results/rhos/sigma0pt2/sim40.csv");
write.csv(sim45, "notebook/sim_results/rhos/sigma0pt2/sim45.csv");



sim10 <- read.csv(file = "sim_results/rhos/sim10.csv");
sim15 <- read.csv(file = "sim_results/rhos/sim15.csv");
sim20 <- read.csv(file = "sim_results/rhos/sim20.csv");
sim30 <- read.csv(file = "sim_results/rhos/sim30.csv");
sim35 <- read.csv(file = "sim_results/rhos/sim35.csv");
par(mar = c(0.25, 0.25, 0.25, 0.25), mfrow = c(3, 2), oma = c(6, 6, 1, 1));
plot(x = sim10[,1], y = sim10[,20], type = "l", pch = 20, lwd = 2, xaxt = "n",
     ylim = c(-1, 3.5), cex.axis = 1.5);
points(x = sim10[,15], y = sim10[,21], type = "l", lwd = 2, lty = "dashed");
abline(h = 0, col = "red", lty = "dotted", lwd = 0.8);
text(x = -0.6, y = 3.2, labels = "S = 10", cex = 2.5);
plot(x = sim15[,1], y = sim15[,20], type = "l", pch = 20, lwd = 2, xaxt = "n",
     ylim = c(-1, 3.5), yaxt = "n", cex.axis = 1.5);
points(x = sim15[,15], y = sim15[,21], type = "l", lwd = 2, lty = "dashed");
abline(h = 0, col = "red", lty = "dotted", lwd = 0.8);
text(x = -0.6, y = 3.2, labels = "S = 15", cex = 2.5);
plot(x = sim20[,1], y = sim20[,20], type = "l", pch = 20, lwd = 2, xaxt = "n",
     ylim = c(-1, 3.5), cex.axis = 1.5);
points(x = sim20[,15], y = sim20[,21], type = "l", lwd = 2, lty = "dashed");
abline(h = 0, col = "red", lty = "dotted", lwd = 0.8);
text(x = -0.6, y = 3.2, labels = "S = 20", cex = 2.5);
plot(x = sim25[,1], y = sim25[,20], type = "l", pch = 20, lwd = 2, xaxt = "n",
     ylim = c(-1, 3.5), yaxt = "n", cex.axis = 1.25);
points(x = sim25[,15], y = sim25[,21], type = "l", lwd = 2, lty = "dashed");
abline(h = 0, col = "red", lty = "dotted", lwd = 0.8);
text(x = -0.6, y = 3.2, labels = "S = 25", cex = 2.5);
plot(x = sim30[,1], y = sim30[,20], type = "l", pch = 20, lwd = 2, 
     ylim = c(-1, 3.5), cex.axis = 1.5);
points(x = sim30[,15], y = sim30[,21], type = "l", lwd = 2, lty = "dashed");
abline(h = 0, col = "red", lty = "dotted", lwd = 0.8);
text(x = -0.6, y = 3.2, labels = "S = 30", cex = 2.5);
plot(x = sim35[,1], y = sim35[,20], type = "l", pch = 20, lwd = 2, 
     ylim = c(-1, 3.5), yaxt = "n", cex.axis = 1.5);
points(x = sim35[,15], y = sim35[,21], type = "l", lwd = 2, lty = "dashed");
abline(h = 0, col = "red", lty = "dotted", lwd = 0.8);
text(x = -0.6, y = 3.2, labels = "S = 35", cex = 2.5);
mtext(text = "E correlation between A_ij and A_ji", side = 1,
      line = 4, outer = TRUE, cex = 1.5);
mtext(text = "E leading real part eigenvalue", side = 2,
      line = 3.5, outer = TRUE, cex = 1.5);











# Completely linear
par(mar = c(1, 1, 1, 1), mfrow = c(3, 2), oma = c(5, 5, 1, 1));

plot(x = sim10[,1], y = sim10[,20], type = "b", pch = 20, lwd = 2, xaxt = "n",
     ylim = c(-1.5, 3));
points(x = sim10[,15], y = sim10[,21], type = "b", lwd = 2, col = "red");
arrows(x0 = sim10[,1], x1 = sim10[,1], y0 = sim10[,20], y1 = sim10[,22],
       type = 3, angle = 90, length = 0.04)
arrows(x0 = sim10[,1], x1 = sim10[,1], y0 = sim10[,20], y1 = sim10[,23],
       type = 3, angle = 90, length = 0.04)
arrows(x0 = sim10[,15], x1 = sim10[,15], y0 = sim10[,21], y1 = sim10[,24],
       type = 3, angle = 90, length = 0.04, col = "red")
arrows(x0 = sim10[,15], x1 = sim10[,15], y0 = sim10[,21], y1 = sim10[,25],
       type = 3, angle = 90, length = 0.04, col = "red")

plot(x = sim15[,1], y = sim15[,20], type = "b", pch = 20, lwd = 2, xaxt = "n",
     yaxt = "n", ylim = c(-0.5, 0.5));
points(x = sim15[,15], y = sim15[,21], type = "b", lwd = 2);

plot(x = sim20[,1], y = sim20[,20], type = "b", pch = 20, lwd = 2, xaxt = "n",
     ylim = c(-0.5, 0.5));
points(x = sim20[,15], y = sim20[,21], type = "b", lwd = 2);

plot(x = sim25[,1], y = sim25[,20], type = "b", pch = 20, lwd = 2, xaxt = "n",
     yaxt = "n", ylim = c(-0.5, 0.5));
points(x = sim25[,15], y = sim25[,21], type = "b", lwd = 2);

plot(x = sim30[,1], y = sim30[,20], type = "b", pch = 20, lwd = 2, 
     ylim = c(-0.5, 0.5));
points(x = sim30[,15], y = sim30[,21], type = "b", lwd = 2);

plot(x = sim35[,1], y = sim35[,20], type = "b", pch = 20, lwd = 2, yaxt = "n",
     ylim = c(-0.5, 0.5));
points(x = sim35[,15], y = sim35[,21], type = "b", lwd = 2);

mtext(text = "E correlation between A_ij and A_ji", side = 1,
      line = 2, outer = TRUE, cex = 1.5);
mtext(text = "E leading real part eigenvalue", side = 2,
      line = 2, outer = TRUE, cex = 1.5);



plot(x = sim30[,1], y = sim30[,20], type = "b", pch = 20, lwd = 2, xaxt = "n",
     ylim = c(-1.5, 3));
points(x = sim30[,15], y = sim30[,21], type = "b", lwd = 2, col = "red");
arrows(x0 = sim30[,1], x1 = sim30[,1], y0 = sim30[,20], y1 = sim30[,22],
       type = 3, angle = 90, length = 0.04)
arrows(x0 = sim30[,1], x1 = sim30[,1], y0 = sim30[,20], y1 = sim30[,23],
       type = 3, angle = 90, length = 0.04)
arrows(x0 = sim30[,15], x1 = sim30[,15], y0 = sim30[,21], y1 = sim30[,24],
       type = 3, angle = 90, length = 0.04, col = "red")
arrows(x0 = sim30[,15], x1 = sim30[,15], y0 = sim30[,21], y1 = sim30[,25],
       type = 3, angle = 90, length = 0.04, col = "red")
abline(h = 0, col = "blue", lty = "dotted")

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










rand_fishing <- function(S, iters, count_stop, int_type = 0, C = 1, 
                         sigma = 0.4){
    tot_res <- NULL;
    iter    <- iters;
    count   <- 0;
    while(iter > 0){
        A0_dat   <- rnorm(n = S * S, mean = 0, sd = sigma);
        A0       <- matrix(data = A0_dat, nrow = S, ncol = S);
        gam1     <- runif(n = S, min = 0, max = 2);
        diag(A0) <- -1;
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
            count    <- count + 1;
            print(count);
        }
        iter <- iter - 1;
        if(count >= count_stop){
            break;
        }
    }
    if(is.null(tot_res) == FALSE){
        colnames(tot_res) <- c("A0_eig", "A1_eig", "A0_rho", "A1_rho", 
                               "A0_comp", "A1_comp");
    }
    return(tot_res);
}



fish <- rand_fishing(S = 30, iters = 10000000, count_stop = 5)




################################################################################
################################################################################
################################################################################
################################################################################

A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = 0);
gam1     <- runif(n = 1000, min = 0, max = 2);
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);
A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 0.2, 0.2));
plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     ylim = c(-24, 24), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
     asp = 1);
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7,  col = "dodgerblue4");


A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = -0.5);
gam1     <- runif(n = 1000, min = 0, max = 2);
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);
A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     ylim = c(-24, 24), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
     asp = 1);
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7, col = "dodgerblue4");

mtext(side = 1, "Real", outer = TRUE, line = 3, cex = 2);
mtext(side = 2, "Imaginary", outer = TRUE, line = 2.5, cex = 2);


################################################################################
################################################################################
################################################################################
################################################################################

A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = 0);
gam1     <- runif(n = 1000, min = 0, max = 2);
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);
A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

par(mfrow = c(2, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 0.2, 0.2));
plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     ylim = c(-24, 24), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
     asp = 1, xaxt = "n");
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7,  col = "dodgerblue4");


A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = -0.2);
gam1     <- runif(n = 1000, min = 0, max = 2);
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);
A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     ylim = c(-24, 24), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
     asp = 1, xaxt = "n", yaxt = "n");
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7, col = "dodgerblue4");


A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = -0.4);
gam1     <- runif(n = 1000, min = 0, max = 2);
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);
A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     ylim = c(-24, 24), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
     asp = 1);
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7, col = "dodgerblue4");


A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = -0.6);
gam1     <- runif(n = 1000, min = 0, max = 2);
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);
A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     ylim = c(-24, 24), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
     asp = 1, yaxt = "n");
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7, col = "dodgerblue4");


mtext(side = 1, "Real", outer = TRUE, line = 3, cex = 2);
mtext(side = 2, "Imaginary", outer = TRUE, line = 2.5, cex = 2);


################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################


A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = 0);
gam1     <- c(rep(0.05, 500), rep(1.95, 500));
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);
A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

par(mfrow = c(1, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 0.2, 0.2));
plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     ylim = c(-28, 28), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
     asp = 1);
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7,  col = "dodgerblue4");




A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = -0.5);
gam1     <- runif(n = 1000, min = 0, max = 2); # c(rep(0.05, 500), rep(1.95, 500));
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);

for(i in 1:dim(A0)[1]){
    for(j in 1:dim(A0)[2]){
        if(i < j){
            A0[i, j] <- A0[j, i];
            A1[i, j] <- A1[j, i];
        }
    }
}

A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     xlim = c(-20, 20), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5)#, 
     #asp = 1);
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7, col = "dodgerblue4");

mtext(side = 1, "Real", outer = TRUE, line = 3, cex = 2);
mtext(side = 2, "Imaginary", outer = TRUE, line = 2.5, cex = 2);


normpdf <- function(mn, sigma){
    x <- seq(from = mn - sigma*4, to = mn + sigma*4, by = 0.0001);
    y <- (1/(sigma * (sqrt(2 * pi)))) * exp(-0.5 * ((x - mn)/sigma)^2);
    return(list(x = x, y = y));
}




sim <- rand_rho_var(S = 1000, rhos = seq(from = -0.95, to = 0.95, by = 0.05), 
                    iters = 100);






A0       <- build_rho_mat(S = 1000, sigma = 0.4, rho = -0.75);
gam1     <- runif(n = 1000, min = 0, max = 2);
A1       <- A0 * gam1;
A0       <- A0 * mean(gam1);
A0_r     <- Re(eigen(A0)$values);
A0_i     <- Im(eigen(A0)$values);
A1_r     <- Re(eigen(A1)$values);
A1_i     <- Im(eigen(A1)$values);

par(mfrow = c(1, 1), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 0.2, 0.2));
plot(x = A1_r, y = A1_i, pch = 4, cex = 0.7, col = "firebrick", 
     ylim = c(-28, 28), xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
     asp = 1);
points(x = A0_r, y = A0_i, pch = 4, cex = 0.7,  col = "dodgerblue4");











sim1000 <- read.csv("../sim1000.csv");
plot(x = sim1000[,1], y = sim1000[,20], type = "l", pch = 20, lwd = 2, 
     xaxt = "n", ylim = c(-1, 3.5), cex.axis = 1.5);
points(x = sim1000[,15], y = sim1000[,21], type = "l", lwd = 2, lty = "dashed");




simt <- rand_rho_var(S = 25, iters = 1000, rhos = seq(from = -0.95, to = 0.95, by = 0.05));


plot(x = simt[,14], y = simt[,20], type = "b", pch = 20, lwd = 2, 
     ylim = c(-1, 3.5), cex.axis = 1.5);
points(x = simt[,15], y = simt[,21], type = "b", lwd = 2, lty = "dashed");


################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################

fxeig <- rand_rho_var(S = 16, iters = 1000, int_type = 0, sigma = (1/4), C = 1,
                      rhos = seq(from = -0.95, to = 0.95, by = 0.05));
                         

plot(x = fxeig[,1], y = fxeig[,20], type = "l", pch = 20, lwd = 2, 
     cex.axis = 1.5, asp = 1);
points(x = fxeig[,15], y = fxeig[,21], type = "l", lwd = 2, lty = "dashed");




fxeig2 <- rand_rho_var(S = 25, iters = 1000, int_type = 0, sigma = (1/5), C = 1,
                      rhos = seq(from = -0.95, to = 0.95, by = 0.05));


plot(x = fxeig2[,1], y = fxeig2[,20], type = "l", pch = 20, lwd = 2, 
     cex.axis = 1.5, asp = 1);
points(x = fxeig2[,15], y = fxeig2[,21], type = "l", lwd = 2, lty = "dashed");







fxeig3 <- rand_rho_var(S = 36, iters = 1000, int_type = 0, sigma = (1/6), C = 1,
                       rhos = seq(from = -0.95, to = 0.95, by = 0.05));


plot(x = fxeig3[,1], y = fxeig3[,20], type = "l", pch = 20, lwd = 2, 
     cex.axis = 1.5, asp = 1);
points(x = fxeig3[,15], y = fxeig3[,21], type = "l", lwd = 2, lty = "dashed");




fxeig4 <- rand_rho_var(S = 100, iters = 100, int_type = 0, sigma = (1/10), C = 1,
                       rhos = seq(from = -0.95, to = 0.95, by = 0.05));

plot(x = fxeig4[,1], y = fxeig4[,20], type = "l", pch = 20, lwd = 2, 
     cex.axis = 1.5, asp = 1);
points(x = fxeig4[,15], y = fxeig4[,21], type = "l", lwd = 2, lty = "dashed");




fxeig5 <- rand_rho_var(S = 100, iters = 10000, int_type = 0, sigma = (1/10), C = 1,
                       rhos = seq(from = -0.95, to = 0.95, by = 0.05));

plot(x = fxeig5[,1], y = fxeig5[,20], type = "l", pch = 20, lwd = 2, 
     cex.axis = 1.5, asp = 1);
points(x = fxeig5[,15], y = fxeig5[,21], type = "l", lwd = 2, lty = "dashed");


# The only things that are changing, really, along the line is rho. Perhaps this
# has some sort of nonlinear response to gamma? Positive varation and negative
# rho causes the switch?



fxeig6 <- rand_rho_var(S = 1600, iters = 10, int_type = 0, sigma = (1/40), C = 1,
                       rhos = seq(from = -0.95, to = 0.95, by = 0.05));


fxeig6 <- read.csv("../fxeig6.csv");

plot(x = fxeig6[,1], y = fxeig6[,20], type = "l", pch = 20, lwd = 2, 
     cex.axis = 1.5, asp = 1);
points(x = fxeig6[,15], y = fxeig6[,21], type = "l", lwd = 2, lty = "dashed");

points(x = 0, y = 0, cex = 2, col = "red");


points(x = fxeig6[,15], y = fxeig6[,24], type = "l", col = "blue", 
       lty = "dashed");
points(x = fxeig6[,15], y = fxeig6[,25], type = "l", col = "blue", 
       lty = "dashed");

################################################################################
################################################################################
################################################################################
################################################################################

mat_cov <- function(mat){
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
    rho <- cov(ut, lt);
    return(rho);
}

mat_cor <- function(mat){
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

mat_var <- function(mat){
    diag(mat) <- NA;
    moffs     <- as.vector(mat);
    moffs     <- moffs[!is.na(moffs)];
    mvar      <- var(moffs);
    return(mvar);
}



S      <- 1000;
v_odA0 <- NULL;
v_odA1 <- NULL;
cvodA0 <- NULL;
cvodA1 <- NULL;
crodA0 <- NULL;
crodA1 <- NULL;
v_gamm <- NULL;
for(i in 1:1000){
    A0_dat      <- rnorm(n = S * S, mean = 0, sd = 1/S);
    A0          <- matrix(data = A0_dat, nrow = S, ncol = S);
    gam_dat     <- runif(n = S, min = 0, max = 2);
    gamma       <- matrix(data = 0, nrow = S, ncol = S);
    diag(gamma) <- gam_dat;
    diag(A0)    <- -1 
    A1          <- gamma %*% A0;
    A0          <- A0 * mean(gamma);
    v_odA0[i]   <- mat_var(A0);
    v_odA1[i]   <- mat_var(A1);
    cvodA0[i]   <- mat_cov(A0);
    cvodA1[i]   <- mat_cov(A1);
    crodA0[i]   <- mat_cor(A0);
    crodA1[i]   <- mat_cor(A1);
    v_gamm[i]   <- var(gam_dat);
    print(i);
}

# It appears that E[ln(var(A0)) / ln(var(A1))] \approx 2


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

# Relate to zero law of biology? Variance in the correlation will inherently 
# increase, meaning that there is a higher probability that a negative
# correlation will occur and lead to stability. There is a bound at 1, but 
# multiplying by a vector makes more situations with a -1, assuming the
# initial correlation of M_{i,j} and M_{j,i} is not uniform?


xx  <- dat$A0_rho;
yy  <- (1 + 1/3)*(dat$A1_rho);
dif <- dat$A0_rho - (1 + 1/3)*(dat$A1_rho);

plot(x = xx, y = yy, asp = 1, xlim = c(-1, 1), ylim = c(-1, 1));
cor.test(xx, yy);









get_varchange <- function(S, bs){
  sigma   <- 1/(1*sqrt(S));
  a       <- 0;
  res_mat <- NULL;
  for(b in bs){
      A_dat  <- rnorm(n = S * S, mean = 0, sd = sigma);
      A_mat  <- matrix(data = A_dat, nrow = S);
      C_dat  <- rbinom(n = S * S, size = 1, prob = C);
      C_mat  <- matrix(data = C_dat, nrow = S, ncol = S);
      A_mat  <- A_mat * C_mat;
      gammas <- rpois(n = S, lambda = b); # runif(n = S, min = a, max = b);
      vargam <- b; # (1/12) * (b-a)^2;
      diag(A_mat) <- 0;
      A1     <- gammas * A_mat;
      A0     <- A_mat;
      A0_e   <- eigen(A0)$values;
      A0_r   <- Re(A0_e);
      A0_i   <- Im(A0_e);
      A1_e   <- eigen(A1)$values;
      A1_r   <- Re(A1_e);
      A1_i   <- Im(A1_e);
      A_test <- A0;
      diag(A_test) <- NA;
      A0_var <-var(as.vector(A_test), na.rm = TRUE);
      A_test <- A1;
      diag(A_test) <- NA;
      A1_var <-var(as.vector(A_test), na.rm = TRUE);
      hackradii <- (range(A1_r)[2] - range(A1_r)[1]) / 2
      keepit <- c(sigma, sigma^2, b, vargam, A0_var, A1_var, hackradii);
      res_mat <- rbind(res_mat, keepit);
      print(b);
  }
  colnames(res_mat) <- c("sigma", "sigma^2", "b", "vargam", "A0_var", "A1_var",
                         "hackradii");
  return(res_mat);
}


dat <- get_varchange(S = 1024, bs = seq(from = 0.1, to = 48, by = 1));


plot(x = dat[,4], dat[,7], xlab = "Variation of gamma", 
     ylab = "Radius of eigenvalue circle");
var_form <- 2 * SCsigma * sqrt(dat[,4]);
points(x = dat[,4], y = var_form, type = "l");






tt <- 144 - exp(dat[,4])
plot(x = dat[,4], tt, xlab = "Variation of gamma", 
     ylab = "Radius of eigenvalue circle");













################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
S     <- 2500;
C     <- 1;
dval  <- 0;
sigma <- 1/(1*sqrt(S));
a     <- 0;
b     <- 2;
mn    <- 0;

A_dat  <- rnorm(n = S * S, mean = mn, sd = sigma);
A_mat  <- matrix(data = A_dat, nrow = S);
C_dat  <- rbinom(n = S * S, size = 1, prob = C);
C_mat  <- matrix(data = C_dat, nrow = S, ncol = S);
A_mat  <- A_mat * C_mat;
gammas <- runif(n = S, min = a, max = b);
vargam <- (1/12) * (b-a)^2;

diag(A_mat) <- dval;
A1     <- gammas * A_mat;
A0     <- A_mat;
A0     <- A0 * mean(gammas);
A0_e   <- eigen(A0)$values;
A0_r   <- Re(A0_e);
A0_i   <- Im(A0_e);
A1_e   <- eigen(A1)$values;
A1_r   <- Re(A1_e);
A1_i   <- Im(A1_e);

plot(A1_r, A1_i, pch = 4, cex = 0.7, xlab = "", ylab = "", cex.lab = 1.3, 
     cex.axis = 1.5, asp = 1, col = "firebrick", yaxt = "n");

joint_var <- ((sigma^2) * vargam) + ((sigma^2) * mean(gammas)^2) + (vargam * (dval*mn)^2);
crad      <- sqrt(joint_var) * sqrt(S);
vl        <- seq(from = 0, to = 2*pi, by = 0.001);
A1x1a     <- crad * cos(vl) + mean(diag(A1));
A1y1a     <- crad * sin(vl);
points(x = A1x1a, y = A1y1a, type = "l", lwd = 3, col = "dodgerblue4");

points(A0_r, A0_i, pch = 4, cex = 0.7, xlab = "", ylab = "", cex.lab = 1.3, 
     cex.axis = 1.5, asp = 1, col = "dodgerblue4", yaxt = "n");
A0x1a     <- mean(gammas) * sigma * sqrt(S) * cos(vl) + mean(diag(A0));
A0y1a     <- mean(gammas) * sigma * sqrt(S) * sin(vl);
points(x = A0x1a, y = A0y1a, type = "l", lwd = 3, col = "dodgerblue4");

calc_var <- A1;
diag(calc_var) <- NA;
calc_var  <- as.numeric(calc_var);
cvar <- var(calc_var[!is.na(calc_var)]);


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# Note that sigma * sqrt(SC) = 1 in the case of the below.
rhos <- read.csv("../notebook/sim_results/rhos/sigma0pt2/sim25.csv");

v_gammas <-  (1/12) * (2-0)^2;
covars   <- rhos$A0_rho * sigma * sigma;

t1 <- rhos$A0_rho;
t2 <- rhos$A1_rho;

joint_var <- ((sigma^2) * v_gammas) + ((sigma^2) * 1) + (v_gammas * 0);
joint_sd  <- sqrt(joint_var);


plot(x = t2, y = t1, asp = 1)

var(v1)*var(v2) + var(v1)*(mean(v2)^2) + var(v2)*(mean(v1)^2);


#
#
# Try to figure out how the variances for the two gamma case add up. 
#  - Need to look at many examples to see if the variance calculation is correct
#  - Maybe try to use 3 gammas, just to triple check this all makes sense.
#
# Surely this is proportional?
#



A0_eigs128 <- NULL;
A1_eigs128 <- NULL;
for(i in 1:100000){
    Snew  <- 128;
    sigma <- 2 / sqrt(256);
    A0_new_dat <- rnorm(n = Snew, mean = 0, sd = sigma);
    A0_new_mat <- matrix(data = A0_new_dat, nrow = Snew, ncol = Snew);
    gammas_new <- c(rep(0.05, Snew/2), rep(1.95, Snew/2));
    mu_gam_new <- mean(gammas_new);
    diag(A0_new_mat) <- -1 * dval;
    A1_new    <- gammas_new * A0_new_mat;
    A0_new    <- mu_gam_new * A0_new_mat;
    A0_new_e   <- eigen(A0_new)$values;
    A0_new_r   <- Re(A0_new_e);
    A0_new_i   <- Im(A0_new_e);
    A1_new_e   <- eigen(A1_new)$values;
    A1_new_r   <- Re(A1_new_e);
    A1_new_i   <- Im(A1_new_e);
    A0_eigs128[i] <-  max(A0_new_r);
    A1_eigs128[i] <-  max(A1_new_r);
    if(i %% 10000 == 0 ){
        print(i);
    }
}

A0_eigs256 <- NULL;
A1_eigs256 <- NULL;
for(i in 1:100000){
    Snew  <- 256;
    sigma <- 2 / sqrt(256);
    A0_new_dat <- rnorm(n = Snew, mean = 0, sd = sigma);
    A0_new_mat <- matrix(data = A0_new_dat, nrow = Snew, ncol = Snew);
    gammas_new <- c(rep(0.05, Snew/2), rep(1.95, Snew/2));
    mu_gam_new <- mean(gammas_new);
    diag(A0_new_mat) <- -1 * dval;
    A1_new    <- gammas_new * A0_new_mat;
    A0_new    <- mu_gam_new * A0_new_mat;
    A0_new_e   <- eigen(A0_new)$values;
    A0_new_r   <- Re(A0_new_e);
    A0_new_i   <- Im(A0_new_e);
    A1_new_e   <- eigen(A1_new)$values;
    A1_new_r   <- Re(A1_new_e);
    A1_new_i   <- Im(A1_new_e);
    A0_eigs256[i] <-  max(A0_new_r);
    A1_eigs256[i] <-  max(A1_new_r);
    if(i %% 10000 == 0 ){
        print(i);
    }
}








