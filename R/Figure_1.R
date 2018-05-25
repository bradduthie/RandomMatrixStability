#' Find a stabilised system
#' 
#' This function finds a random matrix (M) for which the complex system it
#' represents is unstable before the addition of variation in component
#' response rate, but stable after its addition.
#'
#'@return If successful, a list that includes a matrix (A0) that is unstable in 
#'the absence of component response rate variation and a matrix (A1) that is
#'is identical to A0 but includes variation in component response rate and is
#'stable.
#'@param S Size of the complex system
#'@param C Connectedness of complex system
#'@param Osd Standard deviation of interaction strength
#'@param iters Number of random matrices to simulate
#'@examples
#'find_bgamma(S = 200, C = 0.05, Osd = 0.4, iters = 2);
#'@export
find_bgamma <- function(S = 200, C = 0.05, Osd = 0.4, iters = 10000){
    while(iters > 0){
        A_dat  <- rnorm(n = S * S, mean = 0, sd = Osd);
        A_mat  <- matrix(data = A_dat, nrow = S);
        C_dat  <- rbinom(n = S * S, size = 1, prob = C);
        C_mat  <- matrix(data = C_dat, nrow = S, ncol = S);
        A_mat  <- A_mat * C_mat;
        gammas <- c(rep(1.95, S/2), rep(0.05, S/2))
        mu_gam <- mean(gammas);
        diag(A_mat) <- -1;
        A1     <- gammas * A_mat;
        A0     <- mu_gam * A_mat;
        A0_e   <- eigen(A0)$values;
        A0_r   <- Re(A0_e);
        A0_i   <- Im(A0_e);
        A1_e   <- eigen(A1)$values;
        A1_r   <- Re(A1_e);
        A1_i   <- Im(A1_e);
        if(max(A0_r) >= 0 & max(A1_r) < 0){
            return(list(A0 = A0, A1 = A1));
            break;
        }
        print(iters);
        iters <- iters - 1;
    }
}

#' Find stable bimodal systems
#' 
#' Produces a two panel plot for a matrix (M) before and after the addition of 
#' varying component response rate (gamma).
#'
#'@return A plot as in Figure 1 of the manuscript, showing distribution of 
#'eigenvalues before and after the addition of variation in component response
#'rate
#'@param A0 Matrix (M) before addition of varying component response rate
#'@param A1 Matrix (M) after additoin of varying component response rate
#'@examples
#'load(bi_pr_st);
#'plot_Fig_1(A0 = A0, A1 = A1);
#'@export
plot_Fig_1 <- function(A0, A1){
    S_val       <- dim(A0)[1];
    A0_e        <- eigen(A0)$values;
    A0_r        <- Re(A0_e);
    A0_i        <- Im(A0_e);
    A1_e        <- eigen(A1)$values;
    A1_r        <- Re(A1_e);
    A1_i        <- Im(A1_e);
    A0_vm       <- A0;
    diag(A0_vm) <- NA;
    A0vec       <- as.vector(t(A0_vm));
    A0vec       <- A0vec[is.na(A0vec) == FALSE];
    A1_vm       <- A1;
    diag(A1_vm) <- NA;
    A1vec       <- as.vector(t(A1_vm));
    A1vec       <- A1vec[is.na(A1vec) == FALSE];
    fhalf       <- 1:(0.5*length(A1vec));
    shalf       <- (0.5*length(A1vec)+1):length(A1vec);
    par(mfrow = c(1, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 0, 0));
    plot(A0_r, A0_i, xlim = c(-3.7, 0.3), ylim = c(-2, 2), pch = 4, cex = 0.7,
         xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, asp = 1);
    vl   <- seq(from = 0, to = 2*pi, by = 0.001);
    A0x0 <- sqrt(S_val) * sd(A0vec) * cos(vl) + mean(diag(A0));
    A0y0 <- sqrt(S_val) * sd(A0vec) * sin(vl);
    text(x = -3.5, y = 2.25, labels = "a", cex = 2);
    points(x = A0x0, y = A0y0, type = "l", lwd = 3, col = "grey");
    points(A0_r, A0_i, pch = 4, cex = 0.7);
    plot(A1_r, A1_i, xlim = c(-3.7, 0.3), ylim = c(-2, 2), pch = 4, cex = 0.7,
         xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, asp = 1, 
         col = "black", yaxt = "n");
    vl <- seq(from = 0, to = 2*pi, by = 0.001);
    A0x1a <- sqrt(0.5*S_val) * sd(A1vec[fhalf]) * cos(vl) + mean(diag(A1)[1:(0.5*S_val)]);
    A0y1a <- sqrt(S_val) * sd(A1vec[fhalf]) * sin(vl);
    points(x = A0x1a, y = A0y1a, type = "l", lwd = 3, col = "grey");
    A0x1b <- sqrt(0.5*S_val) * sd(A1vec[shalf]) * cos(vl) + 
        mean( diag(A1)[( (0.5*S_val) + 1 ):S_val] );
    A0y1b <- sqrt(0.5*S_val) * sd(A1vec[shalf]) * sin(vl);
    points(x = A0x1b, y = A0y1b, type = "l", lwd = 3, col = "grey");
    points(A1_r[1:S_val], A1_i[1:S_val],pch = 4, cex = 0.7);   
    text(x = -3.5, y = 2.25, labels = "b", cex = 2);
    mtext(side = 1, "Real", outer = TRUE, line = 3, cex = 2);
    mtext(side = 2, "Imaginary", outer = TRUE, line = 2.5, cex = 2);
}

#' Find stable bimodal systems
#' 
#' Produces random matrices and determines whether or not they are stable using 
#' eigenanalysis given two different component response rates of 1.95 and 0.05 
#' over some number of iterations
#'
#'@return A table showing whether matrices are unstable (0) or stable (1) before
#'(column 1) and after (column 2) including variation in component response rate
#'in eigenanalysis for some number (rows) of random matrices.
#'@param S Size of the complex system
#'@param C Connectedness of complex system
#'@param Osd Standard deviation of interaction strength
#'@param iters Number of random matrices to simulate
#'@examples
#'dat <- stab_bgamma(iters = 4);
#'@export
stab_bgamma <- function(S = 200, C = 0.05, Osd = 0.4, iters = 10000){
    ress     <- matrix(data = 0, nrow = iters, ncol = 2);
    A0_count <- 0;
    A1_count <- 0;
    while(iters > 0){
        A_dat  <- rnorm(n = S * S, mean = 0, sd = Osd);
        A_mat  <- matrix(data = A_dat, nrow = S);
        C_dat  <- rbinom(n = S * S, size = 1, prob = C);
        C_mat  <- matrix(data = C_dat, nrow = S, ncol = S);
        A_mat  <- A_mat * C_mat;
        gammas <- c(rep(1.95, S/2), rep(0.05, S/2))
        mu_gam <- mean(gammas);
        diag(A_mat) <- -1;
        A1     <- gammas * A_mat;
        A0     <- mu_gam * A_mat;
        A0_e   <- eigen(A0)$values;
        A0_r   <- Re(A0_e);
        A0_i   <- Im(A0_e);
        A1_e   <- eigen(A1)$values;
        A1_r   <- Re(A1_e);
        A1_i   <- Im(A1_e);
        if(max(A0_r) < 0){
            ress[iters, 1] <- 1;
            A0_count       <- A0_count + 1;
        }
        if(max(A1_r) < 0){
            ress[iters, 2] <- 1;
            A1_count       <- A1_count + 1;
        }
        print(c(iters, A0_count, A1_count));
        iters <- iters - 1;
    }
    return(ress);
}