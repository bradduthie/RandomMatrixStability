#' Plot eigenvalues of a large system
#' 
#' Produces a two panel plot for a matrix (M) before and after the addition of 
#' uniform variation in component response rate (gamma).
#'
#'@return A plot as in Figure 2  of the manuscript, showing distribution of 
#'eigenvalues before and after the addition of uniform variation in component 
#'response rate
#'@examples
#'plot_Fig_2();
#'@export
plot_Fig_2 <- function(){
    A_comp <- NULL;
    A_dat  <- rnorm(n = 1000000, mean = 0, sd = 0.4);
    A_mat  <- matrix(data = A_dat, nrow = 1000);
    C_dat  <- rbinom(n = 1000 * 1000, size = 1, prob = 1);
    C_mat  <- matrix(data = C_dat, nrow = 1000, ncol = 1000);
    A_mat     <- A_mat * C_mat;
    gammas <- runif(n = 1000, min = 0, max = 2);
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
    A0_vm       <- A0;
    diag(A0_vm) <- NA;
    A0vec       <- as.vector(A0_vm);
    A0vec       <- A0vec[is.na(A0vec) == FALSE];
    A1_vm       <- A1;
    diag(A1_vm) <- NA;
    A1vec       <- as.vector(A1_vm);
    A1vec       <- A1vec[is.na(A1vec) == FALSE];
    par(mfrow = c(1, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 0, 0));
    plot(A0_r, A0_i, xlim = c(-16.5, 15.5), ylim = c(-16.5,15.5), pch = 4, 
         cex = 0.7, xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
         asp = 1, col = "dodgerblue4");
    vl <- seq(from = 0, to = 2*pi, by = 0.001);
    x0 <- sqrt(1000) * sd(A0vec) * cos(vl) + mean(diag(A0));
    y0 <- sqrt(1000) * sd(A0vec) * sin(vl);
    x1 <- sqrt(1000) * sd(A1vec) * cos(vl) + mean(diag(A1));
    y1 <- sqrt(1000) * sd(A1vec) * sin(vl);
    text(x = -15.5, y = 19, labels = "a", cex = 2);
    points(x = x0, y = y0, type = "l", lwd = 3, col = "dodgerblue4");
    points(x = x1, y = y1, type = "l", col = "red", lwd = 3, lty = "dashed");
    plot(A1_r, A1_i, xlim = c(-16.5, 15.5), ylim = c(-16.5,15.5), pch = 4, 
         cex = 0.7, xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, 
         asp = 1, col = "firebrick", yaxt = "n");
    text(x = -15.5, y = 19, labels = "b", cex = 2);
    points(x = x1, y = y1, type = "l", col = "firebrick", lwd = 3);
    points(x = x0, y = y0, type = "l", lwd = 3, lty = "dashed");
    mtext(side = 1, "Real", outer = TRUE, line = 3, cex = 2);
    mtext(side = 2, "Imaginary", outer = TRUE, line = 2.5, cex = 2);
}