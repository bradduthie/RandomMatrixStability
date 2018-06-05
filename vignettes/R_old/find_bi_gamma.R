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


plot_eigs <- function(A0, A1){
    S           <- dim(A0)[1];
    A0_e        <- eigen(A0)$values;
    A0_r        <- Re(A0_e);
    A0_i        <- Im(A0_e);
    A1_e        <- eigen(A1)$values;
    A1_r        <- Re(A1_e);
    A1_i        <- Im(A1_e);
    A0_vm       <- A0; 
    diag(A0_vm) <- NA;
    A0vec       <- as.vector(A0_vm);
    A0vec       <- A0vec[is.na(A0vec) == FALSE];
    A1_vm       <- A1;
    diag(A1_vm) <- NA;
    A1vec       <- as.vector(A1_vm);
    A1vec       <- A1vec[is.na(A1vec) == FALSE];
    xL0         <- min(c(A0_r, A1_r));
    xL1         <- max(c(A0_r, A1_r, 1));
    yL0         <- min(c(A0_i, A1_i));
    yL1         <- max(c(A0_i, A1_i, 1));
    par(mfrow = c(1, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 0, 0));
    plot(A0_r, A0_i, xlim = c(xL0, xL1), ylim = c(yL0, yL1), pch = 4, cex = 0.7,
         xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, asp = 1);
    vl <- seq(from = 0, to = 2*pi, by = 0.001);
    x0 <- sqrt(S) * sd(A0vec) * cos(vl) + mean(diag(A0));
    y0 <- sqrt(S) * sd(A0vec) * sin(vl);
    x1 <- sqrt(S) * sd(A1vec) * cos(vl) + mean(diag(A1));
    y1 <- sqrt(S) * sd(A1vec) * sin(vl);
    text(x = -xL0, y = yL1, labels = "a", cex = 2);
    points(x = x0, y = y0, type = "l", lwd = 3);
    abline(v = 0);
    plot(A1_r, A1_i, xlim = c(xL0, xL1), ylim = c(yL0, yL1), pch = 4, cex = 0.7,
         xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.5, asp = 1, 
         col = "red", yaxt = "n");
    text(x = -xL0, y = yL1, labels = "b", cex = 2);
    points(x = x0, y = y0, type = "l", lwd = 3, lty = "dashed");
    abline(v = 0)
    mtext(side = 1, "Real", outer = TRUE, line = 3, cex = 2);
    mtext(side = 2, "Imaginary", outer = TRUE, line = 2.5, cex = 2);
}




