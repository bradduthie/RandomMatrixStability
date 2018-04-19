make_gammas <- function(nn = 10, distribution = 1, mn = 1, sdd = 1){
    if(distribution == 0){
        dat          <- rep(x = mn, times = nn);
    }
    if(distribution == 1){
        mval         <- sdd * sqrt(12);
        dat          <- runif(n = nn, min = 0, max = mval);
    }
    if(distribution == 2){
        dat          <- rexp(n = nn);
        dat          <- sdd * ((dat - mean(dat)) / sd(dat));
        dat          <- dat - min(dat) + min(abs(dat));
    }
    if(distribution == 3){
        dat          <- rbeta(n = nn, shape1 = 0.5, shape2 = 0.5);
        dat          <- sdd * ((dat - mean(dat)) / sd(dat));
        dat          <- dat - min(dat) + min(abs(dat));
    }
    if(distribution == 4){
        dat          <- rgamma(n = nn, shape = 2, scale = 2);
        dat          <- sdd * ((dat - mean(dat)) / sd(dat));
        dat          <- dat - min(dat) + min(abs(dat));
    }
    return(dat);
}

show_gammas <- function(nn = 1000000, sdd = 1, mn = 10){
    y1 <- make_gammas(nn = nn, distribution = 1, sdd = sdd, mn = mn);
    y2 <- make_gammas(nn = nn, distribution = 2, sdd = sdd, mn = mn);
    y3 <- make_gammas(nn = nn, distribution = 3, sdd = sdd, mn = mn);
    y4 <- make_gammas(nn = nn, distribution = 4, sdd = sdd, mn = mn);
    par(mfrow = c(2, 2), oma = c(6, 6, 1, 1), mar = c(4, 0.5, 0.5, 0.5));
    h1 <- hist(y1, breaks = 1000, plot = FALSE);
    hist(y1, breaks = 1000, yaxt = "n", main = "", xlab = "", cex.axis = 1.5,
         ylim = c(0, max(h1$counts)*1.2));
    mtext("a", adj = 0, line = -2.3, cex = 2, 
          at = par("usr")[1]+0.90*diff(par("usr")[1:2]));
    box();
    h2 <- hist(y2, breaks = 1000, plot = FALSE);
    hist(y2, breaks = 1000, yaxt = "n", main = "", xlab = "", cex.axis = 1.5,
         ylim = c(0, max(h2$counts)*1.2));
    mtext("b", adj = 0, line = -2.3, cex = 2, 
          at = par("usr")[1]+0.90*diff(par("usr")[1:2]));
    box();
    h3 <- hist(y3, breaks = 1000, plot = FALSE);
    hist(y3, breaks = 1000, yaxt = "n", main = "", xlab = "", cex.axis = 1.5,
         ylim = c(0, max(h3$counts)*1.2));
    mtext("c", adj = 0, line = -2.3, cex = 2, 
          at = par("usr")[1]+0.90*diff(par("usr")[1:2]));
    box();
    h4 <- hist(y4, breaks = 1000, plot = FALSE);
    hist(y4, breaks = 1000, yaxt = "n", main = "", xlab = "", cex.axis = 1.5,
         ylim = c(0, max(h4$counts)*1.2));
    mtext("d", adj = 0, line = -2.3, cex = 2, 
          at = par("usr")[1]+0.90*diff(par("usr")[1:2]));
    box();
    mtext(expression(paste("Component ",gamma," value")),
          outer=TRUE,side=1,line=3.0,cex=2);
    mtext(expression(paste("Relative frequency")),
          outer=TRUE,side=2,line=2.5,cex=2);
}

eigen_cloud <- function(sp, iters, int_type = 0, gamma_sd = 1){
    co0  <- NULL;
    co1  <- NULL;
    co2  <- NULL;
    co3  <- NULL;
    co4  <- NULL;
    iter <- iters;
    while(iter > 0){
        A0_dat   <- rnorm(n = sp * sp, mean = 0, sd = 0.4);
        A0       <- matrix(data = A0_dat, nrow = sp, ncol = sp);
        A0       <- species_interactions(mat = A0, type = int_type);
        diag(A0) <- -1;
        gam0     <- make_gammas(nn = sp, distribution = 0, sdd = gamma_sd);
        gam1     <- make_gammas(nn = sp, distribution = 1, sdd = gamma_sd);
        gam2     <- make_gammas(nn = sp, distribution = 2, sdd = gamma_sd);
        gam3     <- make_gammas(nn = sp, distribution = 3, sdd = gamma_sd);
        gam4     <- make_gammas(nn = sp, distribution = 4, sdd = gamma_sd);
        A1       <- A0 * gam1;
        A2       <- A0 * gam2;
        A3       <- A0 * gam3;
        A4       <- A0 * gam4;
        A0       <- A0 * gam0;
        Re_A0    <- Re(eigen(A0)$values);
        Im_A0    <- Im(eigen(A0)$values);
        Re_A1    <- Re(eigen(A1)$values);
        Im_A1    <- Im(eigen(A1)$values);
        Re_A2    <- Re(eigen(A2)$values);
        Im_A2    <- Im(eigen(A2)$values);
        Re_A3    <- Re(eigen(A3)$values);
        Im_A3    <- Im(eigen(A3)$values);
        Re_A4    <- Re(eigen(A4)$values);
        Im_A4    <- Im(eigen(A4)$values);
        e_coord0 <- cbind(Re_A0, Im_A0);
        e_coord1 <- cbind(Re_A1, Im_A1);
        e_coord2 <- cbind(Re_A2, Im_A2);
        e_coord3 <- cbind(Re_A3, Im_A3);
        e_coord4 <- cbind(Re_A4, Im_A4);
        co0      <- rbind(co0, e_coord0);
        co1      <- rbind(co1, e_coord1);
        co2      <- rbind(co2, e_coord2);
        co3      <- rbind(co3, e_coord3);
        co4      <- rbind(co4, e_coord4);
        iter     <- iter - 1;
    }
    return(list(d0 = co0, d1 = co1, d2 = co2, d3 = co3, d4 = co4));
} 

eig_plot <- function(ecl){
    xlims <- c(min(ecl$d1[,1]), max(ecl$d1[,1]));
    ylims <- c(min(ecl$d1[,2]), max(ecl$d1[,2]));
    tpink <- "#FFC0CB14"
    par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), oma = c(5, 5, 1, 1));
    # Uniform distribution (a)
    plot(x = ecl$d0[,1], y = ecl$d0[,2], pch = 4, col = "pink", cex = 0.3, 
         xlim = xlims, ylim = ylims, cex.axis = 1.5, xaxt = "n", yaxt = "n");
    points(x = ecl$d1[,1], y = ecl$d1[,2], pch = 4, col = "black", cex = 0.3);
    points(x = ecl$d0[,1], y = ecl$d0[,2], pch = 4, col = tpink, cex = 0.3);
    axis(side = 2, at = c(-4, -2, 0, 2, 4), cex.axis = 1.5);
    abline(h = 0, lwd = 2);
    abline(v = 0, lwd = 2);
    box(lwd = 2);
    text(x = -8, y = 4.5, labels = "a", cex = 3);
    # Exponential distribution (b)
    plot(x = ecl$d0[,1], y = ecl$d0[,2], pch = 4, col = "pink", cex = 0.3, 
         xlim = xlims, ylim = ylims, cex.axis = 1.5, xaxt = "n", yaxt = "n");
    points(x = ecl$d2[,1], y = ecl$d2[,2], pch = 4, col = "black", cex = 0.3);
    points(x = ecl$d0[,1], y = ecl$d0[,2], pch = 4, col = tpink, cex = 0.3);
    abline(h = 0, lwd = 2);
    abline(v = 0, lwd = 2);
    box(lwd = 2);
    text(x = -8, y = 4.5, labels = "b", cex = 3);
    # U-shaped distribution (c)
    plot(x = ecl$d0[,1], y = ecl$d0[,2], pch = 4, col = "pink", cex = 0.3, 
         xlim = xlims, ylim = ylims, cex.axis = 1.5, xaxt = "n", yaxt = "n");
    points(x = ecl$d3[,1], y = ecl$d3[,2], pch = 4, col = "black", cex = 0.3);
    points(x = ecl$d0[,1], y = ecl$d0[,2], pch = 4, col = tpink, cex = 0.3);
    axis(side = 2, at = c(-4, -2, 0, 2, 4), cex.axis = 1.5);
    axis(side = 1, at = c(-12, -9, -6, -3, 0, 3, 6), cex.axis = 1.5);
    abline(h = 0, lwd = 2);
    abline(v = 0, lwd = 2);
    box(lwd = 2);
    text(x = -8, y = 4.5, labels = "c", cex = 3);
    # Poisson-like distribution (d)
    plot(x = ecl$d0[,1], y = ecl$d0[,2], pch = 4, col = "pink", cex = 0.3, 
         xlim = xlims, ylim = ylims, cex.axis = 1.5, xaxt = "n", yaxt = "n");
    points(x = ecl$d4[,1], y = ecl$d4[,2], pch = 4, col = "black", cex = 0.3);
    points(x = ecl$d0[,1], y = ecl$d0[,2], pch = 4, col = tpink, cex = 0.3);
    axis(side = 1, at = c(-12, -9, -6, -3, 0, 3, 6), cex.axis = 1.5);
    abline(h = 0, lwd = 2);
    abline(v = 0, lwd = 2);
    box(lwd = 2);
    text(x = -8, y = 4.5, labels = "d", cex = 3);
    # Outer margin
    mtext(side = 1, text = "Real", cex = 2, outer = TRUE, line = 2);
    mtext(side = 2, text = "Imaginary", cex = 2, outer = TRUE, line = 2);
}
