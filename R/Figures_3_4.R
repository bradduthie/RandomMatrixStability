#' Plot comparison of gamma = 1 versus var(gamma)
#' 
#' Produces a plot comparing random matrix results
#'
#'@return A plot as in Figures 3 of the manuscript, showing the stability of
#'random matrices on a log scale in simulations where component response rate is
#'constant versus simulations in which it varies
#'@param dat Table of simulation results
#'@param S_s Maximum size of S to use
#'@export
plot_stables <- function(dat, S_s = 32){
    Ns          <- 1:S_s;
    par(oma = c(6, 6, 1, 6), mar = c(0.5, 0.5, 0.5, 0.5));
    #=================================
    bar_dat                      <- t(cbind(dat[Ns,3], dat[Ns,5]));
    log_bar_dat                  <- log(bar_dat);
    log_bar_dat[log_bar_dat < 0] <- 0; 
    barplot(log_bar_dat, beside = TRUE, col = c("dodgerblue4", "firebrick"),
            names.arg = dat[Ns,1], ylim = c(0, 16), xlab = "",
            ylab = "Ln number of stable communities", cex.lab = 1, 
            cex.axis = 1.25, xlim = c(1, 94), cex.names = 1, yaxt = "n");
    axis(side = 2, at = c(0, 2, 4, 6, 8, 10, 12, 14), cex.axis = 1.5);
    box(lwd = 2);
    par(new = TRUE);
    y1     <- dat[1:S_s,6] / (dat[1:S_s,5]);
    x1     <- seq(from = 2.132, to = 15.1112, length = S_s);
    plot(x = x1, y = y1, xaxt = "n", yaxt = "n", lwd = 2, ylim = c(0, 1.1),
         xlab = "", ylab = "", type = "b", xlim = c(2, 15), pch = 20, 
         cex = 1, col = "black", yaxs="i");
    points(x = x1, y = y1, lwd = 2, type = "l", col = "black");
    axis(side = 4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = 1.5);
    legend("topleft", c(expression(paste(sigma[gamma]^2," = 0")), 
                        expression(paste(sigma[gamma]^2," = 1/3"))), 
           pch=15, col=c("dodgerblue4","firebrick"), cex = 1.5, horiz = TRUE);
    #=================================
    mtext(side = 1, text = "System size (S)", cex = 2, outer = TRUE, 
          line = 3.0);
    mtext(side = 2, text = "Ln number of stable systems", cex = 2, 
          outer = TRUE, line = 3.5);
    mtext(side = 4, text = expression(
        paste("Pr. of systems stable due to Var(",gamma,")")),
        cex = 2, outer = TRUE, line = 3.5);
}


#' Plot comparison of gamma = 1 versus var(gamma) (genetic algorithm)
#' 
#' Produces a plot comparing random matrix results
#'
#'@return A plot as in Figures 4 of the manuscript, showing the stability of
#'random matrices on a log scale in simulations where component response rate is
#'constant versus simulations in which it varies in a targetted way as a
#'consequence of a genetic algorithm search
#'@param dat Table of simulation results
#'@param S_s Maximum size of S to use
#'@export
plot_stable_4 <- function(dat, S_s = 40){
    Ns          <- 1:S_s;
    par(oma = c(6, 6, 1, 6), mar = c(0.5, 0.5, 0.5, 0.5));
    #=================================
    bar_dat                      <- t(cbind(dat[Ns,3], dat[Ns,5]));
    log_bar_dat                  <- log(bar_dat);
    log_bar_dat[log_bar_dat < 0] <- 0; 
    barplot(log_bar_dat, beside = TRUE, col = c("dodgerblue4", "firebrick"),
            names.arg = dat[Ns,1], ylim = c(0, 16), xlab = "",
            ylab = "Ln number of stable communities", cex.lab = 1, 
            cex.axis = 1.25, xlim = c(1, 94), cex.names = 1, yaxt = "n");
    axis(side = 2, at = c(0, 2, 4, 6, 8, 10, 12, 14), cex.axis = 1.5);
    box(lwd = 2);
    par(new = TRUE);
    y1     <- dat[1:S_s,6] / (dat[1:S_s,5]);
    x1     <- seq(from = 2.132, to = 15.1112, length = S_s);
    plot(x = x1, y = y1, xaxt = "n", yaxt = "n", lwd = 2, ylim = c(0, 1.1),
         xlab = "", ylab = "", type = "b", xlim = c(2, 15), pch = 20, 
         cex = 1, col = "black", yaxs="i");
    points(x = x1, y = y1, lwd = 2, type = "l", col = "black");
    axis(side = 4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), cex.axis = 1.5);
    legend("topleft", c(expression(paste(sigma[gamma]^2," = 0")), 
                        expression(paste(sigma[gamma]^2," = 0"))), 
           pch=15, col=c("dodgerblue4","firebrick"), cex = 1.5, horiz = TRUE);
    #=================================
    mtext(side = 1, text = "System size (S)", cex = 2, outer = TRUE, 
          line = 3.0);
    mtext(side = 2, text = "Ln number of stable systems", cex = 2, 
          outer = TRUE, line = 3.5);
    mtext(side = 4, text = expression(
        paste("Pr. of systems stable due to Var(",gamma,")")),
        cex = 2, outer = TRUE, line = 3.5);
}
