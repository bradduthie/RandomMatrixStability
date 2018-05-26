#' Make a distribution of gamma values
#' 
#' Makes one of five distributions of gamma values, including a constant value
#' (0), a uniform distribution (1), an exponential distribution (2), a beta
#' distribution (3), or a gamma distribution (4).
#'
#'@return A vector with elements sampled froma set distribution
#'@param nn The number of elements in the returned vector
#'@param distribution The type of distribution to return
#'@param mn The mean of the vector to be returned
#'@param sdd The standard devation of the vector to be returned (where 
#'applicable)
#'@examples
#'gammas <- make_gammas(nn = 8);
#'@export
make_gammas <- function(nn = 10, distribution = 1, mn = 1, sdd = 1){
    if(distribution == 0){
        dat          <- rep(x = mn, times = nn);
    }
    if(distribution == 1){
        mval         <- 2;
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

#' Produce figure of gamma distributions
#' 
#' Produces a four panel histogram of gamma distributions from the make_gammas
#' funciton, including a uniform distribution, exponential distribution,
#' beta distribution, and gamma distribution
#'
#'@return A four panel plot showing histograms of distributions
#'@param nn The number of samples from distributions to make the histogram
#'@param sdd The standard devation of the histograms to be made
#'@param mn The mean of the histograms to be made
#'@examples
#'show_gammas(nn = 1000);
#'@export
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
