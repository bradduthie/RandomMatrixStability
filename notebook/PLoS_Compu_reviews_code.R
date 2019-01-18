# Reviewer 2 specific comment 1
SD_mn_row_vals <- function(S = 10, iters = 1000, sigma = 0.4){
    sims  <- rep(x = 0, times = iters);
    for (i in 1:iters){
        A0_dat   <- rnorm(n = S * S, mean = 0, sd = sigma);
        A0       <- matrix(data = A0_dat, nrow = S); 
        diag(A0) <- -1;
        rowmns   <- apply(X = A0, MAR = 1, FUN = mean);
        sims[i]  <- sd(rowmns);
    }
    return(mean(sims));
}
eg_run <- SD_mn_row_vals(S = 50);


# Reviewer 2 specific comment 4
S       <- 30;
A0_mx   <- NULL;
A1_mx   <- NULL;
iter    <- 10000;
while(iter > 0){
    r_vec    <- rnorm(n =S, mean = 0, sd = 0.4);
    A0_dat   <- rnorm(n =S *S, mean = 0, sd = sigma);
    A0       <- matrix(data = A0_dat, nrow =S, 
                       ncol =S);
    C_dat    <- rbinom(n =S *S, size = 1, prob = C);
    C_mat    <- matrix(data = C_dat, nrow =S, 
                       ncol =S);
    A0       <- A0 * C_mat;
    diag(A0) <- -1;
    gam1     <- runif(n =S, min = 0, max = 2);
    A1       <- A0 * gam1;
    A0       <- A0 * mean(gam1);
    A0_stb   <- max(Re(eigen(A0)$values));
    A1_stb   <- max(Re(eigen(A1)$values));
    iter     <- iter - 1;
    A0_mx[iter] <- A0_stb;
    A1_mx[iter] <- A1_stb;
}
unstable_A0 <- A0_mx[A0_mx > 0];
unstable_A1 <- A1_mx[A0_mx > 0];