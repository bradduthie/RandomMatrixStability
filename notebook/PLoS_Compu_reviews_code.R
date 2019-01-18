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




# Reviewer 2 specific comment 3
Ec_calc <- function(A0, E1, E2){
    S   <- dim(A0)[1];
    smm <- 0;
    for(i in 1:S){
        for(j in 1:S){
            if(i != j){
                smm <- smm + A0[i,j] * A0[j, i] - (E1*E1 / E2*E2);
            }
        }
    }
    EC <- (1 / S*(S - 1) * E2*E2) * smm;
    return(EC);
}
#----------------------------
S         <- 28;
A0rE1     <- NULL;
A0rE2     <- NULL;
A0rEC     <- NULL;
A1rE1     <- NULL;
A1rE2     <- NULL;
A1rEC     <- NULL;
A0_eig    <- NULL;
A1_eig    <- NULL;
iter      <- 10000;
while(iter > 0){
    r_vec     <- rnorm(n =S, mean = 0, sd = 0.4);
    A0_dat    <- rnorm(n =S *S, mean = 0, sd = sigma);
    A0        <- matrix(data = A0_dat, nrow =S, 
                        ncol =S);
    C_dat     <- rbinom(n =S *S, size = 1, prob = C);
    C_mat     <- matrix(data = C_dat, nrow =S, 
                        ncol =S);
    A0        <- A0 * C_mat;
    diag(A0)  <- -1;
    gam1      <- runif(n =S, min = 0, max = 2);
    A1        <- A0 * gam1;
    A0        <- A0 * mean(gam1);
    A0_stb    <- max(Re(eigen(A0)$values));
    A1_stb    <- max(Re(eigen(A1)$values));
    Ad0       <- A0;
    diag(Ad0) <- NA;
    Ad1       <- A1;
    diag(Ad1) <- NA;
    A0E1      <- mean(Ad0, na.rm = TRUE);
    A0E2      <- var(as.vector(Ad0), na.rm = TRUE);
    A0EC      <- Ec_calc(A0, E1, E2);
    A1E1      <- mean(Ad1, na.rm = TRUE);
    A1E2      <- var(as.vector(Ad1), na.rm = TRUE);
    A1EC      <- Ec_calc(A1, E1, E2);
    A0rE1     <- c(A0rE1, A0E1);
    A0rE2     <- c(A0rE2, A0E2);
    A0rEC     <- c(A0rEC, A0EC);
    A1rE1     <- c(A1rE1, A1E1);
    A1rE2     <- c(A1rE2, A1E2);
    A1rEC     <- c(A1rEC, A1EC);
    A0_eig    <- c(A0_eig, A0_stb);
    A1_eig    <- c(A1_eig, A1_stb);
    iter      <- iter - 1;
}



resilA0  <- A0_eig^-1;
resilA1  <- A1_eig^-1;
ch_resil <- resilA1 - resilA0;





