rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, eps_max = 999){
    tot_res <- NULL;
    fea_res <- NULL;
    nea_res <- NULL;
    for(i in 2:max_sp){
        nn             <- i;
        A1_stt         <- 0;
        A2_stt         <- 0;
        A1_fet         <- 0;
        A2_fet         <- 0;
        iter           <- iters;
        tot_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 2);
        fea_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 2);
        nea_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 2);
        while(iter > 0){
            r_vec    <- rnorm(n = nn, mean = 0, sd = rmx);
            A1_dat   <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
            A1       <- matrix(data = A1_dat, nrow = nn, ncol = nn);
            A1       <- species_interactions(mat = A1, type = int_type);
            diag(A1) <- -1;
            epsil    <- runif(n = nn, min = 1, max = eps_max);
            eps_dat  <- rep(x = epsil, times = nn);
            eps_mat  <- matrix(data = epsil, nrow = nn, ncol = nn, 
                               byrow = FALSE);
            avg_mat  <- matrix(data = mean(1:eps_max), nrow = nn, ncol = nn, 
                               byrow = FALSE);
            A2       <- A1 * eps_mat;
            A1       <- A1 * avg_mat;
            #r_vec    <- r_vec * epsil / 100;
            A1_stb   <- max(Re(eigen(A1)$values)) < 0;
            A2_stb   <- max(Re(eigen(A2)$values)) < 0;
            A1_fea   <- min(-1*solve(A1) %*% r_vec) > 0;
            A2_fea   <- min(-1*solve(A2) %*% r_vec) > 0;
            A1_nea   <- max(-1*solve(A1) %*% r_vec) < 0;
            A2_nea   <- max(-1*solve(A2) %*% r_vec) < 0;
            if(A1_stb == TRUE){
                tot_res[[i-1]][iter, 1] <- 1;
            }
            if(A2_stb == TRUE){
                tot_res[[i-1]][iter, 2] <- 1;
            }
            if(A1_fea == TRUE){
                fea_res[[i-1]][iter, 1] <- 1;
            }
            if(A2_fea == TRUE){
                fea_res[[i-1]][iter, 2] <- 1;
            }
            if(A1_nea == TRUE){
                nea_res[[i-1]][iter, 1] <- 1;
            }
            if(A2_nea == TRUE){
                nea_res[[i-1]][iter, 2] <- 1;
            }
            iter    <- iter - 1;
        }
        print(i);
    }
    all_res <- summarise_randmat(tot_res = tot_res, fea_res = fea_res, 
                                 nea_res = nea_res);
    return(all_res);
}


species_interactions <- function(mat, type = 0){
    if(type == 1){
        mat[mat > 0] <- -1*mat[mat > 0];
    }
    if(type == 2){
        mat[mat < 0] <- -1*mat[mat < 0];
    }
    if(type == 3){
        for(i in 1:dim(mat)[1]){
            for(j in 1:dim(mat)[2]){
                if(mat[i, j] * mat[j, i] > 0){
                    mat[j, i] <- -1 * mat[j, i];
                }
            }
        }
    }
    return(mat);
}


summarise_randmat <- function(tot_res, fea_res, nea_res){
    sims    <- length(tot_res);
    all_res <- matrix(data = 0, nrow = sims, ncol = 13);
    for(i in 1:sims){
        unstables <- tot_res[[i]][,1] == FALSE & tot_res[[i]][,2] == FALSE;
        stables   <- tot_res[[i]][,1] == TRUE  & tot_res[[i]][,2] == TRUE;
        unstabled <- tot_res[[i]][,1] == TRUE  & tot_res[[i]][,2] == FALSE;
        stabled   <- tot_res[[i]][,1] == FALSE & tot_res[[i]][,2] == TRUE;
        non_feas  <- fea_res[[i]][,1] == FALSE & fea_res[[i]][,2] == FALSE;
        feasibl   <- fea_res[[i]][,1] == TRUE  & fea_res[[i]][,2] == TRUE;
        unfeased  <- fea_res[[i]][,1] == TRUE  & fea_res[[i]][,2] == FALSE;
        feased    <- fea_res[[i]][,1] == FALSE & fea_res[[i]][,2] == TRUE;
        non_neas  <- nea_res[[i]][,1] == FALSE & nea_res[[i]][,2] == FALSE;
        neasibl   <- nea_res[[i]][,1] == TRUE  & nea_res[[i]][,2] == TRUE;
        unneased  <- nea_res[[i]][,1] == TRUE  & nea_res[[i]][,2] == FALSE;
        neased    <- nea_res[[i]][,1] == FALSE & nea_res[[i]][,2] == TRUE;
        all_res[i, 1]  <- i + 1;
        all_res[i, 2]  <- sum(unstables);
        all_res[i, 3]  <- sum(stables);
        all_res[i, 4]  <- sum(unstabled);
        all_res[i, 5]  <- sum(stabled);
        all_res[i, 6]  <- sum(non_feas);
        all_res[i, 7]  <- sum(feasibl);
        all_res[i, 8]  <- sum(unfeased);
        all_res[i, 9]  <- sum(feased);
        all_res[i, 10] <- sum(non_neas);
        all_res[i, 11] <- sum(neasibl);
        all_res[i, 12] <- sum(unneased);
        all_res[i, 13] <- sum(neased);
    }
    return(all_res);
}


# Get:
# (1) proportion stable made unstable
# (2) proportion unstable made stable
# (3) proportion feasible made infeasible
# (4) proportion infeasible made feasible


random      <- read.csv(file = "sim_results/random_community.csv");
competition <- read.csv(file = "sim_results/competitor_community.csv");
mutualism   <- read.csv(file = "sim_results/mutualist_community.csv");
pred_prey   <- read.csv(file = "sim_results/predator_prey_community.csv");



dat <- pred_prey;

pr_destab <- dat[,4] / (dat[,4] + dat[,3]);
pr_stabld <- dat[,5] / (dat[,5] + dat[,2]);
pr_defeas <- dat[,8] / (dat[,8] + dat[,7]);
pr_feased <- dat[,9] / (dat[,9] + dat[,6]);
Sp        <- dat[,1];

plot(dat[,1], pr_destab, type = "l", lwd = 2, ylim = c(0, 1));
points(dat[,1], pr_stabld, type = "l", lwd = 2, col = "red");
points(dat[,1], pr_defeas, type = "l", lwd = 2, col = "black", lty = "dashed");
points(dat[,1], pr_feased, type = "l", lwd = 2, col = "red", lty = "dashed");


pr_destab <- dat[,4] / (dat[,4] + dat[,3]);
pr_stabld <- dat[,5] / (dat[,5] + dat[,2]);
pr_defeas <- dat[,8] / (dat[,8] + dat[,7]);
pr_feased <- dat[,9] / (dat[,9] + dat[,6]);



# ==============================================================================
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1));
# ==============================================================================
dat      <- random;
stab_nov <- dat[,3] + dat[,4];
stab_var <- dat[,3] + dat[,5];
feas_nov <- dat[,7] + dat[,8];
feas_var <- dat[,7] + dat[,9];
plot(x = Sp, y = stab_nov / 100000, type = "l", lwd = 2, ylim = c(0,1),
     cex.lab = 1.5, cex.axis = 1.5, xlab = "Species Number (S)",
     ylab = "Proportion stable or feasible");
points(x = Sp, y = stab_var / 100000, type = "l", lwd = 2, lty = "dashed");
points(x = Sp, y = feas_nov / 100000, type = "l", lwd = 2, col = "red");
points(x = Sp, y = feas_var / 100000, type = "l", lwd = 2, col = "red",
       lty = "dashed");
legend(x = 40, y = 0.9, col = c("black", "black", "red", "red"), cex = 1.2,
       legend = c(expression(paste("Stable, no Var(",gamma,")")),
                             expression(paste("Stable, Var(",gamma,")")),
                             expression(paste("Feasible, no Var(",gamma,")")),
                             expression(paste("Feasible, Var(",gamma,")"))),
       lty = c("solid", "dashed", "solid", "dashed"), lwd = 2);
text(x = 60, y = 0.98, labels = "Random", cex = 2);
# ==============================================================================
dat      <- competition;
stab_nov <- dat[,3] + dat[,4];
stab_var <- dat[,3] + dat[,5];
feas_nov <- dat[,7] + dat[,8];
feas_var <- dat[,7] + dat[,9];
plot(x = Sp, y = stab_nov / 100000, type = "l", lwd = 2, ylim = c(0,1),
     cex.lab = 1.5, cex.axis = 1.5, xlab = "Species Number (S)",
     ylab = "Proportion stable or feasible");
points(x = Sp, y = stab_var / 100000, type = "l", lwd = 2, lty = "dashed");
points(x = Sp, y = feas_nov / 100000, type = "l", lwd = 2, col = "red");
points(x = Sp, y = feas_var / 100000, type = "l", lwd = 2, col = "red",
       lty = "dashed");
text(x = 60, y = 0.98, labels = "Competitive", cex = 2);
# ==============================================================================
dat      <- mutualism;
stab_nov <- dat[,3] + dat[,4];
stab_var <- dat[,3] + dat[,5];
feas_nov <- dat[,7] + dat[,8];
feas_var <- dat[,7] + dat[,9];
plot(x = Sp, y = stab_nov / 100000, type = "l", lwd = 2, ylim = c(0,1),
     cex.lab = 1.5, cex.axis = 1.5, xlab = "Species Number (S)",
     ylab = "Proportion stable or feasible");
points(x = Sp, y = stab_var / 100000, type = "l", lwd = 2, lty = "dashed");
points(x = Sp, y = feas_nov / 100000, type = "l", lwd = 2, col = "red");
points(x = Sp, y = feas_var / 100000, type = "l", lwd = 2, col = "red",
       lty = "dashed");
text(x = 60, y = 0.98, labels = "Mutualist", cex = 2);
# ==============================================================================
dat      <- pred_prey;
stab_nov <- dat[,3] + dat[,4];
stab_var <- dat[,3] + dat[,5];
feas_nov <- dat[,7] + dat[,8];
feas_var <- dat[,7] + dat[,9];
plot(x = Sp, y = stab_nov / 100000, type = "l", lwd = 2, ylim = c(0,1),
     cex.lab = 1.5, cex.axis = 1.5, xlab = "Species Number (S)",
     ylab = "Proportion stable or feasible");
points(x = Sp, y = stab_var / 100000, type = "l", lwd = 2, lty = "dashed");
points(x = Sp, y = feas_nov / 100000, type = "l", lwd = 2, col = "red");
points(x = Sp, y = feas_var / 100000, type = "l", lwd = 2, col = "red",
       lty = "dashed");
text(x = 60, y = 0.98, labels = "Predator-prey", cex = 2);
# ==============================================================================





mat <- matrix(data = rnorm(n = 4, mean = 0, sd = 0.4), nrow = 2);
diag(mat) <- -1;
plot(mat, xlim = c(-0.5, 0.5), ylim = c(-1.25, 1.25), asp = 1, pch = 20)
points(Re(eigen(mat)$vector));
abline(h = 0, lty = "dotted", lwd = 0.8);
abline(v = 0, lty = "dotted", lwd = 0.8);






plot(x = 0, y = 0, type = "n", xlim = c(-0.3, 0.3), ylim = c(-2.5, 2.5),
     xlab = "Real component", ylab = "Imaginary component", cex.lab = 1.25);
abline(h = 0, lty = "dotted", lwd = 0.8);
abline(v = 0, lty = "dotted", lwd = 0.8);

iter     <- 10000;
while(iter > 0){
    mat1       <- matrix(data = rnorm(n = 16, mean = 0, sd = 0.4), nrow = 4);
    diag(mat1) <- -1;
    mat2       <- mat1;
    mat2[1,]   <- mat2[1,] * runif(n = 1, min = 1, max = 1000);
    mat1_e     <- eigen(mat1)$values;
    mat2_e     <- eigen(mat2)$values;
    stab1      <- max(Re(mat1_e)) < 0;
    stab2      <- max(Re(mat2_e)) < 0;
    if(stab1 == TRUE & stab2 == FALSE){
        arrows(x0 = Re(mat1_e), x1 = Re(mat2_e), y0 = Im(mat1_e), y1 = Im(mat2_e),
               length = 0.05, lwd = 2, col = "red");
        #break;
    }
    if(stab1 == FALSE & stab2 == TRUE){
        arrows(x0 = Re(mat1_e), x1 = Re(mat2_e), y0 = Im(mat1_e), y1 = Im(mat2_e),
               length = 0.05, lwd = 2, col = "blue");
        #break;
    }
    iter <- iter - 1;
}



get_rand_eigs <- function(int_type = 0, iter = 1000){
    nn             <- 20;
    eigres         <- matrix(data = 0, ncol = 82, nrow = iter);
    eigres[,1]     <- int_type;
    while(iter > 0){
        r_vec     <- runif(n = nn, min = 1, max = rmx);
        A1_dat    <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
        A1        <- matrix(data = A1_dat, nrow = nn, ncol = nn);
        A1        <- species_interactions(mat = A1, type = int_type);
        diag(A1)  <- -1;
        epsil     <- runif(n = nn, min = 1, max = eps_max);
        eps_dat   <- rep(x = epsil, times = nn);
        eps_mat   <- matrix(data = epsil, nrow = nn, ncol = nn, 
                            byrow = FALSE);
        avg_mat   <- matrix(data = mean(1:eps_max), nrow = nn, ncol = nn, 
                            byrow = FALSE);
        A2        <- A1 * eps_mat;
        A1        <- A1 * avg_mat;
        mat1_e    <- eigen(A1)$values;
        mat2_e    <- eigen(A2)$values;
        real_mat1 <- Re(mat1_e);
        imag_mat1 <- Im(mat1_e);
        real_mat2 <- Re(mat2_e);
        imag_mat2 <- Im(mat2_e);
        stab1     <- max(real_mat1) < 0;
        stab2     <- max(real_mat2) < 0;
        if(stab1 == FALSE & stab2 == FALSE){
            eigres[iter,2] <- 1;
        }
        if(stab1 == TRUE & stab2 == TRUE){
            eigres[iter,2] <- 2;
        }
        if(stab1 == TRUE & stab2 == FALSE){
            eigres[iter,2] <- 3;
        }
        if(stab1 == FALSE & stab2 == TRUE){
            eigres[iter,2] <- 4;
        }
        eigres[iter,3:22]  <- real_mat1;
        eigres[iter,23:42] <- imag_mat1;
        eigres[iter,43:62] <- real_mat2;
        eigres[iter,63:82] <- imag_mat2;
        iter    <- iter - 1;
    }
    return(eigres);
}

plot_eig_set <- function(eigres){
    par(mfrow = c(2,2), mar = c(5, 5, 1, 1));
    x0 <- min(cbind(eigres[,3:22], eigres[,43:62]));
    x1 <- max(cbind(eigres[,3:22], eigres[,43:62]));
    plot(x = 0, y = 0, type = "n", xlim = c(x0, x1), ylim = c(-1500, 1500),
         xlab = "Real component", ylab = "Imaginary component", cex.lab = 1.25);
    abline(h = 0, lty = "dotted", lwd = 0.8);
    abline(v = 0, lty = "dotted", lwd = 0.8);
    if(1 %in% eigres[,2]){
        row <- which(eigres[,2] == 1)[1];
        m1r <- eigres[row,3:22];
        m1i <- eigres[row,23:42];
        m2r <- eigres[row,43:62];
        m2i <- eigres[row,63:82];
        points(x = m1r, y = m1i, pch = 20, cex = 1.5, col = "black");
        points(x = m2r, y = m2i, pch = 20, cex = 1.5, col = "red");
    }
    text(x = x0 + (x1 - x0)*0.2, y = 1500, "Unstable", cex = 1.5);
    plot(x = 0, y = 0, type = "n", xlim = c(x0, x1), ylim = c(-1500, 1500),
         xlab = "Real component", ylab = "Imaginary component", cex.lab = 1.25);
    abline(h = 0, lty = "dotted", lwd = 0.8);
    abline(v = 0, lty = "dotted", lwd = 0.8);
    if(2 %in% eigres[,2]){
        row <- which(eigres[,2] == 2)[1];
        m1r <- eigres[row,3:22];
        m1i <- eigres[row,23:42];
        m2r <- eigres[row,43:62];
        m2i <- eigres[row,63:82];
        points(x = m1r, y = m1i, pch = 20, cex = 1.5, col = "black");
        points(x = m2r, y = m2i, pch = 20, cex = 1.5, col = "red");
    }
    text(x = x0 + (x1 - x0)*0.2, y = 1500, "Stable", cex = 1.5);
    plot(x = 0, y = 0, type = "n", xlim = c(x0, x1), ylim = c(-1500, 1500),
         xlab = "Real component", ylab = "Imaginary component", cex.lab = 1.25);
    abline(h = 0, lty = "dotted", lwd = 0.8);
    abline(v = 0, lty = "dotted", lwd = 0.8);
    if(3 %in% eigres[,2]){
        row <- which(eigres[,2] == 3)[1];
        m1r <- eigres[row,3:22];
        m1i <- eigres[row,23:42];
        m2r <- eigres[row,43:62];
        m2i <- eigres[row,63:82];
        points(x = m1r, y = m1i, pch = 20, cex = 1.5, col = "black");
        points(x = m2r, y = m2i, pch = 20, cex = 1.5, col = "red");
    }
    text(x = x0 + (x1 - x0)*0.2, y = 1500, "Destabilised", cex = 1.5);
    plot(x = 0, y = 0, type = "n", xlim =c(x0, x1), ylim = c(-1500, 1500),
         xlab = "Real component", ylab = "Imaginary component", cex.lab = 1.25);
    abline(h = 0, lty = "dotted", lwd = 0.8);
    abline(v = 0, lty = "dotted", lwd = 0.8);
    if(4 %in% eigres[,2]){
        row <- which(eigres[,2] == 4)[1];
        m1r <- eigres[row,3:22];
        m1i <- eigres[row,23:42];
        m2r <- eigres[row,43:62];
        m2i <- eigres[row,63:82];
        points(x = m1r, y = m1i, pch = 20, cex = 1.5, col = "black");
        points(x = m2r, y = m2i, pch = 20, cex = 1.5, col = "red");
    }
    text(x = x0 + (x1 - x0)*0.2, y = 1500, "Stabilised", cex = 1.5);
}

get_rand_feas <- function(int_type = 0, iter = 1000){
    nn             <- 20;
    eigres         <- matrix(data = 0, ncol = 42, nrow = iter);
    eigres[,1]     <- int_type;
    while(iter > 0){
        r_vec     <- runif(n = nn, min = 1, max = rmx);
        A1_dat    <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
        A1        <- matrix(data = A1_dat, nrow = nn, ncol = nn);
        A1        <- species_interactions(mat = A1, type = int_type);
        diag(A1)  <- -1;
        avg_mat   <- matrix(data = mean(1:eps_max), nrow = nn, ncol = nn, 
                            byrow = FALSE);
        A1        <- A1 * avg_mat;
        mat1_e    <- eigen(A1)$values;
        real_mat1 <- Re(mat1_e);
        imag_mat1 <- Im(mat1_e);
        stab1     <- max(real_mat1) < 0;
        A1_fea    <- min(-1*solve(A1) %*% r_vec) > 0;
        if(stab1 == TRUE & A1_fea == FALSE){
            eigres[iter,2] <- 1;
        }
        if(stab1 == TRUE & A1_fea == TRUE){
            eigres[iter,2] <- 2;
        }
        eigres[iter,3:22]  <- real_mat1;
        eigres[iter,23:42] <- imag_mat1;
        iter    <- iter - 1;
    }
    return(eigres);
}


plot_feas_set <- function(eigres){
    par(mfrow = c(1,2), mar = c(5, 5, 1, 1));
    x0 <- min(eigres[,3:22]);
    x1 <- max(eigres[,3:22]);
    plot(x = 0, y = 0, type = "n", xlim = c(x0, x1), ylim = c(-1500, 1500),
         xlab = "Real component", ylab = "Imaginary component", cex.lab = 1.25);
    abline(h = 0, lty = "dotted", lwd = 0.8);
    abline(v = 0, lty = "dotted", lwd = 0.8);
    if(1 %in% eigres[,2]){
        row <- which(eigres[,2] == 1)[1];
        m1r <- eigres[row,3:22];
        m1i <- eigres[row,23:42];
        points(x = m1r, y = m1i, pch = 20, cex = 1.5, col = "black");
    }
    text(x = x0 + (x1 - x0)*0.2, y = 1500, "Stable", cex = 1.5);
    plot(x = 0, y = 0, type = "n", xlim = c(x0, x1), ylim = c(-1500, 1500),
         xlab = "Real component", ylab = "Imaginary component", cex.lab = 1.25);
    abline(h = 0, lty = "dotted", lwd = 0.8);
    abline(v = 0, lty = "dotted", lwd = 0.8);
    if(2 %in% eigres[,2]){
        row <- which(eigres[,2] == 2)[1];
        m1r <- eigres[row,3:22];
        m1i <- eigres[row,23:42];
        points(x = m1r, y = m1i, pch = 20, cex = 1.5, col = "black");
    }
    text(x = x0 + (x1 - x0)*0.2, y = 1500, "Feasible", cex = 1.5);
}


eres <- get_rand_eigs(int_type = 0, iter = 10000);
plot_eig_set(eres);






get_var_change <- function(int_type = 0, iter = 1000){
    m_changes      <- NULL;
    nn             <- 20;
    eigres         <- matrix(data = 0, ncol = 82, nrow = iter);
    eigtests       <- matrix(data = 0, nrow = iter, ncol = 21);
    eigres[,1]     <- int_type;
    while(iter > 0){
        r_vec     <- runif(n = nn, min = 1, max = rmx);
        A1_dat    <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
        A1        <- matrix(data = A1_dat, nrow = nn, ncol = nn);
        A1        <- species_interactions(mat = A1, type = int_type);
        diag(A1)  <- -1;
        epsil     <- runif(n = nn, min = 1, max = eps_max);
        eps_dat   <- rep(x = epsil, times = nn);
        eps_mat   <- matrix(data = epsil, nrow = nn, ncol = nn, 
                            byrow = FALSE);
        avg_mat   <- matrix(data = mean(1:eps_max), nrow = nn, ncol = nn, 
                            byrow = FALSE);
        A2        <- A1 * eps_mat;
        A1        <- A1 * avg_mat;
        mat1_e    <- eigen(A1)$values;
        mat2_e    <- eigen(A2)$values;
        real_mat1 <- Re(mat1_e);
        imag_mat1 <- Im(mat1_e);
        real_mat2 <- Re(mat2_e);
        imag_mat2 <- Im(mat2_e);
        stab1     <- max(real_mat1) < 0;
        stab2     <- max(real_mat2) < 0;
        mchange   <- solve(t(A1) %*% A1) %*% t(A1) %*% A2;
        m_changes[[iter]] <- mchange;
        if(stab1 == FALSE & stab2 == FALSE){
            eigres[iter,2]    <- 1;
            eigtests[iter,1]  <- 1;
        }
        if(stab1 == TRUE & stab2 == TRUE){
            eigres[iter,2] <- 2;
            eigtests[iter,1]  <- 2;
        }
        if(stab1 == TRUE & stab2 == FALSE){
            eigres[iter,2] <- 3;
            eigtests[iter,1]  <- 3;
        }
        if(stab1 == FALSE & stab2 == TRUE){
            eigres[iter,2] <- 4;
            eigtests[iter,1]  <- 4;
        }
        eigtests[iter,2:21] <- eigen(mchange)$values;
        eigres[iter,3:22]   <- real_mat1;
        eigres[iter,23:42]  <- imag_mat1;
        eigres[iter,43:62]  <- real_mat2;
        eigres[iter,63:82]  <- imag_mat2;
        iter    <- iter - 1;
    }
    return(list(eigres = eigres, mchanges = m_changes, eigtests = eigtests));
}

plot_var_change <- function(dat){
    eigres   <- dat$eigres;
    mchanges <- dat$mchanges;
    par(mfrow = c(2,2), mar = c(2, 2, 1, 1));
    if(1 %in% eigres[,2]){
        use_it <- which(eigres[,2] == 1)[1];
        plot(mchanges[[use_it]], xlab = "", ylab = "");
        abline(h = 0, lty = "dotted", lwd = 0.8);
        abline(v = 0, lty = "dotted", lwd = 0.8);
    }
    if(2 %in% eigres[,2]){
        use_it <- which(eigres[,2] == 2)[1];
        plot(mchanges[[use_it]], xlab = "", ylab = "");
        abline(h = 0, lty = "dotted", lwd = 0.8);
        abline(v = 0, lty = "dotted", lwd = 0.8);
    }
    if(3 %in% eigres[,2]){
        use_it <- which(eigres[,2] == 3)[1];
        plot(mchanges[[use_it]], xlab = "", ylab = "");
        abline(h = 0, lty = "dotted", lwd = 0.8);
        abline(v = 0, lty = "dotted", lwd = 0.8);
    }
    if(4 %in% eigres[,2]){
        use_it <- which(eigres[,2] == 4)[1];
        plot(mchanges[[use_it]], xlab = "", ylab = "");
        abline(h = 0, lty = "dotted", lwd = 0.8);
        abline(v = 0, lty = "dotted", lwd = 0.8);
    }
}












dat      <- random;


Sp       <- dat[,1];
stab_nov <- dat[,3]  + dat[,4];
stab_var <- dat[,3]  + dat[,5];
feas_nov <- dat[,7]  + dat[,8];
feas_var <- dat[,7]  + dat[,9];
#neas_nov <- dat[,11] + dat[,12];
#neas_var <- dat[,11]  + dat[,13];
plot(x = Sp, y = stab_nov / 100000, type = "l", lwd = 2, ylim = c(0,0.2),
     cex.lab = 1.5, cex.axis = 1.5, xlab = "Species Number (S)",
     ylab = "Proportion stable or feasible");
points(x = Sp, y = stab_var / 100000, type = "l", lwd = 2, lty = "dashed");
points(x = Sp, y = feas_nov / 100000, type = "l", lwd = 2, col = "red");
points(x = Sp, y = feas_var / 100000, type = "l", lwd = 2, col = "red",
       lty = "dashed");
#points(x = Sp, y = neas_nov / 100000, type = "l", lwd = 2, col = "blue");
#points(x = Sp, y = neas_var / 100000, type = "l", lwd = 2, col = "blue",
#       lty = "dashed");
legend(x = 40, y = 0.9, col = c("black", "black", "red", "red"), cex = 1.2,
       legend = c(expression(paste("Stable, no Var(",gamma,")")),
                  expression(paste("Stable, Var(",gamma,")")),
                  expression(paste("Feasible, no Var(",gamma,")")),
                  expression(paste("Feasible, Var(",gamma,")"))),
       lty = c("solid", "dashed", "solid", "dashed"), lwd = 2);
text(x = 60, y = 0.98, labels = "Random", cex = 2);













