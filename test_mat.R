rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 10, eps_max = 1000){
    tot_res <- NULL;
    fea_res <- NULL;
    for(i in 2:max_sp){
        nn             <- i;
        A1_stt         <- 0;
        A2_stt         <- 0;
        A1_fet         <- 0;
        A2_fet         <- 0;
        iter           <- iters;
        tot_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 2);
        fea_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 2);
        while(iter > 0){
            r_vec    <- runif(n = nn, min = 1, max = rmx);
            A1_dat   <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
            A1       <- matrix(data = A1_dat, nrow = nn, ncol = nn);
            A1       <- species_interactions(mat = A1, type = int_type);
            diag(A1) <- -1;
            epsil    <- runif(n = nn, min = 1, max = eps_max);
            eps_dat  <- rep(x = epsil, times = nn);
            eps_mat  <- matrix(data = epsil, nrow = nn, ncol = nn, 
                               byrow = FALSE);
            A2       <- A1 * eps_mat;
            A1_stb   <- max(Re(eigen(A1)$values)) < 0;
            A2_stb   <- max(Re(eigen(A2)$values)) < 0;
            A1_fea   <- min(-1*solve(A1) %*% r_vec) > 0;
            A2_fea   <- min(-1*solve(A2) %*% r_vec) > 0;
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
            iter    <- iter - 1;
        }
        print(i);
    }
    all_res <- summarise_randmat(tot_res = tot_res, fea_res = fea_res);
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


summarise_randmat <- function(tot_res, fea_res){
    sims    <- length(tot_res);
    all_res <- matrix(data = 0, nrow = sims, ncol = 9);
    for(i in 1:sims){
        unstables <- tot_res[[i]][,1] == FALSE & tot_res[[i]][,2] == FALSE;
        stables   <- tot_res[[i]][,1] == TRUE  & tot_res[[i]][,2] == TRUE;
        unstabled <- tot_res[[i]][,1] == TRUE  & tot_res[[i]][,2] == FALSE;
        stabled   <- tot_res[[i]][,1] == FALSE & tot_res[[i]][,2] == TRUE;
        non_feas  <- fea_res[[i]][,1] == FALSE & fea_res[[i]][,2] == FALSE;
        feasibl   <- fea_res[[i]][,1] == TRUE  & fea_res[[i]][,2] == TRUE;
        unfeased  <- fea_res[[i]][,1] == TRUE  & fea_res[[i]][,2] == FALSE;
        feased    <- fea_res[[i]][,1] == FALSE & fea_res[[i]][,2] == TRUE;
        all_res[i, 1] <- i + 1;
        all_res[i, 2] <- sum(unstables);
        all_res[i, 3] <- sum(stables);
        all_res[i, 4] <- sum(unstabled);
        all_res[i, 5] <- sum(stabled);
        all_res[i, 6] <- sum(non_feas);
        all_res[i, 7] <- sum(feasibl);
        all_res[i, 8] <- sum(unfeased);
        all_res[i, 9] <- sum(feased);
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











