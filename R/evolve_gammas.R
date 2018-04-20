Evo_rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, C = 1){
    tot_res <- NULL;
    fea_res <- NULL;
    for(i in 2:max_sp){
        nn             <- i;
        A1_stt         <- 0;
        A2_stt         <- 0;
        A1_fet         <- 0;
        A2_fet         <- 0;
        iter           <- iters;
        tot_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 3);
        fea_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 2);
        while(iter > 0){
            r_vec    <- rnorm(n = i, mean = 0, sd = rmx);
            A0_dat   <- rnorm(n = i * i, mean = 0, sd = 0.4);
            A0       <- matrix(data = A0_dat, nrow = i, ncol = i);
            A0       <- species_interactions(mat = A0, type = int_type);
            C_dat    <- rbinom(n = i * i, size = 1, prob = C);
            C_mat    <- matrix(data = C_dat, nrow = i, ncol = i);
            A0       <- A0 * C_mat;
            diag(A0) <- -1;
            gam0     <- make_gammas(nn = i, distribution = 0, sdd = 1);
            gam1     <- make_gammas(nn = i, distribution = 1, sdd = 1);
            A1       <- A0 * gam1;
            A0       <- A0 * gam0;
            A0_stb   <- max(Re(eigen(A0)$values)) < 0;
            A1_stb   <- rand_mat_ga(A1);
            A0_fea   <- min(-1*solve(A0) %*% r_vec) > 0;
            A1_fea   <- min(-1*solve(A1) %*% r_vec) > 0;
            if(A0_stb == TRUE){
                tot_res[[i-1]][iter, 1] <- 1;
            }
            if(A1_stb == TRUE){
                tot_res[[i-1]][iter, 2] <- 1;
            }
            if(A0_fea == TRUE){
                fea_res[[i-1]][iter, 1] <- 1;
            }
            if(A1_fea == TRUE){
                fea_res[[i-1]][iter, 2] <- 1;
            }
            iter    <- iter - 1;
        }
        print(i);
    }
    all_res <- summarise_randmat(tot_res = tot_res, fea_res = fea_res);
    return(all_res);
}


summarise_randmat <- function(tot_res, fea_res){
    sims    <- length(tot_res);
    all_res <- matrix(data = 0, nrow = sims, ncol = 10);
    for(i in 1:sims){
        unstables <- tot_res[[i]][,1] == FALSE & tot_res[[i]][,2] == FALSE;
        stables   <- tot_res[[i]][,1] == TRUE  & tot_res[[i]][,2] == TRUE;
        unstabled <- tot_res[[i]][,1] == TRUE  & tot_res[[i]][,2] == FALSE;
        stabled   <- tot_res[[i]][,1] == FALSE & tot_res[[i]][,2] == TRUE;
        non_feas  <- fea_res[[i]][,1] == FALSE & fea_res[[i]][,2] == FALSE;
        feasibl   <- fea_res[[i]][,1] == TRUE  & fea_res[[i]][,2] == TRUE;
        unfeased  <- fea_res[[i]][,1] == TRUE  & fea_res[[i]][,2] == FALSE;
        feased    <- fea_res[[i]][,1] == FALSE & fea_res[[i]][,2] == TRUE;
        foundd    <- tot_res[[i]][,3] == TRUE;
        all_res[i, 1]  <- i + 1;
        all_res[i, 2]  <- sum(unstables);
        all_res[i, 3]  <- sum(stables);
        all_res[i, 4]  <- sum(unstabled);
        all_res[i, 5]  <- sum(stabled);
        all_res[i, 6]  <- sum(non_feas);
        all_res[i, 7]  <- sum(feasibl);
        all_res[i, 8]  <- sum(unfeased);
        all_res[i, 9]  <- sum(feased);
        all_res[i, 10] <- sum(foundd);
    }
    return(all_res);
}

rand_mat_ga <- function(A1, max_it = 20, converg = 0.01){
    nn       <- dim(A1)[1];
    rind     <- runif(n = nn*1000, min = 0, max = 1);
    inds     <- matrix(data = rind, nrow = 1000, ncol = nn);
    lastf    <- -10;
    ccrit    <- 10;
    find_st  <- 0;
    iter     <- max_it;
    while(iter > 0 & find_st < 1 & ccrit > converg){
        ivar  <- rep(x = 0, length = dim(inds)[1]);
        ifit  <- rep(x = 0, length = dim(inds)[1]);
        isst  <- rep(x = 0, length = dim(inds)[1]);
        for(i in 1:dim(inds)[1]){
            ifit[i] <- -1*max(Re(eigen(inds[i,]*A1)$values));
            ivar[i] <- var(inds[i,]);
            isst[i] <- max(Re(eigen(inds[i,]*A1)$values)) < 0;
        }
        most_fit <- order(ifit, decreasing = TRUE)[1:20];
        parents  <- inds[most_fit,];
        new_gen  <- matrix(data = t(parents), nrow = 1000, ncol = nn, 
                           byrow = TRUE);
        mu_dat   <- rbinom(n = nn*1000, size = 1, prob = 0.2);
        mu_dat2  <- rnorm(n = nn*1000, mean = 0, sd = 0.02);
        mu_dat2[mu_dat2 < 0] <- -mu_dat2[mu_dat2 < 0];
        mu_dat2[mu_dat2 > 1] <- 1;
        mu_dat3  <- mu_dat * mu_dat2;
        mu_mat   <- matrix(data = mu_dat3, nrow = 1000, ncol = nn);
        new_gen  <- new_gen + mu_mat;
        new_gen  <- crossover(inds = new_gen, pr = 0.1);
        inds     <- new_gen;
        find_st  <- max(isst);
        newf     <- mean(ifit);
        ccrit    <- newf - lastf;
        lastf    <- newf;
        iter     <- iter - 1;
    }
    if(find_st == 1){
        s_row <- which(isst == 1)[1];
        writt <- c(nn, inds[s_row,]);
        cat(writt, file = "Evo_out.txt", append = TRUE);
        cat("\n", file = "Evo_out.txt", append = TRUE);
    }
    return(find_st);
}


crossover <- function(inds, pr = 0.1){
    crossed <- floor(dim(inds)[1] * pr);
    cross1  <- sample(x = 1:dim(inds)[1], size = crossed);
    cross2  <- sample(x = 1:dim(inds)[1], size = crossed);
    for(i in 1:length(cross1)){
        fromv   <- sample(x = 1:dim(inds)[2], size = 1);
        tov     <- sample(x = 1:dim(inds)[2], size = 1);
        temp                   <- inds[cross1[i],fromv:tov];
        inds[cross1[i],fromv:tov] <- inds[cross2[i],fromv:tov];
        inds[cross2[i],fromv:tov] <- temp;
    }
    return(inds);
}


