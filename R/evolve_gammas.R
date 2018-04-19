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
        tot_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 3);
        fea_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 2);
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
            A3_stb   <- rand_mat_ga(A1);
            A1       <- A1 * avg_mat;
            
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
            if(A3_stb == 1){
                tot_res[[i-1]][iter, 3] <- 1;
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
        #print(c(iter, newf, lastf, mean(ivar), find_st));
        lastf    <- newf;
        iter     <- iter - 1;
    }
    #if(find_st == 1){
    #findit (the stable one)
    #fileConn<-file("output.txt")
    #writeLines(c("Hello","World"), fileConn)
    #close(fileConn)
    #}
    return(find_st);
}