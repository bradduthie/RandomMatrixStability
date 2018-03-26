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
        all_res[i, 1] <- i;
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






