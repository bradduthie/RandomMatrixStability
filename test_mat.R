################################################################################
################################################################################
################################################################################

summarise_randmat <- function(tot_res, fea_res){
    sims    <- length(tot_res);
    all_res <- matrix(data = 0, nrow = sims, ncol = 37);
    for(i in 1:sims){
        all_res[i, 1]  <- i + 1;
        # Stable and unstable
        all_res[i, 2]  <- sum(tot_res[[i]][,1] == FALSE);
        all_res[i, 3]  <- sum(tot_res[[i]][,1] == TRUE);
        all_res[i, 4]  <- sum(tot_res[[i]][,2] == FALSE);
        all_res[i, 5]  <- sum(tot_res[[i]][,2] == TRUE);
        all_res[i, 6]  <- sum(tot_res[[i]][,3] == FALSE);
        all_res[i, 7]  <- sum(tot_res[[i]][,3] == TRUE);
        all_res[i, 8]  <- sum(tot_res[[i]][,4] == FALSE);
        all_res[i, 9]  <- sum(tot_res[[i]][,4] == TRUE);
        all_res[i, 10] <- sum(tot_res[[i]][,5] == FALSE);
        all_res[i, 11] <- sum(tot_res[[i]][,5] == TRUE);
        # Stabilised and destabilised
        all_res[i, 12] <- sum(tot_res[[i]][,1] == FALSE & 
                                  tot_res[[i]][,2] == TRUE);
        all_res[i, 13] <- sum(tot_res[[i]][,1] == FALSE & 
                                  tot_res[[i]][,3] == TRUE);
        all_res[i, 14] <- sum(tot_res[[i]][,1] == FALSE & 
                                  tot_res[[i]][,4] == TRUE);
        all_res[i, 15] <- sum(tot_res[[i]][,1] == FALSE & 
                                  tot_res[[i]][,5] == TRUE);
        all_res[i, 16] <- sum(tot_res[[i]][,1] == TRUE & 
                                  tot_res[[i]][,2] == FALSE);
        all_res[i, 17] <- sum(tot_res[[i]][,1] == TRUE & 
                                  tot_res[[i]][,3] == FALSE);
        all_res[i, 18] <- sum(tot_res[[i]][,1] == TRUE & 
                                  tot_res[[i]][,4] == FALSE);
        all_res[i, 19] <- sum(tot_res[[i]][,1] == TRUE & 
                                  tot_res[[i]][,5] == FALSE);
        # Feasible and infeasible
        all_res[i, 20]  <- sum(fea_res[[i]][,1] == FALSE);
        all_res[i, 21]  <- sum(fea_res[[i]][,1] == TRUE);
        all_res[i, 22]  <- sum(fea_res[[i]][,2] == FALSE);
        all_res[i, 23]  <- sum(fea_res[[i]][,2] == TRUE);
        all_res[i, 24]  <- sum(fea_res[[i]][,3] == FALSE);
        all_res[i, 25]  <- sum(fea_res[[i]][,3] == TRUE);
        all_res[i, 26]  <- sum(fea_res[[i]][,4] == FALSE);
        all_res[i, 27]  <- sum(fea_res[[i]][,4] == TRUE);
        all_res[i, 28]  <- sum(fea_res[[i]][,4] == FALSE);
        all_res[i, 29]  <- sum(fea_res[[i]][,4] == TRUE);
        # Feased and defeased
        all_res[i, 30] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,2] == TRUE);
        all_res[i, 31] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,3] == TRUE);
        all_res[i, 32] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,4] == TRUE);
        all_res[i, 33] <- sum(fea_res[[i]][,1] == FALSE & 
                                  fea_res[[i]][,5] == TRUE);
        all_res[i, 34] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,2] == FALSE);
        all_res[i, 35] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,3] == FALSE);
        all_res[i, 36] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,4] == FALSE);
        all_res[i, 37] <- sum(fea_res[[i]][,1] == TRUE & 
                                  fea_res[[i]][,5] == FALSE);  
    }
    cnames <- c("N", "A0_unstable", "A0_stable", "A1_unstable", "A1_stable", 
                "A2_unstable", "A2_stable", "A3_unstable", "A3_stable", 
                "A4_unstable", "A4_stable", "A1_stabilised", "A2_stabilised",
                "A3_stabilised", "A4_stabilised", "A1_destabilised",
                "A2_destabilised", "A3_destabilised", "A4_destabilised", 
                "A0_infeasible", "A0_feasible", "A1_infeasible", "A1_feasible", 
                "A2_infeasible", "A2_feasible", "A3_infeasible", "A3_feasible", 
                "A4_infeasible", "A4_feasible", "A1_made_feasible",
                "A2_made_feasible", "A3_made_feasible", "A4_made_feasible",
                "A1_made_infeasible",  "A2_made_infeasible", 
                "A3_made_infeasible", "A4_made_infeasible");
    colnames(all_res) <- cnames;
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


make_gammas <- function(nn = 10, distribution = 1, mn = 10, sdd = 1){
    if(distribution == 0){
        dat          <- rep(x = mn, times = nn);
    }
    if(distribution == 1){
        dat             <- runif(n = nn, min = 0, max = 2*mn);
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



rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, gamma_sd = 1, 
                         gamma_mn = 10, C = 1){
    tot_res <- NULL;
    fea_res <- NULL;
    for(i in 2:max_sp){
        iter           <- iters;
        tot_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 7);
        fea_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 7);
        while(iter > 0){
            r_vec    <- rnorm(n = i, mean = 0, sd = rmx);
            A0_dat   <- rnorm(n = i * i, mean = 0, sd = 0.4);
            A0       <- matrix(data = A0_dat, nrow = i, ncol = i);
            A0       <- species_interactions(mat = A0, type = int_type);
            C_dat    <- rbinom(n = i * i, size = 1, prob = C);
            C_mat    <- matrix(data = C_dat, nrow = i, ncol = i);
            A0       <- A0 * C_mat;
            diag(A0) <- -1;
            gam0     <- make_gammas(nn = i, distribution = 0, sdd = gamma_sd);
            gam1     <- make_gammas(nn = i, distribution = 1, sdd = gamma_sd);
            gam2     <- make_gammas(nn = i, distribution = 2, sdd = gamma_sd);
            gam3     <- make_gammas(nn = i, distribution = 3, sdd = gamma_sd);
            gam4     <- make_gammas(nn = i, distribution = 4, sdd = gamma_sd);
            A1       <- A0 * gam1;
            A2       <- A0 * gam2;
            A3       <- A0 * gam3;
            A4       <- A0 * gam4;
            A0       <- A0 * gam0;
            A0_stb   <- max(Re(eigen(A0)$values)) < 0;
            A1_stb   <- max(Re(eigen(A1)$values)) < 0;
            A2_stb   <- max(Re(eigen(A2)$values)) < 0;
            A3_stb   <- max(Re(eigen(A3)$values)) < 0;
            A4_stb   <- max(Re(eigen(A4)$values)) < 0;
            A0_fea   <- min(-1*solve(A0) %*% r_vec) > 0;
            A1_fea   <- min(-1*solve(A1) %*% r_vec) > 0;
            A2_fea   <- min(-1*solve(A2) %*% r_vec) > 0;
            A3_fea   <- min(-1*solve(A3) %*% r_vec) > 0;
            A4_fea   <- min(-1*solve(A4) %*% r_vec) > 0;
            if(A0_stb == TRUE){
                tot_res[[i-1]][iter, 1] <- 1;
            }
            if(A1_stb == TRUE){
                tot_res[[i-1]][iter, 2] <- 1;
            }
            if(A2_stb == TRUE){
                tot_res[[i-1]][iter, 3] <- 1;
            }
            if(A3_stb == TRUE){
                tot_res[[i-1]][iter, 4] <- 1;
            }
            if(A4_stb == TRUE){
                tot_res[[i-1]][iter, 5] <- 1;
            }
            if(A0_fea == TRUE){
                fea_res[[i-1]][iter, 1] <- 1;
            }
            if(A1_fea == TRUE){
                fea_res[[i-1]][iter, 2] <- 1;
            }
            if(A2_fea == TRUE){
                fea_res[[i-1]][iter, 3] <- 1;
            }
            if(A3_fea == TRUE){
                fea_res[[i-1]][iter, 4] <- 1;
            }
            if(A4_fea == TRUE){
                fea_res[[i-1]][iter, 5] <- 1;
            }
            iter    <- iter - 1;
        }
        print(i);
    }
    all_res <- summarise_randmat(tot_res = tot_res, fea_res = fea_res);
    return(all_res);
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


################################################################################
################################################################################
################################################################################


image(-1*egi$hmat, col = heat.colors(2000));

make_e_image <- function(ei_ret){
     hmat   <- ei_ret$hmat;
     xcoord <- ei_ret$xcoord;
     ycoord <- ei_ret$ycoord;
     adj_ei <- -1 * log(hmat + 0.001);
     xaxis  <- which(xcoord == 0) / length(xcoord);
     yaxis  <- which(ycoord == 0) / length(ycoord);
     coomat <- matrix(data = 0, nrow = dim(hmat)[1], ncol = dim(hmat)[2]);
     coomat[,xaxis] <- 1;
     coomat[yaxis,] <- 1;
     image(adj_ei, col = heat.colors(2000), axes = FALSE);
     abline(h = xaxis, lwd = 2);
     abline(v = yaxis, lwd = 2);
     box();
}


eigen_image  <- function(eigen_list){
    reigs <- round(eigen_list, digits = 1);
    xcoord <- seq(from = -40, to = 40, by = 0.1);
    ycoord <- seq(from = -40, to = 40, by = 0.1);
    hmat   <- matrix(data = 0, ncol = length(xcoord), nrow = length(ycoord));
    for(xx in 1:length(xcoord)){
        for(yy in 1:length(ycoord)){
            heat           <- reigs[,1] == xcoord[xx] & reigs[,2] == ycoord[yy];
            hmat[yy, xx]   <- sum(heat);
        }
        if(xx %% 100 == 0){
            pct <- round((xx / length(xcoord) * 100), digits = 2);
            print(paste(pct, "percent complete"));
        }
    }
    colnames(hmat) <- as.character(xcoord);
    rownames(hmat) <- as.character(ycoord);
    return(list(xcoord = xcoord, ycoord = ycoord, hmat = hmat));
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
        co1      <- rbind(co1, e_coord0);
        co2      <- rbind(co2, e_coord0);
        co3      <- rbind(co3, e_coord0);
        co4      <- rbind(co4, e_coord0);
        iter     <- iter - 1;
    }
    return(list(d0 = co0, d1 = co1, d2 = co2, d3 = co3, d4 = co4));
} 


################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

# THIS IS WHERE THE RESULTS IN evolve-to-stability.csv WERE GENERATED

# sim <- rand_gen_var(max_sp = 50, iters = 1000)

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

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################

evolve_gamma <- function(nn = 24, iters = 1000){
    iter     <- iters;
    ress     <- NULL;
    countt   <- 1;
    while(iter > 0){
        A1_dat   <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
        A1       <- matrix(data = A1_dat, nrow = nn, ncol = nn);
        A1       <- species_interactions(mat = A1, type = 0);
        diag(A1) <- -1;
        is_stab  <- max(Re(eigen(A1)$values)) < 0;
        if(is_stab == FALSE){
            tmat <- rand_mat_ga(A1, max_it = 100, converg = 0.00001);
            if(tmat$find_st == TRUE & dim(tmat$indist)[1] > 4){
                ress[[countt]] <- tmat$indist;
                countt         <- countt + 1;
            }
        }
        iter <- iter - 1;
        print(c(iter, is_stab));
    }
    return(ress);
}

rand_mat_ga <- function(A1, max_it = 20, converg = 0.01){
    nn       <- dim(A1)[1];
    rind     <- runif(n = nn*1000, min = 0, max = 1);
    inds     <- matrix(data = rind, nrow = 1000, ncol = nn);
    lastf    <- -10;
    ccrit    <- 10;
    find_st  <- 0;
    iter     <- max_it;
    indist   <- NULL;
    while(iter > 0 & find_st < 1 & ccrit > converg){
        ivar  <- rep(x = 0, length = dim(inds)[1]);
        ifit  <- rep(x = 0, length = dim(inds)[1]);
        isst  <- rep(x = 0, length = dim(inds)[1]);
        for(i in 1:dim(inds)[1]){
            ifit[i] <- -1*max(Re(eigen(inds[i,]*A1)$values));
            ivar[i] <- var(inds[i,]);
            isst[i] <- max(Re(eigen(inds[i,]*A1)$values)) < 0;
        }
        mn       <- apply(X = inds, MAR = 1, FUN = mean);
        sdd      <- apply(X = inds, MAR = 1, FUN = sd);
        sk       <- apply(X = inds, MAR = 1, FUN = skew);
        kt       <- apply(X = inds, MAR = 1, FUN = kurt);
        indist   <- rbind(indist, c(mn, sdd, sk, kt));
        most_fit <- order(ifit, decreasing = TRUE)[1:20];
        parents  <- inds[most_fit,];
        new_gen  <- matrix(data = t(parents), nrow = 1000, ncol = nn, 
                           byrow = TRUE);
        mu_dat   <- rbinom(n = nn*1000, size = 1, prob = 0.3);
        mu_dat2  <- rnorm(n = nn*1000, mean = 0, sd = 0.02);
        mu_dat2[mu_dat2 < 0] <- -mu_dat2[mu_dat2 < 0];
        mu_dat2[mu_dat2 > 1] <- 1 - mu_dat2[mu_dat2 > 1] %% 1;
        mu_dat3  <- mu_dat * mu_dat2;
        mu_mat   <- matrix(data = mu_dat3, nrow = 1000, ncol = nn);
        new_gen  <- new_gen + mu_mat;
        new_gen  <- crossover(inds = new_gen, pr = 0.2);
        inds     <- new_gen;
        find_st  <- min(isst);
        newf     <- mean(ifit);
        ccrit    <- newf - lastf;
        lastf    <- newf;
        iter     <- iter - 1;
    }
    return(list(indist = indist, find_st = find_st));
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

skew <- function(vals){
    dat <- (vals - mean(vals)) / sd(vals);
    dat <- dat * dat * dat;
    mea <- mean(dat);
    return(mea);
}

kurt <- function(vals){
    dat <- (vals - mean(vals)) / sd(vals);
    dat <- dat * dat * dat * dat;
    mea <- mean(dat);
    return(mea);
}

   
    
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













#===============================================================================







#===============================================================================


rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, pr_dev = 1){
    tot_res <- NULL;
    fea_res <- NULL;
    nea_res <- NULL;
    for(i in 1:pr_dev){
        nn             <- max_sp;
        A1_stt         <- 0;
        A2_stt         <- 0;
        iter           <- iters;
        tot_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 2);
        fea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 2);
        nea_res[[i]]   <- matrix(data = 0, nrow = iter, ncol = 2);
        while(iter > 0){
            r_vec    <- rnorm(n = nn, mean = 0, sd = rmx);
            A1_dat   <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
            A1       <- matrix(data = A1_dat, nrow = nn, ncol = nn);
            A1       <- species_interactions(mat = A1, type = int_type);
            diag(A1) <- -1;
            epsil    <- rnorm(n = nn, mean = 0, sd = pr_dev+100);
            eps_mat  <- matrix(data = exp(epsil), nrow = nn, ncol = nn, 
                               byrow = FALSE);
            avg_mat  <- matrix(data = mean(exp(epsil)), nrow = nn, ncol = nn, 
                               byrow = FALSE);
            A2       <- A1 * eps_mat;
            A1       <- A1 * avg_mat;
            A1_stb   <- max(Re(eigen(A1)$values)) < 0;
            A2_stb   <- max(Re(eigen(A2)$values)) < 0;
            if(A1_stb == TRUE){
                tot_res[[i]][iter, 1] <- 1;
            }
            if(A2_stb == TRUE){
                tot_res[[i]][iter, 2] <- 1;
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
        all_res[i, 1]  <- i;
        all_res[i, 2]  <- sum(unstables);
        all_res[i, 3]  <- sum(stables);
        all_res[i, 4]  <- sum(unstabled);
        all_res[i, 5]  <- sum(stabled);
        all_res[i, 6]  <- sum(non_feas);
        all_res[i, 7]  <- sum(feasibl);
        all_res[i, 8]  <- sum(unfeased);
        all_res[i, 9]  <- sum(feased);
    }
    return(all_res);
}







###############################################################



rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4){
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
      
      SA1      <- sign(A1);
      Sadj     <- apply(X = -SA1, MARGIN = 1, FUN = sum);
      mnadj    <- 1 - min(Sadj);
      Sadj     <- exp(Sadj + mnadj);

      eps_mat  <- matrix(data = Sadj, nrow = nn, ncol = nn, 
                         byrow = FALSE);
      avg_mat  <- matrix(data = mean(Sadj), nrow = nn, ncol = nn, 
                         byrow = FALSE);
      A2       <- A1 * eps_mat;
      A1       <- A1 * avg_mat;
      A1_stb   <- max(Re(eigen(A1)$values)) < 0;
      A2_stb   <- max(Re(eigen(A2)$values)) < 0;
      #A1_fea   <- min(-1*solve(A1) %*% r_vec) > 0;
      #A2_fea   <- min(-1*solve(A2) %*% r_vec) > 0;
      if(A1_stb == TRUE){
        tot_res[[i-1]][iter, 1] <- 1;
      }
      if(A2_stb == TRUE){
        tot_res[[i-1]][iter, 2] <- 1;
      }
      #if(A1_fea == TRUE){
      #  fea_res[[i-1]][iter, 1] <- 1;
      #}
      #if(A2_fea == TRUE){
      #  fea_res[[i-1]][iter, 2] <- 1;
      #}
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
    all_res[i, 1]  <- i + 1;
    all_res[i, 2]  <- sum(unstables);
    all_res[i, 3]  <- sum(stables);
    all_res[i, 4]  <- sum(unstabled);
    all_res[i, 5]  <- sum(stabled);
    all_res[i, 6]  <- sum(non_feas);
    all_res[i, 7]  <- sum(feasibl);
    all_res[i, 8]  <- sum(unfeased);
    all_res[i, 9]  <- sum(feased);
  }
  return(all_res);
}









###############################################################



SA1      <- sign(A1);
Sadj     <- apply(X = -SA1, MARGIN = 1, FUN = sum);
mnadj    <- 1 - min(Sadj);
Sadj     <- exp(Sadj + mnadj);




rand_cov_res <- function(nn, iters, int_type = 0, rmx = 0.4, eps_max = 999){
    iter           <- iters;
    cols           <- 4 + 2*nn;
    res_mat        <- matrix(data = 0, nrow = iters, ncol = cols);
    while(iter > 0){
        A1_dat   <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
        A1       <- matrix(data = A1_dat, nrow = nn, ncol = nn);
        A1       <- species_interactions(mat = A1, type = int_type);
        diag(A1) <- -1;
        epsil    <- runif(n = nn, min = 1, max = eps_max);
        eps_dat  <- rep(x = epsil, times = nn);
        eps_mat  <- matrix(data = epsil, nrow = nn, ncol = nn, 
                           byrow = FALSE);
        A2       <- A1 * eps_mat;
        ceig_res <- cor_eigens(A2, eps_mat);
        A2_real  <- Re(eigen(A2)$values);
        A2_imag  <- Im(eigen(A2)$values);
        res_mat[iter, 1:4]               <- ceig_res;
        res_mat[iter, 5:(nn+4)]          <- A2_real;
        res_mat[iter, (nn+5):(2*nn + 4)] <- A2_imag;
        iter    <- iter - 1;
    }
    return(res_mat);
}

cor_eigens <- function(A2, eps_mat){
    epses    <- eps_mat[,1];
    SA2      <- sign(A2);
    Sadj     <- apply(X = -SA2, MARGIN = 1, FUN = sum);
    var_ep   <- var(epses);
    var_Sadj <- var(Sadj);
    cov_ep   <- cov(Sadj, epses);
    cor_ep   <- cor(Sadj, epses);
    cor_res  <- c(var_ep, var_Sadj, cov_ep, cor_ep);
    return(cor_res);
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





min_eig <- apply(rand_vals[,5:20], 1, max);
is_stab <- min_eig < 0;





mnadj    <- 1 - min(Sadj);
Sadj     <- Sadj + mnadj;




Re(eigen(exp(Sadj) * A2)$values);



tmat  <- matrix(data = 0, nrow = 10000, ncol = 25);

for(i in 1:dim(tmat)[1]){
    r_vec <- runif(n = 8, min = 1, max = 999);
    stbb  <- max(Re(eigen(r_vec * A2)$values)) < 0;
    tmat[i, 1] <- stbb;
    tmat[i, 2:9] <- r_vec;
    tmat[i, 10:17] <- Re(eigen(r_vec*A2)$values);
    tmat[i, 18:25] <- Im(eigen(r_vec*A2)$values);
}

x0 <- min(tmat[,10:17]);
x1 <- max(tmat[,10:17]);
y0 <- min(tmat[,18:25]);
y1 <- max(tmat[,18:25]);



plot(x = Re(eigen(A2)$values), y = Im(eigen(A2)$values), xlim = c(x0, x1),
     ylim = c(y0, y1));




plot(x = Re(eigen(500*A2)$values), y = Im(eigen(500*A2)$values), xlim = c(-1000000, 100000),
     ylim = c(-500000, 500000), cex = 0.6);
for(i in 1:1000){#dim(tmat)[1]){
    if(tmat[i,1] == 1){
        points(x = tmat[i, 10:17], y = tmat[i, 18:25], col = "red", cex = 0.6);
    }
    if(tmat[i,1] == 0){
        points(x = tmat[i, 10:17], y = tmat[i, 18:25], col = "blue", cex = 0.6);
    }
}
points(x = Re(eigen(500*A2)$values), y = Im(eigen(500*A2)$values), cex = 2, pch = 20);
abline(v = 0, lwd = 2, col = "orange");


stables <- tmat[tmat[,1] == 1,];
unstabl <- tmat[tmat[,1] == 0,];

hist(unstabl[,2:9],  breaks = 1000)


mn_r <- apply(tmat[,2:9],1,max);
mn_e <- apply(tmat[,10:17],1,max);




plot(x = Re(eigen(500*A2)$values), y = Im(eigen(500*A2)$values), xlim = c(-1000000, 10000),
     ylim = c(-2000000, 2000000), cex = 0.6);

        points(x = tmat[17, 10:17], y = tmat[17, 18:25], col = "red", cex = 0.6);

        points(x = tmat[16, 10:17], y = tmat[16, 18:25], col = "blue", cex = 0.6);

points(x = Re(eigen(500*A2)$values), y = Im(eigen(500*A2)$values),cex = 0.6);
abline(v = 0, lwd = 2, col = "orange");






##############################################################################

A1_dat   <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
A1       <- matrix(data = A1_dat, nrow = nn, ncol = nn);
A1       <- species_interactions(mat = A1, type = int_type);
diag(A1) <- -1;
Re(eigen(A1)$values)



tmat  <- matrix(data = 0, nrow = 1000, ncol = 25);

for(i in 1:dim(tmat)[1]){
    r_vec <- runif(n = 8, min = 0, max = 1);
    stbb  <- max(Re(eigen(r_vec * A1)$values)) < 0;
    tmat[i, 1] <- stbb;
    tmat[i, 2:9] <- r_vec;
    tmat[i, 10:17] <- Re(eigen(r_vec*A1)$values);
    tmat[i, 18:25] <- Im(eigen(r_vec*A1)$values);
}


x0 <- min(tmat[,10:17]);
x1 <- max(tmat[,10:17]);
y0 <- min(tmat[,18:25]);
y1 <- max(tmat[,18:25]);

plot(x = Re(eigen(A1)$values), y = Im(eigen(A1)$values), xlim = c(x0, x1), 
     ylim = c(y0, y1), cex = 0.6);
for(i in 1:1000){#dim(tmat)[1]){
    if(tmat[i,1] == 1){
        points(x = tmat[i, 10:17], y = tmat[i, 18:25], col = "red", cex = 0.6);
    }
    if(tmat[i,1] == 0){
        points(x = tmat[i, 10:17], y = tmat[i, 18:25], col = "blue", cex = 0.6);
    }
}
points(x = Re(eigen(500*A2)$values), y = Im(eigen(500*A2)$values), cex = 2, pch = 20);
abline(v = 0, lwd = 2, col = "orange");

stables <- tmat[tmat[,1] == 1,];
unstabl <- tmat[tmat[,1] == 0,];

mn_stabl <- apply(stables[,2:9],1,min);
mn_unstb <- apply(unstabl[,2:9],1,min);

mx_stabl <- apply(stables[,2:9],1,max);
mx_unstb <- apply(unstabl[,2:9],1,max);

av_stabl <- apply(stables[,2:9],1,mean);
av_unstb <- apply(unstabl[,2:9],1,mean);

vr_stabl <- apply(stables[,2:9],1,var);
vr_unstb <- apply(unstabl[,2:9],1,var);



mx_eig <- apply(tmat[,10:17],1,max);
vr_all <- apply(tmat[,2:9],1,var);

mn_all <- apply(tmat[,2:9],1,min);


plot(mn_all,mx_eig, xlab = "Minimum component response rate", ylab = "Maximum real eigenvalue", cex.lab = 1.5, cex.axis = 1.5, cex = 0.6)
abline(h = 0, lwd = 2, lty = "dotted", col = "red")
model <- lm(mx_eig ~ mn_all + I(mn_all^2) + I(mn_all^3))
xx <- seq(from = 0, to =0.6, length = 1000)
yy <- -0.013983 + -0.195314*xx + 1.357474 * xx^2 + -1.878256 * xx^3
points(xx, yy, type = "l", col = "blue")


av_all <- apply(tmat[,2:9],1,mean);

plot(av_all,mx_eig, xlab = "Mean component response rate", ylab = "Maximum real eigenvalue", cex.lab = 1.5, cex.axis = 1.5, cex = 0.6)
abline(h = 0, lwd = 2, lty = "dotted", col = "red")
model <- lm(mx_eig ~ av_all)
xx <- seq(from = 0, to =0.6, length = 1000)
yy <- -0.013983 + -0.195314*xx + 1.357474 * xx^2 + -1.878256 * xx^3
points(xx, yy, type = "l", col = "blue")






################################################################################
# GENETIC ALGORITHM:

nn       <- 24;
A1_dat   <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
A1       <- matrix(data = A1_dat, nrow = nn, ncol = nn);
A1       <- species_interactions(mat = A1, type = int_type);
diag(A1) <- -1;
t1       <- rand_mat_ga(A1, max_it = 100, converg = 0.00001);
hist(t1[1,], xlim = c(0, 1), breaks = 20);

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
        mu_dat   <- rbinom(n = nn*1000, size = 1, prob = 0.3);
        mu_dat2  <- rnorm(n = nn*1000, mean = 0, sd = 0.02);
        mu_dat2[mu_dat2 < 0] <- -mu_dat2[mu_dat2 < 0];
        mu_dat2[mu_dat2 > 1] <- 1 - mu_dat2[mu_dat2 > 1] %% 1;
        mu_dat3  <- mu_dat * mu_dat2;
        mu_mat   <- matrix(data = mu_dat3, nrow = 1000, ncol = nn);
        new_gen  <- new_gen + mu_mat;
        new_gen  <- crossover(inds = new_gen, pr = 0.2);
        inds     <- new_gen;
        find_st  <- max(isst);
        newf     <- mean(ifit);
        ccrit    <- newf - lastf;
        #print(c(iter, newf, lastf, mean(ivar), find_st));
        lastf    <- newf;
        iter     <- iter - 1;
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
      A3_stb   <- rand_mat_ga(A1, max_it = 100, converg = 0.00001);
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
  all_res <- matrix(data = 0, nrow = sims, ncol = 17);
  for(i in 1:sims){
    all_res[i, 1]  <- i + 1;
    all_res[i, 2]  <- sum(tot_res[[i]][,1] == FALSE);
    all_res[i, 3]  <- sum(tot_res[[i]][,1] == TRUE);
    all_res[i, 4]  <- sum(tot_res[[i]][,2] == FALSE);
    all_res[i, 5]  <- sum(tot_res[[i]][,2] == TRUE);
    all_res[i, 6]  <- sum(tot_res[[i]][,3] == FALSE);
    all_res[i, 7]  <- sum(tot_res[[i]][,3] == TRUE);
    all_res[i, 8]  <- sum(tot_res[[i]][,4] == FALSE);
    all_res[i, 9]  <- sum(tot_res[[i]][,4] == TRUE);
    all_res[i, 10] <- sum(fea_res[[i]][,1] == FALSE);
    all_res[i, 11] <- sum(fea_res[[i]][,1] == TRUE);
    all_res[i, 12] <- sum(fea_res[[i]][,2] == FALSE);
    all_res[i, 13] <- sum(fea_res[[i]][,2] == TRUE);
    all_res[i, 14] <- sum(fea_res[[i]][,3] == FALSE);
    all_res[i, 15] <- sum(fea_res[[i]][,3] == TRUE);
    all_res[i, 16] <- sum(fea_res[[i]][,4] == FALSE);
    all_res[i, 17] <- sum(fea_res[[i]][,4] == TRUE);
  }
  return(all_res);
}


make_gammas <- function(nn = 10, distribution = 1, mn = 10, sdd = 1){
    if(distribution == 0){
        dat          <- rep(x = mn, times = nn);
    }
    if(distribution == 1){
        dat          <- rnorm(n = nn, mean = mn, sd = sdd);
        dat[dat < 0] <- -dat[dat < 0];
    }
    if(distribution == 2){
        dat          <- runif(n = nn, min = 0, max = 20);
        dat          <- sdd * (dat / sd(dat));
    }
    if(distribution == 3){
        dat          <- sdd * rexp(n = nn);
    }
    if(distribution == 4){
        dat          <- 1 - (sdd * rexp(n = nn));
        dat          <- dat - min(dat);
        dat          <- dat - mean(dat) + mn;
        dat[dat < 0] <- -dat[dat < 0];
    }
    return(dat);
}



rand_gen_var <- function(max_sp, iters, int_type = 0, rmx = 0.4, eps_max = 999){
    tot_res <- NULL;
    fea_res <- NULL;
    for(i in 2:max_sp){
        nn             <- i;
        A1_stt         <- 0;
        A2_stt         <- 0;
        A1_fet         <- 0;
        A2_fet         <- 0;
        iter           <- iters;
        tot_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 4);
        fea_res[[i-1]] <- matrix(data = 0, nrow = iter, ncol = 4);
        while(iter > 0){
            r_vec    <- rnorm(n = nn, mean = 0, sd = rmx);
            A1_dat   <- rnorm(n = nn * nn, mean = 0, sd = 0.4);
            A1       <- matrix(data = A1_dat, nrow = nn, ncol = nn);
            A1       <- species_interactions(mat = A1, type = int_type);
            diag(A1) <- -1;
            gam0     <- make_gammas(nn = i, distribution = 0);
            gam1     <- make_gammas(nn = i, distribution = 1);
            gam2     <- make_gammas(nn = i, distribution = 2);
            gam3     <- make_gammas(nn = i, distribution = 3);
            gam4     <- make_gammas(nn = i, distribution = 4);
            A2       <- A1 * gam1;
            A3       <- A1 * gam2;
            A4       <- A1 * gam3;
            A1       <- A1 * gam0;
            A1_stb   <- max(Re(eigen(A1)$values)) < 0;
            A2_stb   <- max(Re(eigen(A2)$values)) < 0;
            A3_stb   <- max(Re(eigen(A3)$values)) < 0;
            A4_stb   <- max(Re(eigen(A4)$values)) < 0;
            A1_fea   <- min(-1*solve(A1) %*% r_vec) > 0;
            A2_fea   <- min(-1*solve(A2) %*% r_vec) > 0;
            A3_fea   <- min(-1*solve(A3) %*% r_vec) > 0;
            A4_fea   <- min(-1*solve(A4) %*% r_vec) > 0;
            if(A1_stb == TRUE){
                tot_res[[i-1]][iter, 1] <- 1;
            }
            if(A2_stb == TRUE){
                tot_res[[i-1]][iter, 2] <- 1;
            }
            if(A3_stb == TRUE){
                tot_res[[i-1]][iter, 3] <- 1;
            }
            if(A4_stb == TRUE){
                tot_res[[i-1]][iter, 4] <- 1;
            }
            if(A1_fea == TRUE){
                fea_res[[i-1]][iter, 1] <- 1;
            }
            if(A2_fea == TRUE){
                fea_res[[i-1]][iter, 2] <- 1;
            }
            if(A3_fea == TRUE){
                fea_res[[i-1]][iter, 3] <- 1;
            }
            if(A4_fea == TRUE){
                fea_res[[i-1]][iter, 4] <- 1;
            }
            iter    <- iter - 1;
        }
        print(i);
    }
    all_res <- summarise_randmat(tot_res = tot_res, fea_res = fea_res);
    return(all_res);
}



################################################################################


for(i in 1:dim(tmat)[1]){
    r_vec <- runif(n = 8, min = 0, max = 1);
    stbb  <- max(Re(eigen(r_vec * A1)$values)) < 0;
    tmat[i, 1] <- stbb;
    tmat[i, 2:9] <- r_vec;
    tmat[i, 10:17] <- Re(eigen(r_vec*A1)$values);
    tmat[i, 18:25] <- Im(eigen(r_vec*A1)$values);
}




###############################################################################







plot(x = Re(eigen(500*A2)$values), y = Im(eigen(500*A2)$values), xlim = c(-1000, 1000),
     ylim = c(-5000, 5000), cex = 0.6);
for(i in 1:1000){#dim(tmat)[1]){
    if(tmat[i,1] == 1){
        points(x = tmat[i, 10:17], y = tmat[i, 18:25], col = "red", cex = 0.6);
    }
    if(tmat[i,1] == 0){
        points(x = tmat[i, 10:17], y = tmat[i, 18:25], col = "blue", cex = 0.6);
    }
}
points(x = Re(eigen(500*A2)$values), y = Im(eigen(500*A2)$values), cex = 0.6);
abline(v = 0, lwd = 2, col = "orange");









