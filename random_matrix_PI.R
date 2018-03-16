



random_comm <- function(dm = 2, iter = 10000){
    stabilised_fastevo <- NULL;
    unstabled_fastevo  <- NULL;
    nn           <- dm * dm;
    tot_stabled  <- 0;
    tot_ustabled <- 0;
    tot_eco_stab <- 0;
    stabilised_fastevo <- NULL;
    unstabled_fastevo  <- NULL;
    while(iter > 0){
        stabilised   <- 0;
        unstabilised <- 0;
        Adat <- runif(n = nn, min = -4, max = 4);
        Bdat <- runif(n = nn, min = -4, max = 4);
        Cdat <- runif(n = nn, min = -4, max = 4);
        A    <- matrix(data = Adat, nrow = dm, ncol = dm);
        B    <- matrix(data = Bdat, nrow = dm, ncol = dm);
        C    <- matrix(data = Cdat, nrow = dm, ncol = dm);
        egD  <- 1;
        while(egD >= 0){
             Ddat <- runif(n = nn, min = -4, max = 4);
             D    <- matrix(data = Ddat, nrow = dm, ncol = dm);
             egD  <- max(Re(eigen(D)$values))
        }
        eco_stable       <- max(Re(eigen(A)$values)) < 0;
        slow_evo_maxE    <- max(Re(eigen(A + B %*% solve(D) %*% (-C))$values));
        stable_fast_evo  <- slow_evo_maxE < 0;
        if(eco_stable == FALSE & stable_fast_evo == TRUE){
            list_ele <- length(stabilised_fastevo$A) + 1;
            stabilised_fastevo$A[[list_ele]] <- A;
            stabilised_fastevo$B[[list_ele]] <- B;
            stabilised_fastevo$C[[list_ele]] <- C;
            stabilised_fastevo$D[[list_ele]] <- D;
            stabilised   <- 1;
        }
        if(eco_stable == TRUE & stable_fast_evo == FALSE){
            list_ele <- length(unstabled_fastevo$A) + 1;
            unstabled_fastevo$A[[list_ele]] <- A;
            unstabled_fastevo$B[[list_ele]] <- B;
            unstabled_fastevo$C[[list_ele]] <- C;
            unstabled_fastevo$D[[list_ele]] <- D;
            unstabilised   <- 1;
        }
        if(stabilised == 1){
            tot_stabled <- tot_stabled + 1;
        }
        if(unstabilised == 1){
            tot_ustabled <- tot_ustabled + 1;
        }
        if(eco_stable == TRUE){
            tot_eco_stab <- tot_eco_stab + 1;
        }
        iter <- iter - 1;
    }
    results <- list(ecologically_stable = tot_eco_stab, 
                    stabilised = tot_stabled, destabilised = tot_ustabled,
                    stabilised_all   = stabilised_fastevo, 
                    destabilised_all = unstabled_fastevo);
    return(results);
}


NN       <- 8;
res_tabl <- matrix(data = 0, nrow = (NN-1), ncol = 5);
its      <- 10000;
for(i in 2:NN){
    res              <- random_comm(dm = i, iter = its); 
    res_tabl[i-1, 1] <- i;                               # Species and traits
    res_tabl[i-1, 2] <- res$ecologically_stable;         # Total stable
    res_tabl[i-1, 3] <- res$destabilised;                # Total destabilised
    res_tabl[i-1, 4] <- its - res_tabl[i-1, 2];          # Total unstable
    res_tabl[i-1, 5] <- res$stabilised;                  # Total stabilised
    print(i);
}

# Big idea: More species or traits make it more likely that stable ecology
# will be de-stabilised by evolution, and less likely that unstable ecology
# will be stabilised by evolution.










tA <- NULL;
tB <- NULL;
tC <- NULL;
for(i in 1:10000){
  tA <- c(tA, diag(stabilised_fastevo$A[[i]]));
  tB <- c(tB, diag(stabilised_fastevo$B[[i]]));
  tC <- c(tC, diag(stabilised_fastevo$C[[i]]));
  
}





alphas <- as.vector(stabilised_fastevo$A[[3]])
betas  <- as.vector(stabilised_fastevo$B[[3]])
cesas  <- as.vector(stabilised_fastevo$C[[3]])
cor.test(alphas, betas);
























epsil   <- runif(n = dm, min = 1, max = 1000);
eps_dat <- rep(x = epsil, times = dm);
eps_mat <- matrix(data = c(epsil, epsil), nrow = dm, ncol = dm, byrow = FALSE);



random_comm <- function(dm = 2, iter = 10000){
  stabilised_fastevo <- NULL;
  unstabled_fastevo  <- NULL;
  nn           <- dm * dm;
  tot_stabled  <- 0;
  tot_ustabled <- 0;
  tot_eco_stab <- 0;
  stabilised_fastevo <- NULL;
  unstabled_fastevo  <- NULL;
  while(iter > 0){
    epsil   <- runif(n = dm, min = 1, max = 1000);
    eps_dat <- rep(x = epsil, times = dm);
    eps_mat <- matrix(data = c(epsil, epsil), nrow = dm, ncol = dm, byrow = FALSE);
    stabilised   <- 0;
    unstabilised <- 0;
    Adat <- runif(n = nn, min = -4, max = 4);
    Bdat <- runif(n = nn, min = -4, max = 4);
    Cdat <- runif(n = nn, min = -4, max = 4);
    A    <- matrix(data = Adat, nrow = dm, ncol = dm) * eps_mat;
    B    <- matrix(data = Bdat, nrow = dm, ncol = dm) * eps_mat;
    C    <- matrix(data = Cdat, nrow = dm, ncol = dm) * eps_mat;
    egD  <- 1;
    while(egD >= 0){
      Ddat <- runif(n = nn, min = -4, max = 4);
      D    <- matrix(data = Ddat, nrow = dm, ncol = dm) * eps_mat;
      egD  <- max(Re(eigen(D)$values))
    }
    eco_stable       <- max(Re(eigen(A)$values)) < 0;
    slow_evo_maxE    <- max(Re(eigen(A + B %*% solve(D) %*% (-C))$values));
    stable_fast_evo  <- slow_evo_maxE < 0;
    if(eco_stable == FALSE & stable_fast_evo == TRUE){
      list_ele <- length(stabilised_fastevo$A) + 1;
      stabilised_fastevo$A[[list_ele]] <- A;
      stabilised_fastevo$B[[list_ele]] <- B;
      stabilised_fastevo$C[[list_ele]] <- C;
      stabilised_fastevo$D[[list_ele]] <- D;
      stabilised   <- 1;
    }
    if(eco_stable == TRUE & stable_fast_evo == FALSE){
      list_ele <- length(unstabled_fastevo$A) + 1;
      unstabled_fastevo$A[[list_ele]] <- A;
      unstabled_fastevo$B[[list_ele]] <- B;
      unstabled_fastevo$C[[list_ele]] <- C;
      unstabled_fastevo$D[[list_ele]] <- D;
      unstabilised   <- 1;
    }
    if(stabilised == 1){
      tot_stabled <- tot_stabled + 1;
    }
    if(unstabilised == 1){
      tot_ustabled <- tot_ustabled + 1;
    }
    if(eco_stable == TRUE){
      tot_eco_stab <- tot_eco_stab + 1;
    }
    iter <- iter - 1;
    if(iter %% 1000000 == 0){
        print(c(dm, iter));
    }
  }
  results <- list(ecologically_stable = tot_eco_stab, 
                  stabilised = tot_stabled, destabilised = tot_ustabled,
                  stabilised_all   = stabilised_fastevo, 
                  destabilised_all = unstabled_fastevo);
  return(results);
}

NN       <- 20;
res_tabl <- matrix(data = 0, nrow = (NN-1), ncol = 5);
its      <- 100000000;
for(i in 2:NN){
  res              <- random_comm(dm = i, iter = its); 
  res_tabl[i-1, 1] <- i;                               # Species and traits
  res_tabl[i-1, 2] <- res$ecologically_stable;         # Total stable
  res_tabl[i-1, 3] <- res$destabilised;                # Total destabilised
  res_tabl[i-1, 4] <- its - res_tabl[i-1, 2];          # Total unstable
  res_tabl[i-1, 5] <- res$stabilised;                  # Total stabilised
  print(i);
}
