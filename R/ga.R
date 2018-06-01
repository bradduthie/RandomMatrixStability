#' Run genetic algorithms to find stability
#' 
#' Outer function for investigating the potential for a complex system to be
#' stabilised by component response rates. For a stretch of values of component
#' number (i.e., system size, S), the function runs some number of iterations, 
#' and in each iteration, if a system is not initially found to be stable with
#' a set of component response times, then a genetic algorithm is called in 
#' attempt to find a stable solution.
#'
#'@return A table showing frequency of stable complex systems found
#'@param max_sp The maximum system component size (S) to simulate
#'@param iters The number of iterations (complex systems) tried for each S
#'@param int_type The type of interaction, including random (0), competitive 
#'(1), mutualist (2), or predator-prey (3)
#'@param rmx The standard devation of interaction (off-diagonal) elements in a
#'complex system
#'@param C The connectedness of a complex system (i.e., the probability that an
#'off-diagonal element does not equal zero)
#'@examples
#'sim <- Evo_rand_gen_var(max_sp = 2, iters = 1);
#'@export
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
            gam1     <- runif(n = i, min = 0, max = 2);
            A1       <- A0 * gam1;
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


#' Genetic algorithm to find stability
#' 
#' Runs a genetic algorithm in attempt to find a vector (gamma) that, when the 
#' elements of which are multiplied as scalar values on a matrix, A1, which 
#' represents a complex system, causes all real parts of eigenvalues of A1 to
#' be negative, and the complex system to therefore be stable.
#'
#'@return A value indicating whether or not a stabilising vector is (1) or is
#'not (0) found by the genetic algorithm. If so, then the vector is appended to
#'the file "evo_out.txt".
#'@param A1 Square matrix representing a complex system
#'@param max_it Maximum number of generations the genetic algorithm is allowed
#'to run before giving up and declaring no stabilising set of gammas is found
#'@param converg Convergence criteria for the genetic algorithm; if the mean 
#'fitness from one generation to the next (in terms of improvement shifting
#'leading Real parts of  eigenvalues to the left) is below this, then the 
#'algorithm gives up and declares no stabilising set of gammas is found
#'@examples
#'A_mat <- matrix(data = rnorm(n = 16), nrow = 4);
#'sim   <- rand_mat_ga(A1 = A_mat, max_it = 10);
#'@export
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
        mu_dat2[mu_dat2 > 2] <- 2;
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
        cat(writt, file = "evo_out.txt", append = TRUE);
        cat("\n", file = "evo_out.txt", append = TRUE);
    }
    return(find_st);
}

# This function is used within the genetic algorithm
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

#' Summarise random matrix results from genetic algorithm
#' 
#' Takes the list output of the Evo_rand_gen_var function and summarises it into
#' a useable format.
#'
#'@return A table of stability and feasibility results.
#'@param tot_res The tot_res list output from the rand_gen_var function
#'@param fea_res The fea_res list output from the rand_gen_var function
#'@examples
#'sim_ga      <- Evo_rand_gen_var(from = 2, to = 2, iters = 1);
#'sum_rand_ga <- summarise_randmat_ga(sim_ga$tot_res, sim_ga$fea_res);
#'@export
summarise_randmat_ga <- function(tot_res, fea_res){
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


#' Find the component response rate vectors of highest size
#' 
#' Takes the output of the genetic algorithm showing the component response rate
#' (gamma) vectors found to be stabilising and returns the ones that are from 
#' the largest systems (i.e., highest S)
#'
#'@return A list of vectors of component response rates (gammas)
#'@param evo_out Output of component response values
#'@param size The number of gamma vectors to be returned
#'@export
get_top_evo_out <- function(evo_out, size){
    highest <- max(evo_out);
    gammas  <- NULL;
    while(size > 0){
        pos <- which(evo_out == highest);
        len <- length(pos);
        nli <- NULL;
        for(i in 1:len){
            start          <- pos[i] + 1;
            end            <- pos[i] + highest;
            gammas[[size]] <- evo_out[start:end];
            size           <- size - 1;
            if(size == 0){
                break;
            }
        }
        highest <- highest - 1;
    }
    return(gammas);
}

#' Plot component response rate distributions
#' 
#' Plots the component response rate distributions from the largest systems that
#' the genetic algorithm finds to be stable. At the moment, only a six panel
#' plot is allowed
#'
#'@return A list of vectors of component response rates (gammas)
#'@param evo_out Output of component response values
#'@export
plot_evo_out <- function(evo_out){
    gammas <- get_top_evo_out(evo_out, size = 9);
    par(mfrow = c(3, 3), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5, 5, 1, 1));
    hist(gammas[[1]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         xaxt = "n", xlim = c(0, 1.1));
    hist(gammas[[2]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         xaxt = "n", yaxt = "n", xlim = c(0, 1.1));
    hist(gammas[[3]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         xaxt = "n", yaxt = "n", xlim = c(0, 1.1));
    hist(gammas[[4]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         xaxt = "n", xlim = c(0, 1.1));
    hist(gammas[[5]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         xaxt = "n", yaxt = "n", xlim = c(0, 1.1));
    hist(gammas[[6]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         xaxt = "n", yaxt = "n", xlim = c(0, 1.1));
    hist(gammas[[7]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         xlim = c(0, 1.1));
    hist(gammas[[8]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         yaxt = "n", xlim = c(0, 1.1));
    hist(gammas[[9]], main = "", breaks = 20, col = "grey", ylim = c(0, 6),
         yaxt = "n", xlim = c(0, 1.1));
    mtext(side = 1, outer = TRUE, line = 3, cex = 1.5,
          text = expression(paste("Component response rate (",gamma,")")));
    mtext(side = 2, text = "Frequency", outer = TRUE, line = 3, cex = 1.5);
}







