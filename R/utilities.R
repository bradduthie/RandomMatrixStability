#' Summarise the correlation of matrix elements
#' 
#' Returns the correlation between off-diagonal elements A_{ij} and A_{ji}
#'
#'@return Correlation between matrix off-diagonal elements
#'@param mat The tot_res list output from the rand_gen_var function
#'@examples
#'eg_mat  <- matrix(dat = rnorm(n = 16), nrow = 4);
#'sum_rand <- mat_rho(mat = eg_mat);
#'@export
mat_rho <- function(mat){
    if(dim(mat)[1] != dim(mat)[2]){
        stop("Error: Not a square matrix");
    }
    S  <- dim(mat)[1];
    od <- 0.5 * (S * S - S);
    ut <- rep(x = NA, times = od);
    lt <- rep(x = NA, times = od);
    ct <- 1;
    for(i in 1:dim(mat)[1]){
        for(j in 1:dim(mat)[2]){
            if(i < j){
                ut[ct] <- mat[i, j];
                lt[ct] <- mat[j, i]; 
                ct     <- ct + 1;
            }
        }
    }
    rho <- cor(ut, lt);
    return(rho);
}

#' Build a matrix with a specified correlation structure
#' 
#' Builds a random network with a pre-specified correlation structure between 
#' off-diagonal elements A_{ij} and A_{ji}
#'
#'@return A matrix with a pre-specified correlation structure.
#'@param S The size of the network (number of components)
#'@param sigma The standard deviation of network interaction strengths
#'@param mn The mean interaction strength
#'@param rho The correlation between component interaction strengths (i.e.,
#'between off-diagonal matrix elements A_{ij} and A_{ji}
#'@param dval Self-regulation of network elements (1 by default)
#'@examples
#'eg_mat <- build_rho_mat(S = 8, sigma = 0.4, rho = 0.2);
#'@export
build_rho_mat <- function(S, sigma, rho, mn = 0, dval = 1){
    mat  <- matrix(data = 0, nrow = S, ncol = S);
    ia   <- 0.5 * ((S * S) - S);
    ut   <- rnorm(n = ia, mean = mn, sd = sigma);
    lt   <- rnorm(n = ia, mean = mn, sd = sigma);
    od   <- cbind(ut, lt);
    co   <- var(od);
    ch   <- solve(chol(co));
    nx   <- od %*% ch;
    ms   <- matrix(data = c(1, rho, rho, 1), ncol = 2);
    c2   <- chol(ms);
    el   <- nx %*% c2 * sd(ut) + mean(ut);
    fc   <- 1;
    for(i in 1:S){
        for(j in 1:S){
            if(i < j){
                mat[i, j] <- el[fc, 1];
                mat[j, i] <- el[fc, 2];
                fc        <- fc + 1;
            }
        }
    }
    diag(mat) <- -1 * dval;
    return(mat);
}

#' Get the complexity of matrix
#' 
#' Returns the complexity of a complex system, defined by the number of system
#' components (S), connectance between components (C), and standard deviation of
#' interaction strengths (sigma) such that complexity equals sigma * sqrt(SC)
#'
#'@return The complexity of the matrix mat
#'@param mat The matrix to be valuated
#'@examples
#'eg_mat  <- matrix(dat = rnorm(n = 16), nrow = 4);
#'mat_cmp <- get_complexity(mat = eg_mat);
#'@export
get_complexity <- function(mat){
    S         <- length(diag(mat));
    mat_gz    <- sum(mat != 0);
    trace_gz  <- sum(diag(mat) != 0);
    offdiagC  <- mat_gz - trace_gz;
    offdiags  <- (dim(mat)[1] * dim(mat)[2]) - S;
    calc_C    <- offdiagC / offdiags;
    diag(mat) <- 0;
    sigma     <- sd(mat[mat != 0], na.rm = TRUE);
    complx    <- sigma * sqrt(S * calc_C);
    return(complx);
}

#' Get the connectance of a network
#' 
#' Returns the connectance of a network (square matrix), defined as the 
#' proportion of non-zero off-diagonal elements.
#'
#'@return The complexity of the matrix mat
#'@param mat The matrix to be valuated
#'@examples
#'eg_mat  <- matrix(dat = rnorm(n = 16), nrow = 4);
#'mat_cmp <- get_C(mat = eg_mat);
#'@export
get_C <- function(mat){
    c1  <- sum(mat != 0);
    tot <- (dim(mat)[1] * dim(mat)[2]) - length(diag(mat));
    return(c1 / tot);
}

#' Visualise a network
#' 
#' Returns a plot of a network (square matrix)
#'
#'@return The complexity of the matrix mat
#'@param mat The matrix to be valuated
#'@examples
#'eg_mat  <- matrix(dat = rnorm(n = 100), nrow = 10);
#'mat_cmp <- visualise_network(mat = eg_mat);
#'@export
visualise_network <- function(mat){
    N  <- dim(mat)[1];
    rd <- seq(from = 0, to = 2*pi, length = N + 1);
    yy <- sin(rd);
    xx <- cos(rd);
    par(bty = "n");
    plot(x = xx, y = yy, pch = 20, cex = 1.5, 
         xaxt = "n", yaxt = "n", xlab = "", ylab = "");
    for(i in 1:N){
        for(j in 1:N){
            if(i > j & mat[i,j] > 0){
                lines(x = c(xx[i], xx[j]), y = c(yy[i], yy[j]), 
                      lwd = 0.5);
            }
        }
    }
}

add_C_stats <- function(sim){
    real_Cs <- sim$real_Cs;
    all_res <- sim$all_res;
    res_mns <- matrix(data = 0, nrow = length(sim$real_Cs), ncol = 2);
    for(i in 1:dim(res_mns)[1]){
        mn_vals       <- apply(X = real_Cs[[i]], MARGIN = 2, FUN = mean);
        res_mns[i, 1] <- mn_vals[1];
        res_mns[i, 2] <- mn_vals[3];
    }
    new_all_res <- cbind(all_res, res_mns[,2]);
    colnames(new_all_res)[dim(new_all_res)[2]] <- "C";
    return(new_all_res);
}

#' Get variance of interaction strengths
#' 
#' Returns the variance of interaction stengths (sigma^2).
#'
#'@return The variance of off-diagonal elements of 'mat'
#'@param mat The matrix to be valuated
#'@examples
#'eg_mat  <- matrix(dat = rnorm(n = 16), nrow = 4);
#'mat_cmp <- get_sigma_sqd(mat = eg_mat);
#'@export
get_sigma_sqd <- function(mat){
    diag(mat) <- 0;
    mat_offs  <- mat[mat != 0];
    var_offs  <- var(mat_offs);
    return(var_offs);
}


simpleboot <- function(freqs , repli = 1000, alpha = 0.05){    
    vals  <- NULL;                                                          
    i     <- 0;                                         
    while(i < repli){                                                   
        boot  <- sample(x = freqs, size = length(freqs), replace = TRUE);           
        strap <- mean(boot);                                                
        vals  <- c(vals,strap);                                             
        i     <- i + 1;                                                     
    }                                                                       
    vals   <- sort(x = vals, decreasing = FALSE);                            
    lowCI  <- vals[round((alpha*0.5)*repli)];                               
    highCI <- vals[round((1-(alpha*0.5))*repli)];                           
    CIs    <- c(lowCI,highCI);                                          
    return(CIs);                                                            
}  
