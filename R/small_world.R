

create_swn <- function(N = 100, K = 20, beta = 0.05){
    mat <- matrix(data = 0, nrow = N, ncol = N);
    diag(mat) <- 1;
    if(K %% 2 != 0){
        warning("K needs to be an even integer");
    }
    if(K > dim(mat)[1]){
        warning("K is too high");
    }
    mat <- add_sw_edges(mat, K);
    mat <- rewire_sw_edges(mat, K, beta);
    return(mat);
}

add_sw_edges <- function(mat, K){
    N <- dim(mat)[1];
    for(i in 1:N){
        for(j in 1:N){
            if(i > j){
                neig1 <- abs(i - j);
                neig2 <- abs(i - (N+j));
                neig  <- min(c(neig1, neig2));
                if(0 < neig & neig <= (K/2)){
                    mat[i, j] <- 1;
                    mat[j, i] <- 1;
                }
            }
       }
    }
    return(mat);
}

rewire_sw_edges <- function(mat, K, beta){
    N <- dim(mat)[1];
    for(i in 1:N){
        for(j in 1:N){
            if(i != j){
                samp_beta <- runif(n = 1) < beta;
                if(samp_beta){ # Rewire node
                    mat[i, j] <- 0;
                    mat[j, i] <- 0;
                    k         <- sample_k(i, N, K);
                    mat[i, k] <- 1;
                    mat[k, i] <- 1;
                }
            }
        }
    }
    return(mat);
}

sample_k <- function(i, N, K){
    tsamp <- 1:N;
    lo    <- i - (K/2);
    hi    <- i + (K/2);
    span1 <- seq(from = lo, to = hi, by = 1);
    span2 <- span1 + N;
    span  <- c(span1, span2);
    nouse <- span[span %in% tsamp];
    tsamp <- tsamp[-nouse];
    getit <- sample(tsamp, size = 1);
    return(getit);
}

