# Build a scale free network
find_sfree <- function(S){
    tverts   <- 0;
    while(tverts < S){
        mat      <- sfree(S = 2*S);
        verts    <- apply(X = mat, MARGIN = 1, FUN = sum);
        overts   <- which(verts > 0);
        cmat     <- mat[overts, overts];
        tverts   <- dim(cmat)[1];
    }    
    if(dim(cmat)[1] > S){
        nmat <- cmat[1:S, 1:S];
    }
    return(nmat);
}

# Need to seed with some m such that m can be a fixed node number
sfree <- function(S, m){
    mat <- matrix(data = 0, nrow = S, ncol = S);
    mat <- seed_mat(mat, m);
    for(i in m:S){
        new_p <- get_p(mat = mat[1:i, 1:i], m = m);
        mat   <- add_p(mat = mat, new_p = new_p);
    }
    return(mat);
}

seed_mat <- function(mat, m){
    mat[1:m, 1:m]  <- 1;
    diag(mat[1:m]) <- 0;
    return(mat);
}

# This function gets links for the new node
get_p <- function(mat, m){
    msize    <- dim(mat)[1];
    edges    <- apply(X = mat, MARGIN = 1, FUN = sum);
    prs      <- edges / sum(edges);
    ch       <- sample(x = 1:msize, size = m, prob = prs);
    newe     <- rep(x = 0, times = msize);
    newe[ch] <- 1;
    return(newe);
}

add_p <- function(mat, new_p){
    row             <- length(new_p);
    mat[row, 1:row] <- new_p;
    mat[1:row, row] <- new_p;
    return(mat);
}



