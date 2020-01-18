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

sfree <- function(S){
    mat <- matrix(data = 0, nrow = S, ncol = S);
    mat <- seed_mat(mat);
    for(i in 3:S){
        new_p <- get_p(mat = mat[1:i, 1:i]);
        mat   <- add_p(mat = mat, new_p = new_p);
    }
    return(mat);
}

seed_mat <- function(mat){
    mat[1, 2] <- 1;
    mat[2, 1] <- 1;
    return(mat);
}

# This function gets links for the new node
get_p <- function(mat){
    msize <- dim(mat)[1];
    edges <- apply(X = mat, MARGIN = 1, FUN = sum);
    prs   <- edges / sum(edges);
    newe  <- rbinom(n = msize, size = 1, prob = prs);
    return(newe);
}

add_p <- function(mat, new_p){
    row             <- length(new_p);
    mat[row, 1:row] <- new_p;
    mat[1:row, row] <- new_p;
    return(mat);
}



