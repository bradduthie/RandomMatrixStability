# NOTE FOR RESEARCH 27 OCT 2019
# SEE EQN 2.4-2.9 OF AHMADIAN ET AL 2015
# PARTICULARLY 2.8, IN WHICH THE SPECTRAL DENSITY IS DEPENDENT UP ON THE MATRIX
# L. THE DISTRIBUTION APPEARS TO DEPEND ON L^-1 IN A WAY (TRACE) THAT MIGHT 
# CAUSE THE STABILITY. NEED TO FIGURE OUT HOW THIS WORKS FOR N < 10 AND N > 10?
N  <- 4000;
z  <- 10;
M  <- matrix(data = rnorm(N * N, sd = 1/N), nrow = N, ncol = N);
Mz <- z - M;
ml <- Mz %*% t(Mz);
ms <- solve(ml);
tr <- (1/N) * sum(diag(ms));
tr;


Rg <- Re(eigen(M)$values);
Ig <- Im(eigen(M)$values);
plot(x = Rg, y = Ig);

# Maybe instead, define L as {N, N, N, N}, so the SD is preserved? Then
# use the tools above?

# Another check: Show unusually high values can be shunted into the very
# negative group.
N <- 800;
L <- matrix(data = 0, nrow = N, ncol = N);
J <- matrix(data = rnorm(N*N, sd = sqrt(1/N)), nrow = N, ncol = N);
diag(J)  <- 1;
diag(L)  <- -1 * c(rep(x = 0.01, times = N/2), rep(x = 1.98, times = N/2));
M1       <- L %*% J
eM1      <- eigen(M1)$values;
rM1      <- Re(eM1);
iM1      <- Im(eM1);
plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", xlim = c(-0.015, 0.015),
     ylim = c(-0.015, 0.015));
abline(v = 0, col = "red");
plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-1.5, -1.5), ylim = c(-1.5, 1.5));
rm1 <- tail(sort(rM1));

L        <- matrix(data = 0, nrow = N, ncol = N);
diag(L)  <- -0.995;
M2       <- L %*% J;
eM2      <- eigen(M2)$values;
rM2      <- Re(eM2);
iM2      <- Im(eM2);
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", xlim = c(-1.5, 1.5),
     ylim = c(-1.5, 1.5));
abline(v = 0, col = "red");
rm2 <- tail(sort(rM2));

rm1;
rm2;

plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "x", 
     xlim = c(-1.5, -1.5), ylim = c(-1.5, 1.5));
points(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "x", xlim = c(-1.5, 1.5),
     ylim = c(-1.5, 1.5), col = "blue");
abline(v = 0, col = "red");


plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "x", xlim = c(-0.15, 0.15),
     ylim = c(-0.15, 0.15));
points(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "x", col = "blue");
abline(v = 0, col = "red");

################################################################################
################################################################################
################################################################################


# Let's compare the eigenvalue distributions of the upper-left of M1, versus
# the whole eigenvalues to see what the rest adds.
M_ul  <- M1[1:400, 1:400];
em_ul <- eigen(M_ul)$values;
r_ul  <- Re(em_ul);
i_ul  <- Im(em_ul);

# Zoom in
plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02));
abline(v = 0, col = "red");
points(x = r_ul, y = i_ul, asp = 1, cex = 0.8, pch = "x", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02), col = "orange");
# The cloud is clearly bigger for the eigenvalues of the full matrix.

# Now add in the rest of it.
M_ll  <- M1[401:800, 401:800];
em_ll <- eigen(M_ll)$values;
r_ll  <- Re(em_ll);
i_ll  <- Im(em_ll);

plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-3.5, 0.5), ylim = c(-2, 2));
abline(v = 0, col = "red");
points(x = r_ul, y = i_ul, asp = 1, cex = 0.8, pch = "x", 
       xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02), col = "orange");
points(x = r_ll, y = i_ll, asp = 1, cex = 0.8, pch = "x", 
       xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02), col = "orange");
# The size of the bigger cloud does not appear to have changed much.
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sqrt(1/N) * 0.01 * sqrt(N/2) * cos(vl) + mean(diag(M_ul));
y0 <- sqrt(1/N) * 0.01 * sqrt(N/2) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3, col = "orange");

################################################################################
################################################################################
################################################################################
# SOME PROGRESS
################################################################################
################################################################################
################################################################################
# The following code will make a random complex system with the following
# parameter values:
# S       = 800
# sigma^2 = 0.1870829
# sigma   = 0.035
# C       = 1
# Make the random matrix with no variation in gamma
S        <- 800;
sigma    <- 0.036;
M        <- matrix(data = rnorm(S*S, sd = sigma), nrow = S, ncol = S);
diag(M)  <- -1;
M1       <- M;
L        <- matrix(data = 0, nrow = S, ncol = S);
diag(L)  <- rep(x = 1, times = S);
eM1      <- eigen(M1)$values;
rM1      <- Re(eM1);
iM1      <- Im(eM1);
plot(x = rM1, y = iM1, asp = 1, cex = 0.8, pch = "+", xlim = c(-2.5, 0.5),
     ylim = c(-1.5, 1.5), xlab = "Real", ylab = "Imaginary", cex.lab = 1.25);
abline(v = 0, col = "red");

# Now make the one with half slow components and half fast components
M2       <- M;
L        <- matrix(data = 0, nrow = S, ncol = S);
diag(L)  <- c(rep(x = 0.01, times = S/2), rep(x = 1.99, times = S/2));
diag(M2) <- -1;
M2       <- L %*% M2;
eM2      <- eigen(M2)$values;
rM2      <- Re(eM2);
iM2      <- Im(eM2);
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", xlim = c(-3.5, 0.5),
     ylim = c(-2, 2));
abline(v = 0, col = "red");
# Zoom in on that small circle
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02));
abline(v = 0, col = "red");

# See what the eigenvalue distributions would be if we only had two matrices
# Upper left, and lower right.
M3       <- M2[1:(S/2), 1:(S/2)];
eM3      <- eigen(M3)$values;
rM3      <- Re(eM3);
iM3      <- Im(eM3);
# Look at that circle
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02),
     xlab = "Real", ylab = "Imaginary", cex.lab = 1.25);
abline(v = 0, col = "red");
points(x = rM3, y = iM3, asp = 1, cex = 0.8, pch = "x", 
     xlim = c(-0.03, 0.01), ylim = c(-0.02, 0.02), col = "orange");
# This corresponds to the expected diameter
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sigma * 0.01 * sqrt(S/2) * cos(vl) + mean(diag(M3));
y0 <- sigma * 0.01 * sqrt(S/2) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3, col = "orange");
# Check out what happens when you change S
# Something is increasing the size of the circle.
# Solve this analytically.

# See what both look like if separate matrices, upper left and lower right.
M4       <- M2[(S/2 + 1):S, (S/2 + 1):S];
eM4      <- eigen(M4)$values;
rM4      <- Re(eM4);
iM4      <- Im(eM4);
# Look at that circle
plot(x = rM2, y = iM2, asp = 1, cex = 0.8, pch = "+", xlim = c(-3.5, 0.5),
     ylim = c(-2, 2), xlab = "Real", ylab = "Imaginary", cex.lab = 1.25);
abline(v = 0, col = "red");
# Add the old M3 ones
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x0 <- sigma * 0.01 * sqrt(S/2) * cos(vl) + mean(diag(M3));
y0 <- sigma * 0.01 * sqrt(S/2) * sin(vl);
points(x = x0, y = y0, type = "l", lwd = 3, col = "orange");
# Add the new M4 ones
vl <- seq(from = 0, to = 2*pi, by = 0.001);
x1 <- sigma * 1.99 * sqrt(S/2) * cos(vl) + mean(diag(M4));
y1 <- sigma * 1.99 * sqrt(S/2) * sin(vl);
points(x = x1, y = y1, type = "l", lwd = 3, col = "orange");

# Need to figure out where that variation is going.
