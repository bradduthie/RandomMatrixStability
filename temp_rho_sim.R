

sim30 <- rand_rho_var(S = 30, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);
write.csv(sim30, "notebook/sim_results/rhos/sigma0pt2/sim30.csv");

sim35 <- rand_rho_var(S = 35, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);
write.csv(sim35, "notebook/sim_results/rhos/sigma0pt2/sim35.csv");

sim40 <- rand_rho_var(S = 40, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);
write.csv(sim40, "notebook/sim_results/rhos/sigma0pt2/sim40.csv");

sim45 <- rand_rho_var(S = 45, rhos = seq(from = -0.9, to = 0.9, by = 0.05), 
                      iters = 10000, sigma = 0.2);
write.csv(sim45, "notebook/sim_results/rhos/sigma0pt2/sim45.csv");
