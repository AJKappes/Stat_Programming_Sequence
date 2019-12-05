# STAT 536 Project
# Gibbs sampling for hierarchical normal model 
remove(list = objects())

#### conditional posterior distributions ------------------
theta_gen <- function(mu, n_j, y_j, tau2, sig2) {
  
  # conditional dist is normal
  theta_j <- 1/(1/tau2 + n_j/sig2)*(mu/tau2 + mean(y_j)*n_j/sig2)
  V_j <- 1/(1/tau2 + n_j/sig2)
  theta <- rnorm(1, theta_j, V_j)
  return(theta)
  
}

mu_gen <- function(theta_vec, tau2, G = J) {
  
  # conditional dist is normal
  mu_hat <- mean(theta_vec)
  s2 <- tau2/G
  mu <- rnorm(1, mu_hat, s2)
  return(mu)
  
}

sig2_gen <- function(theta_vec, dl = dlist, n = N) {
  
  # conditional dist is inv chisq
  sums <- c()
  for (i in 1:length(dl)) {
    
    f <- sum((dl[[i]] - theta_vec[i])^2)
    sums <- c(sums, f)
    
  }
  
  s2_hat <- 1/n*sum(sums)
  # x ~ inv-chisq(n, sig2) <=> x ~ invgamma(n/2, n*s2/2)
  s2 <- 1/rgamma(1, n/2, scale = n*s2_hat/2)
  return(s2)
  
}

tau2_gen <- function(theta_vec, mu, G = J) {
  
  # conditional dist is inv chisq
  t2_hat <- 1/(G - 1)*sum((theta_vec - mu)^2)
  # x ~ inv-chisq(n, sig2) <=> x ~ invgamma(n/2, n*s2/2)
  t2 <- 1/rgamma(1, (G - 1)/2, scale = (G - 1)*t2_hat/2)
  return(t2)
  
}

#### Data prep and initialization ------------------------
# diets A, B, C, D
# measurment is coagulation time in seconds for blood drawn from
#   24 random animals allocated to 4 different diets
A <- c(62, 60, 63, 59)
B <- c(63, 67, 71, 64, 65, 66)
C <- c(68, 66, 71, 67, 68, 68)
D <- c(56, 62, 60, 61, 63, 64, 63, 59)
dlist <- list(A, B, C, D)

J <- length(dlist)
N <- sum(sapply(dlist, length))

thetas_m <- sapply(dlist, mean)
mu_m <- mean(thetas_m)

nsim <- 1E+3
sim_params <- c('theta1', 'theta2', 'theta3', 'theta4',
                'mu', 'sig2', 'tau2')
sim_data <- setNames(data.frame(matrix(0, nrow = nsim, ncol = length(sim_params))),
                     sim_params)

#### gibbs sampling --------------------------------------
for (i in 1:nsim) {
  
  sig2 <- sig2_gen(thetas_m)
  tau2 <- tau2_gen(thetas_m, mu_m)
  mu_mp1 <- mu_gen(thetas_m, tau2)
  thetas_mp1 <- c()
  for (j in 1:length(dlist)) {
    
    n_j <- length(dlist[[j]])
    y_j <- dlist[[j]]
    theta_j <- theta_gen(mu_mp1, n_j, y_j, tau2, sig2)
    thetas_mp1 <- c(thetas_mp1, theta_j)
    
  }
  
  sim_vals <- c(thetas_mp1, mu_mp1, sig2, tau2)
  sim_data[i, sim_params] <- sim_vals
  thetas_m <- thetas_mp1
  mu_m <- mu_mp1
  
}

#### results table ---------------------------------------
mu_hats <- colMeans(sim_data)
quants <- sapply(sim_data, function(x) quantile(x, c(0.025, 0.25, 0.75, 0.975)))
results_df <- data.frame('2.5%' = quants[1, ],
                         '25%' = quants[2, ],
                         mean = mu_hats,
                         '75%' = quants[3, ],
                         '97.5%' = quants[4, ],
                         check.names = FALSE) %>% round(digits = 5)
results_df[6:7, ] <- sapply(results_df[6:7, ], function(x) exp(sqrt(x)))
rownames(results_df)[6:7] <- c('sigma', 'tau')

print(results_df)
