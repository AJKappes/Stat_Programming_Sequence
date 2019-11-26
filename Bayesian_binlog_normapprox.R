# Bayesian binomial-logistic model for multivariate X
# Binomal GLM with logistic link
# Taylor series expansion about center XBeta for realizations z and variance

# Model is used to test the hypothesis:
## Rural Kenya livestock holders substitute ex-post livestock treatment
## with ex-ante vaccination 

remove(list = objects())
library(ggplot2)
library(latex2exp)

#### observed data ###

df_bovine <- read.csv('~/research/africa/Nutrient_Demand/nd_R/df_bovine.csv')
y <- matrix(df_bovine$AnimalVaccinated)
y[which(y == 77)] <- 0
X <- matrix(cbind(1,
                  df_bovine$AnimalTreated,
                  df_bovine$livincome,
                  df_bovine$OffFarmNetIncome),
            nrow = nrow(df_bovine))
# change 'unkown' value to 0
X[which(X[, 2] == 99), 2] <- 0


### define functions and conduct optimization ###

# pseudodata functions resulting from expansion
z_gen <- function(x, b, y, n = length(y)) {
  
  eta <- x%*%b
  z <- eta + (1 + exp(eta))^2/exp(eta)*(y/n - exp(eta)/(1 + exp(eta)))
  return(z)
  
}

sig2_gen <- function(x, b, n = length(y)) {
  
  eta <- x%*%b
  s <- 1/n*(1 + exp(eta))^2/exp(eta)
  return(s)
  
}

beta_approx <- function(x, z, w) {
  
  params <- solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%z
  # params <- coef(lm(z ~ x - 1, weights = w))
  return(params)
  
}

pred_prob <- function(x, b){
  
  eta <- x%*%b
  p <- exp(eta)/(1 + exp(eta))
  return(p)
  
}

# iterative normal center approximation
err <- 1
tol <- 1E-11
itr <- 5
max_itr <- 50
beta_m <- matrix(rep(0, dim(X)[2]), nrow = dim(X)[2])
beta_vals <- data.frame()
while (err > tol & itr <= max_itr) {
  
  z <- z_gen(X, beta_m, y)
  sig2 <- sig2_gen(X, beta_m)
  W <- diag(1/as.vector(sig2))
  # W <- as.vector(sig2)
  beta_mp1 <- beta_approx(X, z, W)
  err <- max(abs(beta_mp1 - beta_m))
  
  beta_vals[itr, c('b0', 'b1', 'b2', 'b3')] <- c(beta_mp1)
  beta_m <- beta_mp1
  
  cat('Iteration', itr, 'max optimization error:', err, '\n')
  
  itr <- itr + 1
  
}

# estimator convergence plot
est_conv_plot <- ggplot(beta_vals, aes(1:nrow(beta_vals))) +
  geom_line(aes(y = b0, colour = 'b_0')) +
  geom_line(aes(y = b1, colour = 'b_1')) +
  geom_line(aes(y = b2, colour = 'b_2')) +
  geom_line(aes(y = b3, colour = 'b_3')) +
  scale_color_manual(name = '',
                     values = c('b_0' = 'black', 'b_1' = 'orange',
                                'b_2' = 'blue', 'b_3' = 'red')) +
  labs(title = 'Estimator Convergence',
       x = 'Iteration',
       y = TeX('$\\beta$')) +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'),
        legend.text = element_text(family = 'serif'))

# predict conditional probabilities
vacc <- c(1, 0)
cond <- c('treatment', 'no treatment')
phat <- c()
for (i in 1:length(vacc)) {
  
  pred_X <- matrix(c(1, vacc[i], mean(X[, 3]), mean(X[, 4])), nrow = 1)
  pred_b <- beta_mp1
  phat[i] <- pred_prob(pred_X, pred_b)
  
  cat('Predicted Pr( vaccinate |', cond[i], ') = ', phat[i], '\n')
  
}

# approximate normal beta posterior distribution
post_var <- function(x, b, n = length(y)) {
  
  eta <- x%*%b
  f <- c(n*(exp(eta)/(1 + exp(eta))^2))
  F <- diag(f)
  V <- solve(t(x)%*%F%*%x)
  return(V)
  
}

# create beta distributions and CIs
beta_var <- diag(post_var(X, beta_mp1))
nsim <- 1E+3
beta_dist <- setNames(data.frame(matrix(nrow = nsim, ncol = dim(beta_vals)[2])),
                      names(beta_vals))
beta_ci <- data.frame(lb = numeric(),
                      ub = numeric())
for (i in 1:length(beta_mp1)) {
  
  beta_dist[, names(beta_dist)[i]] <- rnorm(nsim, beta_mp1[i], sqrt(beta_var[i]))
  beta_ci[i, 'lb'] <- beta_mp1[i] - qt(1 - .05/2, dim(X)[1] - dim(X)[2])*sqrt(beta_var[i])
  beta_ci[i, 'ub'] <- beta_mp1[i] + qt(1 - .05/2, dim(X)[1] - dim(X)[2])*sqrt(beta_var[i])
  
}


