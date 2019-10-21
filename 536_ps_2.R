library(data.table)
library(ggplot2)
library(latex2exp)

# ps 2

# global critical t-value
t_crit <- qt(1-.05/2, 1E+4)

#### q1 ####


# 1d
# target f(x) = 1-|x|
# instrumental g(x) = 1/2 for U(-1,1)
# M = 2
n <- 100

condition_fun <- function(x){
  
  fx <- 1 - abs(x)
  gx <- 1/2
  M <- 2
  
  return(fx/(M*gx))
  
}

y_vec <- c()
i <- 0
while (i <= n) {
  
  u <- runif(1)
  x <- runif(1, -1, 1)
  
  if (u <= condition_fun(x)) {
    
    y_vec[i] <- x
    i <- i + 1
    
  } 
  
  else {
    
    i <- i
    
  }
  
}

y_dt <- data.table(rv_y = y_vec)

ggplot(y_dt, aes(rv_y)) +
  geom_histogram(aes(y = ..density..), alpha = 0.25, bins = 8) +
  geom_density(alpha = 0.35, fill = 'black', colour = 'grey37') +
  xlim(-1, 1) +
  labs(title = 'Accept-Reject Method',
       x = TeX('$y$'),
       y = '') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'))

# 1c
condition_fun <- function(m, u){
  
  if (m < 0) {
    
    x = -1 + sqrt(2*u)
    
  } else {
    
    x = 1 - sqrt(2*u)
    
  }
  
  return(x)
  
}

x_vec <- c()
i <- 0
while (i < n) {
  
  m <- runif(1, -1, 1)
  u <- runif(1)
  x_vec[i] <- condition_fun(m, u)
  i <- i + 1
  
}

x_dt <- data.table(rv_x = x_vec)

ggplot(x_dt, aes(rv_x)) +
  geom_histogram(aes(y = ..density..), alpha = 0.25, bins = 8) +
  geom_density(alpha = 0.35, fill = 'black', colour = 'grey37') +
  xlim(-1, 1) +
  labs(title = 'Inverse Method',
       x = TeX('$x$'),
       y = '') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'))

#### q2 ####
# z ~ N(0,1), p=P(Z>3)

# q2a
# monte carlo integration
n <- 1E+4
x_o <- rnorm(n)
p_hat <- length(x_o[x_o > 3]) / n # what we should converge to => 0.0011

p_ind <- ifelse(x_o > 3, 1, 0) # indicator

p_var <- function(n) {
  
  return(1/(4*n))
  
}

p_conv <- c()
p_v <- c()
p_upper <- c()
p_lower <- c()
for (i in 1:n) {
  
  p_conv[i] <- mean(p_ind[1:i])
  p_v[i] <- p_var(i)
  p_upper[i] <- p_conv[i] + t_crit*sqrt(p_v[i])
  p_lower[i] <- p_conv[i] - t_crit*sqrt(p_v[i])
  
}

q2a_dt <- data.table(p = p_conv, u = p_upper, l = p_lower, itr = seq(1, n, 1))

ggplot(q2a_dt, aes(itr)) +
  geom_line(aes(y = p, colour = 'mean'), alpha = 0.75) +
  geom_line(aes(y = u, colour = 'upper')) +
  geom_line(aes(y = l, colour = 'lower')) +
  scale_color_manual(name = '',
                     values = c('mean' = 'black', 'upper' = 'gray', 'lower' = 'gray')) +
  geom_ribbon(aes(ymin = l, ymax = u), fill = 'gray', alpha = 0.4) +
  ylab(TeX('$\\hat{p}$ convergence')) +
  xlab('Iteration') +
  ggtitle(TeX('$P(Z>3)$, $Z\\sim N(0, 1)$ Estimator Convergence'),
          subtitle = 'Monte Carlo Approximation') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        plot.subtitle = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'),
        legend.position = 'none')

# q2b
# importance sampling
# g(x) shifted exponential e^-(x - s), s shift
# f(x) target standard normal 1/(2pi)*e^(-x^2/2)

IS_fun <- function(x, s) {
  
  gx <- exp(-(x-s)) # instrumental shifted exponential
  fx <- 1/(2*pi)*exp(-x^2/2) # target standard normal
  hx <- ifelse(x > s, 1, 0) # indicator 
  
  return(fx/gx*hx) # IS estimator
  
}

IS_var <- function(n, mu, x, s) {
  
  v <- 1/n*sum((IS_fun(x, s) - mu)^2) # sample variance
  return(v)
  
}

x <- rexp(n) + 3 # gx shifted exponential, shift = 3
IS_vec <- IS_fun(x, 3)
is_conv <- c()
is_upper <- c()
is_lower <- c()
is_v <- c()
for (i in 1:n) {
  
  is_conv[i] <- mean(IS_vec[1:i])
  is_v[i] <- IS_var(i, is_conv[i], x[1:i], 3)
  is_upper[i] <- is_conv[i] + t_crit*sqrt(is_v[i]/i)
  is_lower[i] <- is_conv[i] - t_crit*sqrt(is_v[i]/i)
  
}

q2b_dt <- data.table(p = is_conv, u = is_upper, l = is_lower, itr = seq(1, n, 1))

ggplot(q2b_dt, aes(itr)) +
  geom_line(aes(y = p, colour = 'mean'), alpha = 0.75) +
  geom_line(aes(y = u, colour = 'upper')) +
  geom_line(aes(y = l, colour = 'lower')) +
  scale_color_manual(name = '',
                     values = c('mean' = 'black', 'upper' = 'gray', 'lower' = 'gray')) +
  geom_ribbon(aes(ymin = l, ymax = u), fill = 'gray', alpha = 0.4) +
  ylab(TeX('$\\hat{p}$ convergence')) +
  xlab('Iteration') +
  ggtitle(TeX('$P(Z>3)$, $Z\\sim N(0, 1)$ Estimator Convergence'),
          subtitle = 'Importance Sampling') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        plot.subtitle = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'),
        legend.position = 'none')


#### q3 ####
# x ~ cauchy(0, 1), p=P(0<Z<2)

# q3a
# monte carlo approximation
n <- 1E+3

# cauchy sim with BM method
cauchy_sim <- function(n){
  
  # box-mueller n(0,1) sim
  x1_vec <- c()
  x2_vec <- c()
  
  for (i in 1:n) {
    
    u_vec <- runif(2)
    x1_vec[i] <- sqrt(-2*log(u_vec[1]))*cos(2*pi*u_vec[2])
    x2_vec[i] <- sqrt(-2*log(u_vec[1]))*sin(2*pi*u_vec[2])
    
  }
  cauchy_vec <- x1_vec / x2_vec
  return(cauchy_vec)
}

cauchy_rvs <- cauchy_sim(n)
cauchy_ind <- ifelse(cauchy_rvs > 0 & cauchy_rvs < 2, 1, 0)
p_hat <- length(cauchy_rvs[cauchy_rvs > 0 & cauchy_rvs < 2]) / n # to check convergence p_hat = 0.369

cauchy_var <- function(j, p) {
  
  return(p*(1-p)/j)
  
}

# mc approximation
c_conv <- c()
c_upper <- c()
c_lower <- c()
c_var <- c()
for (i in 1:n) {
  
  c_conv[i] <- mean(cauchy_ind[1:i])
  c_var[i] <- cauchy_var(i, c_conv[i])
  c_upper[i] <- c_conv[i] + t_crit*sqrt(c_var[i])
  c_lower[i] <- c_conv[i] - t_crit*sqrt(c_var[i])
  
}

q3a_dt <- data.table(p = c_conv, u = c_upper, l = c_lower, itr = seq(1, n, 1))

ggplot(q3a_dt, aes(itr)) +
  geom_line(aes(y = p, colour = 'mean'), alpha = 0.75) +
  geom_line(aes(y = u, colour = 'upper')) +
  geom_line(aes(y = l, colour = 'lower')) +
  scale_color_manual(name = '',
                     values = c('mean' = 'black', 'upper' = 'gray', 'lower' = 'gray')) +
  geom_ribbon(aes(ymin = l, ymax = u), fill = 'gray', alpha = 0.4) +
  ylab(TeX('$\\hat{p}$ convergence')) +
  xlab('Iteration') +
  ggtitle(TeX('$P(0<Z<2)$, $X\\sim C(0, 1)$ Estimator Convergence'),
          subtitle = 'Monte Carlo Approximation') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        plot.subtitle = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'),
        legend.position = 'none')

# q3b
# instrumental U(0,2) pdf 
# target C(0,1) pdf 1/(pi*(1+x^2))

target_dens <- function(x) {
  
  return(1/(pi*(1 + x^2)))
  
}

IS_fun <- function(x) {
  
  gx <- 1/2
  fx <- target_dens(x) # target cauchy(0, 1)
  hx <- ifelse(x > 0 & x < 2, 1, 0) # indicator
  
  return(fx/gx*hx) # IS estimator
  
}

IS_var <- function(n, mu, x) {
  
  v <- 1/n*sum((IS_fun(x) - mu)^2) # sample variance
  return(v)
  
} 

x <- runif(n, 0, 2)
IS_vec <- IS_fun(x)
is_conv <- c()
is_v <- c()
is_upper <- c()
is_lower <- c()
for (i in 1:n) {
  
  is_conv[i] <- mean(IS_vec[1:i])
  is_v[i] <- IS_var(i, is_conv[i], x[1:i])
  is_upper[i] <- is_conv[i] + t_crit*sqrt(is_v[i]/i)
  is_lower[i] <- is_conv[i] - t_crit*sqrt(is_v[i]/i)
  
}

q3b_dt <- data.table(p = is_conv, u = is_upper, l = is_lower, itr = seq(1, n, 1))

ggplot(q3b_dt, aes(itr)) +
  geom_line(aes(y = p, colour = 'mean'), alpha = 0.75) +
  geom_line(aes(y = u, colour = 'upper')) +
  geom_line(aes(y = l, colour = 'lower')) +
  scale_color_manual(name = '',
                     values = c('mean' = 'black', 'upper' = 'gray', 'lower' = 'gray')) +
  geom_ribbon(aes(ymin = l, ymax = u), fill = 'gray', alpha = 0.4) +
  ylab(TeX('$\\hat{p}$ convergence')) +
  xlab('Iteration') +
  ggtitle(TeX('$P(0<Z<2)$, $X\\sim C(0, 1)$ Estimator Convergence'),
          subtitle = 'Importance Sampling') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        plot.subtitle = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'),
        legend.position = 'none')

#### q4 ####
# x ~ N(2, 9)
# h(x) = 3x-1
# mc approximate E[3x-1]

# use BM to generate z ~ N(0,1)
# x = mu + sig*z
n <- 1E+3

gen_normal <- function(n, mu, sig2) {
  
  x1_vec <- c()
  x2_vec <- c()
  
  for (i in 1:n) {
    
    u_vec <- runif(2)
    x1_vec[i] <- sqrt(-2*log(u_vec[1]))*cos(2*pi*u_vec[2])
    x2_vec[i] <- sqrt(-2*log(u_vec[1]))*sin(2*pi*u_vec[2])
    
  }
  
  z <- append(x1_vec, x2_vec)
  x <- mu + sqrt(sig2)*z
  
  return(x)
  
}

hx <- function(x) {
  
  return(3*x-1)
  
}

x <- gen_normal(n, mu = 2, sig2 = 9)
Ehx_conv <- c()
hx_v <- c()
Ehx_upper <- c()
Ehx_lower <- c()
for (i in 1:n) {
  
  Ehx_conv[i] <- mean(hx(x[1:i]))
  hx_v[i] <- 3^2*9/i
  Ehx_upper[i] <- Ehx_conv[i] + t_crit*sqrt(hx_v[i])
  Ehx_lower[i] <- Ehx_conv[i] - t_crit*sqrt(hx_v[i])
  
}

q4_dt <- data.table(Ehx = Ehx_conv, u = Ehx_upper, l = Ehx_lower, itr = seq(1, n, 1))

ggplot(q4_dt, aes(itr)) +
  geom_line(aes(y = Ehx, colour = 'mean'), alpha = 0.75) +
  geom_line(aes(y = u, colour = 'upper')) +
  geom_line(aes(y = l, colour = 'lower')) +
  scale_color_manual(name = '',
                     values = c('mean' = 'black', 'upper' = 'gray', 'lower' = 'gray')) +
  geom_ribbon(aes(ymin = l, ymax = u), fill = 'gray', alpha = 0.4) +
  ylab(TeX('$E\\[h(x)\\]$ convergence')) +
  xlab('Iteration') +
  ggtitle(TeX('Expectation Convergence for $h(x)=3x-1$, $X \\sim N(2, 9)$'),
          subtitle = 'Monte Carlo Approximation') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        plot.subtitle = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'),
        legend.position = 'none')

#### q5 ####

# to replicate book results
book_dt <- data.table(y = c(21, 15), n = c(2, 3))
n <- sum(book_dt)
# q5_dt <- data.table(y = c(34, 6), n = c(19, 16))
# n <- sum(q5_dt)

# generate p1, p2, and associated y
p <- runif(2)
p_vec <- c(p[1]*p[2], p[1]*(1-p[2]), p[2]*(1-p[1]), (1-p[1])*(1-p[2]))
y_rv <- as.vector(rmultinom(1, size = n, prob = p_vec))

p_vec_null <- function(y) {
  
  p1 <- (y[1]+y[2])/n
  p2 <- (y[1]+y[3])/n
  p11 <- p1*p2 # restriction
  p12 <- y[2]/n
  p21 <- y[3]/n
  p22 <- y[4]/n
  
  return(c(p11, p12, p21, p22))
}

p_vec_alt <- function(y) {
  
  p11 <- y[1]/n
  p12 <- y[2]/n
  p21 <- y[3]/n
  p22 <- y[4]/n
  
  return(c(p11, p12, p21, p22))
}

lam <- function(y) {
  
  llh_n <- prod(p_vec_null(y)^y)
  llh_a <- prod(p_vec_alt(y)^y)
  return(llh_n/llh_a)
  
}

nsim <- 1E+4
lam_y <- c()
for (i in 1:nsim) {
  
  p <- runif(2)
  p_vec <- c(p[1]*p[2], p[1]*(1-p[2]), p[2]*(1-p[1]), (1-p[1])*(1-p[2]))
  y_rv <- as.vector(rmultinom(1, size = n, prob = p_vec))
  lam_y[i] <- lam(y_rv)
  
}

sim_lam_dt <- data.table(lam = lam_y)

ggplot(sim_lam_dt, aes(lam)) +
  geom_histogram(alpha = 0.5, bins = 100) +
  xlim(0, 50) +
  xlab(TeX('$\\lambda_i(y_i)$')) +
  ylab('Frequency') +
  ggtitle(TeX('Monte Carlo $H_o$ Distribution')) +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'))

alphas <- c(.1, .05, .01)
crit_chi <- c()
crit_monte <- c()
for (i in 1:length(alphas)) {
  
  crit_chi[i] <- round(qchisq(1-alphas[i], 1), 5)
  crit_monte[i] <- round(quantile(-2*log(lam_y), 1-alphas[i]), 4)
  
}

cutoff_dt <- data.table(monte = crit_monte, chi = crit_chi)

#### q6 ####
n <- 1E+4

IS_fun <- function(x) {
  
  # target truncated normal
  fx <- exp((-x^2/2)) * ifelse(x > 1, 1, 0) 
  
  # instrumental shifted exponential
  gx <- exp(-(x-1)) * ifelse(x > 1, 1, 0)
  
  # h function
  hx <- ifelse(x > 1.5, 1, 0)
  
  return(fx/gx*hx)
  
}

IS_var <- function(n, x, mu) {
  
  v <- 1/n*sum((x - mu)^2)
  return(v)
  
}

x <- rexp(n) + 1 # gx shifted exponential, shift = 1
IS_vec <- IS_fun(x)
is_conv <- c()
is_v <- c()
is_upper <- c()
is_lower <- c()
for (i in 1:n) {
  
  is_conv[i] <- mean(IS_vec[1:i])
  is_v[i] <- IS_var(i, is_conv[i], x[1:i])
  is_upper[i] <- is_conv[i] + t_crit*sqrt(is_v[i]/i)
  is_lower[i] <- is_conv[i] - t_crit*sqrt(is_v[i]/i)
  
}

q6_dt <- data.table(p = is_conv, u = is_upper, l = is_lower, itr = seq(1, n, 1))


ggplot(q6_dt, aes(itr)) +
  geom_line(aes(y = p, colour = 'mean'), alpha = 0.75) +
  geom_line(aes(y = u, colour = 'upper')) +
  geom_line(aes(y = l, colour = 'lower')) +
  scale_color_manual(name = '',
                     values = c('mean' = 'black', 'upper' = 'gray', 'lower' = 'gray')) +
  geom_ribbon(aes(ymin = l, ymax = u), fill = 'gray', alpha = 0.4) +
  ylab(TeX('$\\hat{p}$ convergence')) +
  xlab('Iteration') +
  ggtitle(TeX('Expectation Convergence for $P(x > 1.5)$'),
          subtitle = 'Monte Carlo Approximation') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        plot.subtitle = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'),
        legend.position = 'none')

#### q7 ####
n <- 1E+3
pdim <- 3
Sigma <- matrix(c(2, -1, 0.5, -1, 4, 1, 0.5, 1, 5), nrow = 3, ncol = 3)

# check positive def
eigen(Sigma)$values > 0

# given positive def matrix A
# choleski factorization A = LL'
# where L is lower triangle with diag(l1,...lp) > 0
Gam <- t(chol(Sigma))
mu <- c(-2, 4, 3)

# generation approach
# generate p-dimensial z ~ N(0, Ip)
# make transformation x_p = mu_p + L * z_p

# generate z using BM method
sim_x_p <- function(p) {
  
  z1_vec <- c()
  z2_vec <- c()
  
  # generate z
  for (i in 1:p) {
    
    u_vec <- runif(2)
    z1_vec[i] <- sqrt(-2*log(u_vec[1]))*cos(2*pi*u_vec[2])
    z2_vec[i] <- sqrt(-2*log(u_vec[1]))*sin(2*pi*u_vec[2])
    
  }
  
  # extract p-dimensional z
  z_samp <- sample(c(z1_vec, z2_vec), p)
  z <- matrix(z_samp, nrow = p, ncol = 1)
  
  # compute x transformation
  x_pdim <- mu + Gam %*% z
  # return (3x1) dimensional x = (x11, x12, x13)
  return(t(x_pdim))
  
}

# simulate required multivariate sample
X <- c()
for (i in 1:n) {
  
  X <- rbind(X, sim_x_p(pdim))
  
}

mu_vec <- c()
var_vec <- c()
for (i in 1:dim(X)[2]) {
  
  mu_vec[i] <- mean(X[, i])
  var_vec[i] <- var(X[, i])
  
}

# mc approximation of cov(X)
cov_list <- list()
for (i in 1:n) {
  
  cov_list[i] <- list(cov(matrix(X[1:i, ], nrow = i, ncol = pdim)))
  
}

# approximate variances
x_var_dt <- data.table(v_x1 = rep(0, n-1),
                       v_x2 = rep(0, n-1),
                       v_x3 = rep(0, n-1))

for (i in 2:n) {
  
  for (j in 1:pdim) {
    
    x_var_dt[i-1, j] <- cov_list[[i]][j, j]
    
  }
 
}

# approximate covariances
c_x1x2 <- c()
c_x1x3 <- c()
c_x2x3 <- c()
for (i in 2:n) {
  
  c_x1x2 <- c(c_x1x2, cov_list[[i]][1, 2])
  c_x1x3 <- c(c_x1x3, cov_list[[i]][1, 3])
  c_x2x3 <- c(c_x2x3, cov_list[[i]][2, 3])
  
}

# construct cov matrix
cov_vals <- c(mean(c_x1x2), mean(c_x1x3), mean(c_x2x3))
mc_cov <- diag(colMeans(x_var_dt))
mc_cov[lower.tri(mc_cov)] <- cov_vals
mc_cov[upper.tri(mc_cov)] <- mc_cov[lower.tri(mc_cov)]
