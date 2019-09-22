library(ggplot2)
library(latex2exp)

#######################################################
#                 problem set 1                       #
#######################################################

# global rv generations
I <- 1000

#######################################################
#                       q4                            #
#######################################################

# simulate 1000 rvs from Bin(3, 0.25) from Ber(0.25)
# rv z = \sum_{i=3}^n x_i for x ~ Ber(0.25)
z_binom <- function(n, p, I) {
  rv_vec <- c()
  
  for (i in 1:I) {
    u <- ifelse(runif(n) > 1-p, 1, 0)
    rv_vec[i] <- sum(u)
  }
  
  return(data.frame('z' = rv_vec))
}

p <- 0.25
n <- 3

ggplot(z_binom(n, p, I), aes(z)) +
  geom_bar() +
  labs(title = TeX('$z \\sim$ Bin(3, 0.25) 1000 RV Generation'),
       x = TeX('$z$'),
       y = TeX('Frequency')) +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'))

#######################################################
#                       q5                            #
#######################################################

# transformation to generate 1000 gamma(4, 1) rvs
# use U(0,1) to simulate x = -ln(u) ~ Exp(1)
# rv y = -beta \sum_{j=1}^aplha ln(u)
y_gam <- function(a, b, I){
  rv_vec = c()
  
  for (i in 1:I){
    x <- -log(runif(a))
    y <- b*sum(x)
    rv_vec[i] <- y
  }
  
  return(data.frame('y' = rv_vec))
}

alpha <- 4
beta <- 1

ggplot(y_gam(alpha, beta, I), aes(y)) +
  geom_histogram(aes(y = ..density..), alpha = 0.25) +
  geom_density(alpha = 0.35, fill = 'black') +
  labs(title = TeX('$y \\sim$ gamma(4, 1) 1000 RV Generation'),
       x = TeX('$y$'),
       y = '') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'))

#######################################################
#                       q8                            #
#######################################################

# simulate rvs N(3, sig^4 = 2)

# box-mueller transform
# simulate 2 U(0,1) for u1 and u2
# x1 = sqrt(-2ln(u1))cos(2piu2), x2 = sqrt(-2ln(u1))sin(2piu2)
# X1, X2 ~ N(0,1) => returns two iid z rvs
# x = mu + sig^(1/4)z

# q4 global variables
mu <- 3
sig4 <- 2

norm_boxmuell <- function(mu, sig, I){
  t <- proc.time()
  z_vec <- c()
  
  for (i in 1:I) {
    u_vec <- runif(2)
    x_vec <- c(sqrt(-2*log(u_vec[1]))*cos(2*pi*u_vec[2]),
               sqrt(-2*log(u_vec[1]))*sin(2*pi*u_vec[2]))
    z_vec[i] <- ifelse(1 %% 2 == 0, x_vec[2], x_vec[1])
  }
  
  bm_x_vec <- mu + sig^(1/4)*z_vec
  boxmuell_elapsed <- proc.time() - t
  out <- list('df' = data.frame('x' = bm_x_vec), boxmuell_elapsed)
  return(out)
}

ggplot(norm_boxmuell(mu, sig4, I)[['df']], aes(x)) +
  geom_histogram(aes(y = ..density..), alpha = 0.25) +
  geom_density(alpha = 0.35, fill = 'black') +
  labs(title = TeX('$x \\sim$ N(3, $\\sqrt{2}$) 1000 RV Generation'),
       subtitle = TeX('Box-Mueller Transform'),
       x = TeX('$x$'),
       y = '') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        plot.subtitle = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'))

# Marsaglia's polar method
norm_marspm <- function(mu, sig, I) {
  t <- proc.time()
  z_vec <- c()
  
  i <- 0
  while (i < I/2) {
    u_vec <- runif(2, -1, 1)
    s <- sum(u_vec^2)
    if (s < 1) {
      x1 <- sqrt(-2*log(s))*(u_vec[1]/sqrt(s))
      x2 <- sqrt(-2*log(s))*(u_vec[2]/sqrt(s))
      z_vec <- rbind(z_vec, x1, x2)
      i <- i + 1
    } else {
      i <- i
    }
  }
  
  x_marspm_vec <- mu + sig^(1/4)*z_vec
  marspm_elapsed <- proc.time() - t
  out <- list('df' = data.frame('x' = x_marspm_vec), marspm_elapsed)
  return(out)
}

ggplot(norm_marspm(mu, sig4, I), aes(x)) +
  geom_histogram(aes(y = ..density..), alpha = 0.25) +
  geom_density(alpha = 0.35, fill = 'black') +
  labs(title = TeX('$x \\sim$ N(3, $\\sqrt{2}$) 1000 RV Generation'),
       subtitle = TeX('Marsaglia\'s Polar Method'),
       x = TeX('$x$'),
       y = '') +
  theme(plot.title = element_text(family = 'serif', hjust = 0.5),
        plot.subtitle = element_text(family = 'serif', hjust = 0.5),
        axis.title = element_text(family = 'serif'))

# algorithm cost comparison
marspm_t <- round(norm_marspm(mu, sig4, I)[2][[1]][['elapsed']], 5)
bm_t <- round(norm_boxmuell(mu, sig4, I)[2][[1]][['elapsed']], 5)

if (marspm_t < bm_t) {
  cat('Masaglia\'s polar method more is effeicient.\nProcessing time:',
      marspm_t, '<', bm_t)
} else {
  cat('Box-Mueller transform method is more efficient.\nProcessing time:',
      bm_t, '<', marspm_t)
}

