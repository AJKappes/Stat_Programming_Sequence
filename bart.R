# Bayesian Addative Regression Trees #
# Chipman, George, McCulloch (1998) -> CGM98
remove(list = objects())
library(tidyverse)


#### Generate data ----------------------------------
# data generation following Hill (2010)
set.seed(7777)
p <- 0.5
z <- rbinom(200, 1, p)
y <- c()
x <- c()
for (i in 1:length(z)) {
  
  if (z[i] == 1) {
    
    x[i] <- rnorm(1, 40, 10)
    y[i] <- rnorm(1, 90 + exp(0.06*x[i]), 1)
    
  } else {
    
    x[i] <- rnorm(1, 20, 10)
    y[i] <- rnorm(1, 72 + 3*sqrt(x[i]), 1)
    
  }
    
}

set.seed(NULL)
ridx <- which(is.na(y))
data <- data.frame(y = y[-ridx],
                   z = z[-ridx],
                   x = x[-ridx])
remove(x, y, z)

#### Priors and likelihood -------------------------
# back to CGM98

# tree prior
#(alpha, beta) values based on suggestions from CGM98
alpha <- 0.95
beta <- 2
p_split <- function(d, a = alpha, b = beta) {
  
  s <- a*(1 + d)^(-beta)
  # rv <- rbinom(1, 1, s)
  return(s)
  
}

# leaf param prior (mu)
leaf_mu <- function(tn, mu, s2, c) {
  
  mus <- rnorm(tn, mu, sqrt(s2/c))
  return(mus)
  
}

# leaf param hyper
sig2 <- function(v, lambda) {
  
  s2 <- 1/rgamma(v/2, scale = v*lambda/2)
  return(s2)
  
}

# log-likelihood 
loglik <- function(n_vec, y_list, mu_vec, tn, c, v, lambda) {
  
  mu <- mean(mu_vec)
  s_vec <- map2_dbl(y_list, n_vec,
                    function(y, n) (n - 1)*var(y))
  t_vec <- map2_dbl(y_list, n_vec,
                    function(y, n) ((n*c)/(n + c))*(mean(y) - mu)^2)
  
  f <- -sum(n_vec/2)*log(pi) + tn*v/2*log(v*lambda) + tn/2*log(c) +
    -1/2*sum(log(n_vec + c)) + sum(log(gamma((n_vec + v)/2))) - 
    tn*log(gamma(v/2)) - sum((n_vec + v)/2)*sum(log(s_vec + t_vec + v*lambda))
  
  return(f)
  
} 

#### Tree action setup -------------------------------

# initialize z partition nodes and split value
nodes_zl <- list()
nodes_zr <- list()
split_z <- 0.5
# initialize z_left and z_right nodes for > than conditions (go right)
nodes_zlr <- list(NA)
nodes_zrr <- list(NA)
# initialize split value vectors
splits_zl <- c()
splits_zr <- c()

# initiate tree random action
tree_step <- c('grow', 'prune', 'change', 'swap')
draw_action <- function(d) {
  
  # probability of grow
  p <- p_split(d)
  s <- sample(tree_step, 1, prob = c(p, rep((1-p)/3, length(tree_step) - 1)))
  return(s)
  
}

# if action is grow
t_grow <- function() {
  
  # draw initial z node direction, < 0.5 or > 0.5
  s_direct <- sample(c('l', 'r'), 1)
  min_leaf_obs <- 4
  
  if(length(nodes_zl) == 0 | length(nodes_zr) == 0) {
    
    # define intial partition on z 
    nodes_zl[[1]] <<- data %>% filter(z < split_z) %>% .$y
    nodes_zr[[1]] <<- data %>% filter(!y %in% nodes_zl[[1]]) %>% .$y
    # split choosen uniformly, serves as initial splits
    splits_zl[1] <<- sample(data$x, 1)
    splits_zr[1] <<- sample(data$x, 1)
    
  } else {
    
    # ifelse and if_else only handles conditions of length 1
    if (s_direct == 'l') {
      
      
      nodes <- nodes_zl
      nodes_r <- nodes_zlr
      splits <- splits_zl
      
    } else {
      
      nodes <- nodes_zr
      nodes_r <- nodes_zrr
      splits <- splits_zr
      
    }
    
    i <- length(nodes) + 1
      
    # store current split value
    new_split <- data$x %>%
      .[between(., min(.), splits[i - 1])] %>%
      sample(1)
    split <- if_else(length(nodes) <= 1,
                     splits[1],
                     new_split)
    
    # append new split for next iteration
    if (length(nodes) >=1 ) splits[i] <- new_split
    
    # partition data at current split point for data in z_left (< .05)
    split_y <- data %>%
      filter(x < split & y %in% nodes[[1]]) %>%
      .$y
    
    # check if minimum leaf obs is met
    if(length(split_y) > min_leaf_obs) {
      
      # assign data to new node
      nodes[[i]] <- split_y
      # assign data to z_left 'right' node that doesn't satisfy current split rule
      nodes_r[[i]] <- data %>%
        filter(!y %in% unlist(nodes[-1]) &
                 y %in% nodes[[1]]) %>%
        .$y
      
    } else {
      
      splits <- splits[-length(splits)]
      
      cat('Cannot partition.',
          '\nMinimum leaf obs of', min_leaf_obs, 'met.',
          '\n')
      
    }
    
    if (s_direct == 'l') {
      
      nodes_zl <<- nodes
      nodes_zlr <<- nodes_r
      splits_zl <<- splits
      
    } else {
      
      nodes_zr <<- nodes
      nodes_zrr <<- nodes_r
      splits_zr <<- splits
      
    }
    
  }
  
  out <- list(left_nodes = nodes_zl, right_nodes = nodes_zr,
              left_comps = nodes_zlr, right_comps = nodes_zrr,
              left_splits = splits_zl, right_splits = splits_zr,
              direc = s_direct)
  return(out)
  
}

tt <- t_grow()
tt2 <- t_grow()
tt3 <- t_grow()
tt4 <- t_grow()
tt5 <- t_grow()
tt6 <- t_grow()
tt7 <- t_grow()





###########
if (action == 'prune') {
    if(length(nodes_l) == 1 & length(nodes_r) == 1) {
      
      nodes_l <- NULL
      nodes_r <- NULL
      
    }
    
    if (length(nodes_l) > 1 | length(nodes_r) > 1) {
      if (s_direct == 'l') {
        
        nodes_l <- nodes_l[[-c(length(nodes_l))]]
        split_l <- split_l[-c(length(split_l))]
        
      } else {
        
        nodes_r <- nodes_r[[-c(length(nodes_r))]]
        
      }
      
    }
    
  }
  
  if (action == 'change' | action == 'swap') {
    
    cat('Test phase.', action, 'not evaluated yet.')
    
  }





