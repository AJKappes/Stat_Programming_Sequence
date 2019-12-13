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

#### Tree action setup -----------------------------------

# initialize z partition nodes and split value
nodes_zl <- list()
nodes_zr <- list()
split_z <- 0.5
# initialize z_left and z_right nodes for > than conditions (go right)
nodes_zlr <- list(NA)
nodes_zrr <- list(NA)
# initialize split value vectors
splits_zl <- c(NA)
splits_zr <- c(NA)
# initialize leaf size constraint
min_leaf_obs <- 4
# initialize leaf parameter vectors
leaf_params <- setNames(data.frame(matrix(0, nrow = 1, ncol = 4)),
                        c('l', 'lcomp', 'r', 'rcomp'))
# initiate tree random action
tree_step <- c('grow', 'prune', 'change', 'swap')
draw_action <- function(d) {
  
  # probability of grow
  p <- p_split(d)
  s <- sample(tree_step, 1, prob = c(p, rep((1-p)/3, length(tree_step) - 1)))
  return(s)
  
}

#### Tree action functions --------------------------------

# if action is grow
t_grow <- function() {
  
  # draw initial z node direction, < 0.5 or > 0.5
  s_direct <- sample(c('l', 'r'), 1)
  
  if(length(nodes_zl) == 0 | length(nodes_zr) == 0) {
    
    # define intial partition on z 
    nodes_zl[[1]] <<- data %>%
      filter(z < split_z) %>%
      .$y
    nodes_zr[[1]] <<- data %>%
      filter(!y %in% nodes_zl[[1]]) %>%
      .$y
    
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
    
    # define uniform splitting rultes
    if (length(nodes) == 1) {
      
      splits[i] <- sample(data$x, 1)
      
    } else {
      
      splits[i] <- data$x %>%
        .[between(., min(.), splits[i - 1])] %>%
        sample(1)
      
    }
    
    split <- splits[i]
    
    # partition data at current split point for data in z <> .05
    split_y <- data %>%
      filter(x < split & y %in% nodes[[1]]) %>%
      .$y
    
    # check if minimum leaf obs is met
    if(length(split_y) >= min_leaf_obs) {
      
      # assign data to new node
      nodes[[i]] <- split_y
      
      # parition compliment data depending on node length
      if (length(nodes) <= 2) {
        
        split_comp_y <- data %>%
          filter(!y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      } else {
        
        split_comp_y <- data %>%
          filter(between(x, splits[i], splits[i - 1]) &
                   !y %in% unlist(nodes_r) &
                   y %in% nodes[[1]]) %>%
          .$y
        
      }
      
      # assign data to compliment node
      nodes_r[[i]] <- split_comp_y
      
    } else {
      
      splits <- splits[-length(splits)]
      cat('Cannot partition.',
          '\nMinimum leaf obs of', min_leaf_obs, 'not met.',
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

# if action is prune
t_prune <- function() {
  
  s_direct <- sample(c('l', 'r'), 1)
  
  if(length(nodes_zl) == 1  & s_direct == 'l') {
    
    cat('Cannot prune', s_direct, 'any further.')
    
  } else if (length(nodes_zr) == 1  & s_direct == 'r') {
    
    cat('Cannot prune', s_direct, 'any further.')
    
  } else {
    
    if (s_direct == 'l') {
      
      i <- length(nodes_zl)
      nodes_zl <<- nodes_zl[-i]
      nodes_zlr <<- nodes_zlr[-i]
      splits_zl <<- splits_zl[-i]
      
    } else {
      
      i <- length(nodes_zr)
      nodes_zr <<- nodes_zr[-i]
      nodes_zrr <<- nodes_zrr[-i]
      splits_zr <<- splits_zr[-i]
      
    }
    
  }
  
  out <- list(left_nodes = nodes_zl, right_nodes = nodes_zr,
              left_comps = nodes_zlr, right_comps = nodes_zrr,
              left_splits = splits_zl, right_splits = splits_zr,
              direc = s_direct)
  return(out)
  
}

# if action is change
t_change <- function() {
  
  s_direct <- sample(c('l', 'r'), 1)
  
  if(length(nodes_zl) == 1  & s_direct == 'l') {
    
    cat('Cannot change splitting rule for', s_direct,
        '\nCurrent node is the initial node.')
    
  } else if (length(nodes_zr) == 1  & s_direct == 'r') {
    
    cat('Cannot change splitting rule for', s_direct,
        '\nCurrent node is the initial node.')
    
  } else {
    
    if (s_direct == 'l') {
      
      nodes <- nodes_zl
      nodes_r <- nodes_zlr
      splits <- splits_zl
      
    } else {
      
      nodes <- nodes_zr
      nodes_r <- nodes_zrr
      splits <- splits_zr
      
    }
    
    # random node draw conditioned on random side (z <> 0.5)
    j <- sample(2:length(nodes), 1)
    # subset nodes to operate on random node draw
    nodes <- nodes[-c(j + 1:length(nodes))]
    nodes_r <- nodes_r[-c(j + 1:length(nodes_r))]
    splits <- splits[-c(j + 1:length(splits))]
    
    # define operating index
    i <- length(nodes)
    
    # assign new splitting rule at the random draw
    splits[i] <- data$x %>%
      .[between(., min(.), splits[i - 1])] %>%
      sample(1)
    split <- splits[i]
    
    # partition data at current split point for data in z <> .05
    split_y <- data %>%
      filter(x < split & y %in% nodes[[1]]) %>%
      .$y
    
    # check minimum leaf obs condition
    if(length(split_y) >= min_leaf_obs) {
      
      # assign data to new random node
      nodes[[i]] <- split_y
      
      if (length(nodes) <= 2) {
        
        split_comp_y <- data %>%
          filter(!y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      } else {
        
        split_comp_y <- data %>%
          filter(between(x, splits[i], splits[i - 1]) &
                   !y %in% unlist(nodes_r) &
                   y %in% nodes[[1]]) %>%
          .$y
        
      }
      
      # assign data to new random compliment node
      nodes_r[[i]] <- split_comp_y
      
    } else {
      
      # need to figure out what happens to node when cond not met
      
      cat('Cannot partition at random node with new splitting rule.',
          '\nMinimum leaf obs of', min_leaf_obs, 'not met.',
          '\n')
      
    }
      
  }
  
}


# check processes
tt <- list('z node', c(1:4), c(5:9), c(10:14), c(15:19))

