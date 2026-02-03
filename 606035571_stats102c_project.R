
rlaplace <- function(n, mu, lambda) {
  out <- numeric(n)
  accept <- 0
  total <- 0
  i <- 1
  
  M <- if (lambda <= 1) {
    pi * (1 + ((1 - sqrt(1 - lambda^2)) / lambda)^2) * lambda / 2 *
      exp(-(1 - sqrt(1 - lambda^2)))
  } else {
    pi * lambda / 2
  }
  
  while (i <= n) {
    x <- rcauchy(1)
    u <- runif(1)
    total <- total + 1
    
    fx <- (lambda / 2) * exp(-lambda * abs(x - mu))
    gx <- 1 / (pi * (1 + x^2))
    
    if (u <= fx / (M * gx)) {
      out[i] <- x
      accept <- accept + 1
      i <- i + 1
    }
  }
  
  attr(out, "accept") <- accept / total
  out
  
}









gibbs_pumps <- function(n, data, burns, params) {
  if (is.null(names(params)) || any(names(params) == "")) stop("params must be a named list/vector")
  need <- c("alpha", "gamma", "delta")
  if (!all(need %in% names(params))) stop("params must contain alpha, gamma, delta")
  
  alpha <- as.numeric(params[["alpha"]])
  gamma <- as.numeric(params[["gamma"]])
  delta <- as.numeric(params[["delta"]])
  
  Y <- data[[1]]
  t <- data[[2]]
  p <- nrow(data)
  
  theta <- rep((alpha + mean(Y)) / (delta + mean(t) + 1e-8), p)
  beta  <- delta
  
  total <- burns + n
  out <- matrix(NA_real_, nrow = n, ncol = p + 1)
  
  for (iter in 1:total) {
    theta <- rgamma(p, shape = alpha + Y, rate = beta + t)
    beta  <- rgamma(1, shape = gamma + p * alpha, rate = delta + sum(theta))
    
    if (iter > burns) {
      out[iter - burns, 1:p] <- theta
      out[iter - burns, p + 1] <- beta
    }
  }
  
  colnames(out) <- c(paste0("theta", 1:p), "beta")
  out
}


