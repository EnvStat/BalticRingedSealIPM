

logit <- function(x) { return(log(x/(1-x))) } # logit transform
logistic <- function(x) { return(1 / (1 + exp(-x))) } # inverse logit transform

# compute asymptotic population growth rate by finding the root of the Euler-Lotka equation
lambda <- function(phi, b) {
  root <- optimize(function(lambda, params) {
    return( (lambda^5 - phi[6]*lambda^4 - b/2*prod(phi[1:5]))^2 )
    }, c(1, 2), tol=1e-6)$minimum
  
  return(root)
}

# compute stable demographic structure at carrying capacity
# stable structure for exponentially growing population can be computed by setting phi => phi/lambda
stable.structure.eq <- function(phi) {
  a <- 6
  
  # female structure
  l.f <- exp(cumsum(log(c(1, phi[1:(a-1)]))))
  l.f[a] <- l.f[a] / (1-phi[a])
  
  # male structure
  l.m <- exp(cumsum(log(c(1, phi[(a+1):(2*a-1)]))))
  l.m[a] <- l.m[a] / (1-phi[2*a])
  
  l <- c(l.f, l.m) # combined demographic structure
  
  return(l / sum(l)) # return normalized demographic structure
}

