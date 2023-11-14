library(mvtnorm)
library(invgamma)
library(LaplacesDemon)

#====================================================

# beta

beta.sample = function(X, Sigma.B, Sigma.I, Sigma.S, Binv, Sinv, Iinv, sigma.sq, Y.list, z, lambda, mu.beta, T){
  
  L = nrow(Sigma.B) 
  n = nrow(Sigma.S)
  
  lambda = as.vector(lambda)
  
  # Binv = solve(Sigma.B)
  # Sinv = solve(Sigma.S)
  # Iinv = solve(Sigma.I)
  
  Sigma.B.star = Sigma.S %x% Sigma.I %x% solve( Binv + crossprod( t(X)/ sqrt(sigma.sq) ) )
  
  term2 = matrix(0, ncol = 1, nrow = nrow(Sigma.B.star))
  
  for(i in 1:T){
    Y.cur <- as.vector(  t(   Y.list[[i]]   )   )
    term2 = term2 + (1/sigma.sq[i]) * ( ( Sinv %x% Iinv ) %*% (Y.cur - abs(z[i]) * rep(1,n) %x% lambda )   ) %x% X[,i]
  }
  
  mu.beta.star = Sigma.B.star %*% ( term2 + (    (Sinv%*%rep(1,n)) %x% (Iinv %*% mu.beta) %x% (Binv %*% rep(1,L))   ) )
  
  out = mvnfast::rmvn(1, mu.beta.star, Sigma.B.star, isChol = TRUE)
  
  return(out)
}

#====================================================

# mu_beta

mu_beta.sample = function(Sigma.S, Sigma.B, Sigma.I, Binv, Sinv, Iinv, P, z, B.beta){
  
  L = nrow(Sigma.B)
  n = nrow(Sigma.S)
  
  # Binv = solve(Sigma.B)
  # Sinv = solve(Sigma.S)
  # Iinv = solve(Sigma.I)
  
  sig.s.1 = Sinv %*% rep(1,n)
  sig.b.1 = Binv %*% rep(1,L)
  
  Sigma.mu.B = solve( as.vector(t(rep(1,n)) %*% sig.s.1) * as.vector( t(rep(1,L)) %*% sig.b.1 ) * Iinv + 1e-4 * diag(P) )
  mu.B.bar = Sigma.mu.B %*% Iinv %*% t(B.beta) %*% ( sig.s.1 %x% sig.b.1 )
  
  out = matrix(mvnfast::rmvn(1, mu.B.bar, Sigma.mu.B, isChol = TRUE), ncol = 1)
  
  return(out)
}
#mu_beta.sample(Sigma.S, Sigma.B, Sigma.I, P, z, sigma2, B, mu.beta)

#======================================================

# lambda

lambda.sample = function(Sigma.S, z, sigma.sq, Sigma.I, Sinv, Iinv, T, Y.list, mu.star.list){
  
  n = nrow(Sigma.S)
  
  # Sinv = solve(Sigma.S)
  # Iinv = solve(Sigma.I)
  
  sig.s.1 = Sinv %*% rep(1,n)
  
  Sigma.lambda.star = solve( as.numeric((t(rep(1,n)) %*% sig.s.1)) * sum(z^2/sigma.sq) * Iinv + 1e-2 * diag(P) )
  
  term = matrix(0, P, n)
  
  for(i in 1 : T){
    term = term + (abs(z[i])/sigma.sq[i]) * t(  Y.list[[i]] - mu.star.list[[i]]  ) 
  }
  
  mu.lambda.star = Sigma.lambda.star %*% ( Iinv %*% term %*% sig.s.1 )
  
  out = mvnfast::rmvn(1, mu.lambda.star, Sigma.lambda.star, isChol = TRUE)
  return(out)
}
# lambda.sample(Sigma.S, z, sigma, Sigma.I, T, Y.list, mu.star.list)

#=========================================================

# |zt|

zt.sample = function(sigma.sq.t, Sigma.S, Sinv, Iinv, lambda, Sigma.I, Yt, mu.t){
  
  n = nrow(Sigma.S)
  
  lambda = as.vector(lambda)
  
  # Sinv = solve(Sigma.S)
  # Iinv = solve(Sigma.I)
  
  sig.s.1 = Sinv %*% rep(1,n)
  sig.s.1.1 = t(rep(1,n))%*% sig.s.1
  
  sig.lambda = Iinv%*%lambda
  sig.lambda.lamda = t(lambda)%*%sig.lambda
  
  val <- as.numeric(   1/(1+ sig.s.1.1*sig.lambda.lamda   )  )
  
  Sigma2.z.t = sigma.sq.t * val
  mu.z.t = val * (   t(sig.s.1)%*%(Yt-mu.t)%*%sig.lambda   )
  
  out = mu.z.t + abs(rnorm(1,0,Sigma2.z.t))
  return(out)
}
#zt.sample(1, Sigma.S, lambda, Sigma.I, Y.list[[1]], mu.star.list[[1]])

#=====================================================================

# sigmat^2

sigma.sq.sample = function(a, P, Yt, mu.t, zt, lambda, Sigma.S, Sigma.I, Sinv, Iinv){
  
  # Sinv = solve(Sigma.S)
  # Iinv = solve(Sigma.I)
  
  n = nrow(Sigma.S)
  lambda = as.vector(lambda)
  
  lambda.1 = rep(1,n) %x% lambda
  
  term1 = Yt - mu.t - abs(zt)*lambda.1
  
  p1 = (a+n*P+1)/2
  p2 = 0.5*(a + t(term1) %*% (Sinv %x% Iinv) %*% (term1) + zt^2)
  
  out = rinvgamma(1, p1, p2)
  return(out)
}
#sigmat2.sample(a,P, Yt = as.vector(Y.list[[1]]), mu.t = as.vector(mu.star.list[[1]]), zt = z[1], lambda, Sigma.S, Sigma.I)

#===========================================================

# a

a.sample = function(sigma.sq){
  
  x.vals <- seq(0.1, 20, 0.1)
  prob.vals.each <- numeric(length = length(x.vals))
  
  for (i.prob in 1:length(x.vals)){
    a.star <- x.vals[i.prob]
    prob.vals.each[i.prob] <- sum( dinvgamma(sigma.sq, a.star/2, a.star/2, log = TRUE) )  
    #prob.vals.each[i.prob] <- prod(dinvgamma(sigma.sq, a.star/2, a.star/2))  
  }
  
  prob.vals <- exp(prob.vals.each - max(prob.vals.each)) / sum(exp(prob.vals.each - max(prob.vals.each)))
  
  out = sample(x.vals, 1, prob = prob.vals)
  return(out)
}
#a.sample(1)

#===========================================================

# Sigma.I

Sigma.I.sample <- function(Y.list, T, Sigma.B, Sigma.S, Sinv, Binv, P, B.beta, mu.beta, mu.star.list, z, lambda){
  
  L <- nrow(Sigma.B)
  n <- nrow(Sigma.S)
  lambda = as.vector(lambda)
  
  # Sinv = solve(Sigma.S)
  # Binv = solve(Sigma.B)
  
  nu.I <- 0.01 + n * T + n * L 
  
  foo <- diag(0, P)
  
  for (i.t in 1:T){
    term = Y.list[[i.t]] - mu.star.list[[i.t]] - abs(z[i.t]) * rep(1,n) %x% t(lambda)
    foo <- foo + t(term) %*% Sinv %*% ( term ) 
  }
  
  term1 = B.beta - rep(1, n*L) %x% t(mu.beta)
  choo <- t(term1) %*% (Sinv %x% Binv) %*% ( term1 )
  
  Psi.I <- 0.01 * diag(P) + foo + choo
  
  ret <- rinvwishart(nu.I, Psi.I)
  
  return(ret)
  
}
#Sigma.I.sample(Y.list, T, Sigma.B, Sigma.S, P, B, mu.beta, mu.star.list, z, lambda)

#===========================================================

# Sigma.B

Sigma.B.sample <- function(Sigma.S, Sigma.I, Sinv, Iinv, mu.beta, B.beta.star, L){
  
  n <- nrow(Sigma.S)
  P <- nrow(Sigma.I)
  
  # Sinv = solve(Sigma.S)
  # Iinv = solve(Sigma.I)
  
  nu.B <- 0.01 + n*P
  
  term = B.beta.star - t( rep(1, L) ) %x% rep(1,n) %x% mu.beta
  
  Psi.B <- 0.01 * diag(L) +  t(term) %*% (Sinv %x% Iinv) %*% ( term )
  
  ret <- rinvwishart(nu.B, Psi.B)
  
  return(ret)
  
}
#Sigma.B.sample(Sigma.S, Sigma.I, mu.beta, B.beta.star, L)

#===========================================================

# rho, nu, gamma

#in the main file

