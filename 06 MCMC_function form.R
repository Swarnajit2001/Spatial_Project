rm(list = ls())

load("03 We_Did_It.Rdata")

#-------------------------------------------------
#libraries

library(fdrtool)
library(rSPDE)
library(Rcpp)
library(fda)
library(mvnfast)
library(gstat)
library(sp)


#--------------------------------------------------
#all the posterior functions from 05 Posterior.R
source("05 Posterior.R")

T <- nrow(final.data[[1]])
P <- length(final.data)
L <- 5

#--------------------------------------------------
#                 The MCMC chain                 #
#--------------------------------------------------

MCMC.chain.for.a.region <- function(R, B){
  
  #-------------------------------------------------
  #specifications
  
  S.indv = which(S.reg[,3] == R) 
  
  reg.indv <- S.reg[S.indv, ]
  
  N = length(S.indv)
  
  #--------------------------------------------------
  #preparing X
  
  basis <-  create.bspline.basis(c(1901,2022), nbasis = L)
  B.mat <- eval.basis(1901:2022, basis)
  X <- t(B.mat)
  
  ##--------------------------------------------------
  #Sigma.S construction
  
  Sigma.S.create = function(rho, nu, gamma, dist){
    
    mat = matrix(as.vector(matern.covariance(dist, 1/rho, nu, 1)), ncol = N)
    
    out = gamma * mat + (1-gamma) * diag(N)
    return(out) 
  }
  
  #--------------------------------------------------
  #data for the region and Y for beta initiation
  
  indv.data = list()
  Y.for.init.beta <- numeric(length = N * T * P)
  
  for(i.p in 1:P){
    
    foo <- final.data[[i.p]][,S.indv]
    
    indv.data[[i.p]] <- foo
    
    Y.for.init.beta[ ((i.p - 1) * (N * T) + 1) :  (i.p * N * T) ] <- as.vector( t(foo) )
    
  }
  
  #--------------------------------------------------
  #making the data time-wise 
  
  Y.list = list()
  
  for(i in 1:T){
    
    foo = matrix(0, nrow = N, ncol = P)
    for(j in 1:P){
      foo[,j] = indv.data[[j]][i, ]
    }
    
    Y.list[[i]] = foo
    
  }
  
  #--------------------------------------------------
  #distances needed to create Sigma.S
  
  n.indv <- nrow(reg.indv)
  combinations <- expand.grid(1:n.indv, 1:n.indv)
  new_matrix <- matrix(0, nrow = n.indv*n.indv, ncol = 4)
  
  new_matrix[, 1:2] <- as.matrix(reg.indv[combinations[, 1], 1:2 ])
  new_matrix[, 3:4] <- as.matrix(reg.indv[combinations[, 2], 1:2 ])
  
  dist = sqrt((new_matrix[,1]-new_matrix[,3])^2 + (new_matrix[,2]-new_matrix[,4])^2)
  
  D <- max(dist) - min(dist)
  
  #--------------------------------------------------
  ##               MCMC Chain                      ##
  #--------------------------------------------------
  
  #-- -- --
  #delete it
  # B = 200
  #-- -- --
  
  #-----
  #index for B.beta in mu_beta.sample
  ind.B.beta <- numeric()
  
  foo <- numeric()
  for ( i in 1:N){
    foo <- c( foo, ((i-1) * P * L + 1) : ( (i-1) * P * L + L   ) )
  }
  
  for (i in 1:P){
    ind.B.beta <- c(ind.B.beta, (foo + (i-1)*L) ) 
  }
  
  #-----
  #index for B.beta.star in Sigma.B.sample
  
  foo <- numeric()
  for (i in 1:P){
    foo <- c(foo, ( L *(i-1) + 1 ))
  }
  
  foo1 <- numeric()
  for (i in 1:N){
     foo1 <- c(foo1, (foo + (i-1) * L * P ) )
  }
  
  ind.B.beta.star <- numeric()
  for (i in 1:L){
    ind.B.beta.star <- c(ind.B.beta.star, (  foo1 + (i-1))  )
  }
  
  #-----
  #index for mu.beta initiation
  ind.mu.init <- rep(rep(1:P, each = L), N)
  
  #-----
  #storage
  store.beta <- matrix(0, nrow = B, ncol = (N*L*P))
  store.mu.beta <- matrix(0, nrow = B, ncol = P)
  store.lambda <- matrix(0, nrow = B, ncol = P)
  store.zt <- matrix(0, nrow = B, ncol = T)
  store.sigma.sq <- matrix(0, nrow = B, ncol = T)
  store.a <- numeric(length = B)
  store.sigma.I <- array(0, dim = c(P,P,B))
  store.sigma.B <- array(0, dim = c(L,L,B))
  store.rho <- numeric(length = B)
  store.nu <- numeric(length = B)
  store.gamma <- numeric(length = B)
  store.sigma.S <- array(0, dim = c(N,N,B))
  
  #-----
  #initiations
  beta <- as.vector(( diag(N*P) %x% (solve( X %*% t(X) ) %*% X) ) %*% Y.for.init.beta)
  
  #list making for mu.t.star in various other functions
  mu.star.list <- list()
  for (i.t in 1:T){
    x.t <- matrix(X[,i.t], ncol = 1) 
    X.t <- diag(N*P) %x% x.t
    
    mu.t <- t( X.t ) %*% as.vector(beta)
    mu.t.star <- t(matrix(mu.t, ncol = N))
    
    mu.star.list[[i.t]] <- mu.t.star
  }
  
  mu.beta <- numeric(length = P)
  for (i.mu.ini in 1:P)
  {
    mu.beta[i.mu.ini] <- mean(beta[which(ind.mu.init == i.mu.ini)])
  }
  
  lambda <- rep(0, P)
  z <- rep(0, T)
  a <- 10
  sigma.sq <- rinvgamma(T, a/2, a/2)

  mat.Y.for.init.Sigma.S <- matrix(Y.for.init.beta, ncol = P)
  Sigma.I <- var(mat.Y.for.init.Sigma.S)
  
  mat.Y.for.init.Sigma.B <- matrix( beta[ind.B.beta.star], ncol = L)
  Sigma.B <- var(mat.Y.for.init.Sigma.B)
  
  rho.vec <- numeric(length = P)
  for (i.S.S.ini in 1:P){
    foo <- indv.data[[i.S.S.ini]]
    foo.vec <- colMeans(foo)
    
    df <- data.frame(x = reg.indv[ , 1], y = reg.indv[ , 2], values = foo.vec)
    coordinates(df) <- ~x + y
    vv <- variogram(values~1, df)
    m.fit <- fit.variogram(vv, vgm(1,"Exp"))
    rho.vec[i.S.S.ini] <- m.fit$range
    
  }
  rho <- mean(rho.vec)
  
  nu <- 0.5
  gamma <- 0.5
  Sigma.S <- Sigma.S.create(rho, nu, gamma, dist)
  
  #-----
  #initiating the inverses
  Binv = solve(Sigma.B)
  Sinv = solve(Sigma.S)
  Iinv = solve(Sigma.I)
  
  #-- -- --
  #delete it
  # i.B <- 1
  #-- -- --
  
  #-----
  #loops for iteration starts here
  
  for (i.B in 1:B)
  {
    
    set.seed(221453 + i.B)
    
    burgers = proc.time()[3]
    
    #-----
    #beta
    beta <- beta.sample(X, Sigma.B, Sigma.I, Sigma.S, Binv, Sinv, Iinv, sigma.sq, Y.list, z, lambda, mu.beta, T)
    store.beta[i.B,] <- beta
    
    #-----
    #mu.beta
    
    B.beta <- matrix( beta[ind.B.beta], ncol = P )
    
    mu.beta <- mu_beta.sample(Sigma.S, Sigma.B, Sigma.I, Binv, Sinv, Iinv, P, z, B.beta)
    
    store.mu.beta[i.B, ] <- mu.beta
    
    #------
    #lambda
    
    lambda <- lambda.sample(Sigma.S, z, sigma.sq, Sigma.I, Sinv, Iinv, T, Y.list, mu.star.list)
    
    store.lambda[i.B, ] <- lambda
    
    #-----
    #|zt|
    
    z <- numeric(length = T)
    
    for (i.t in 1:T){
      sigma.sq.t <- sigma.sq[i.t]
      Yt <- Y.list[[i.t]]
      mu.t <- mu.star.list[[i.t]]
      
      z[i.t] <- zt.sample(sigma.sq.t, Sigma.S, Sinv, Iinv, lambda, Sigma.I, Yt, mu.t)
    }
    
    store.zt[i.B, ] <- z
    
    #-----
    #sigma2.t
    
    # sigma.sq <- numeric(length = T)
    # for (i.t in 1:T){
    #   Yt <- as.vector(t(Y.list[[i.t]]))
    #   mu.t <- as.vector(t(mu.star.list[[i.t]]))
    #   zt <- z[i.t]
    #   sigma.sq[i.t] <- sigma.sq.smple(a, P, Yt, mu.t, zt, lambda, Sigma.S, Sigma.I, Sinv, Iinv)
    # }
    
    #sigma.sq <- sigma.sq + abs(rnorm(T))
    
    store.sigma.sq[i.B, ] <- sigma.sq

    #-----
    #a

    a <- a.sample(sigma.sq)
    store.a[i.B] <- a
    
    #-----
    #Sigma.I
    
    Sigma.I <- Sigma.I.sample(Y.list, T, Sigma.B, Sigma.S, Sinv, Binv, P, B.beta, mu.beta, mu.star.list, z, lambda)
    
    Iinv = solve(Sigma.I)
    store.sigma.I[ , , i.B] <- Sigma.I
    
    #-----
    #Sigma.B
    
    B.beta.star <- matrix( beta[ind.B.beta.star], ncol = L)
    
    Sigma.B <- Sigma.B.sample(Sigma.S, Sigma.I, Sinv, Iinv, mu.beta, B.beta.star, L)
    
    Binv = solve(Sigma.B)
    store.sigma.B[ , , i.B] <- Sigma.B
    
    #-----
    #rho, nu, gamma
    
    rho.cand <- exp(rnorm(1, log(rho), 1))
    nu.cand <- exp(rnorm(1, log(nu), 1))
    gamma.cand <- invlogit(rnorm(1, logit(gamma), 1))

    Sigma.S.cand <- Sigma.S.create(rho.cand, nu.cand, gamma.cand, dist)

    #the acceptance ratio
    term1 <- 1

    for (i.t in 1:T)
    {

      Yt <- as.vector(t(Y.list[[i.t]]))
      mu.t <- as.vector(t(mu.star.list[[i.t]]))
      z.t <- z[i.t]
      sigma.sq.t <- sigma.sq[i.t]

      mu.rng <- mu.t + abs(z.t) * rep(1,N) %x% lambda
      sigma.rng <- sigma.sq.t * Sigma.S %x% Sigma.I
      sigma.rng.cand <- sigma.sq.t * Sigma.S.cand %x% Sigma.I

      foo <- t(as.matrix(Yt))

      #term <- term1 * exp(dmvnrm_arma_mc(foo, mu.rng, sigma.rng, logd = TRUE) - dmvnrm_arma_mc(foo, mu.rng, sigma.rng.cand, logd = TRUE))
      term1 <- term1 * exp(mvnfast::dmvn(foo, mu.rng, sigma.rng, log = TRUE, isChol = TRUE) - mvnfast::dmvn(foo, mu.rng, sigma.rng.cand, log = TRUE, isChol = TRUE))
      #term1 <- term1 * exp(dMvn(Yt, mu.rng, sigma.rng, log = TRUE) - dMvn(Yt, mu.rng, sigma.rng.cand, log = TRUE))

    }

    term2 <- dnorm(log(nu.cand), -1.2, 1) / dnorm(log(nu), -1.2, 1)

    term3 <- rho.cand * (D - rho.cand) / ( rho * (D - rho) )

    R <- min(1, (term1 * term2 * term3))

    if (runif(1) < R){
      rho <- rho.cand
      nu <- nu.cand
      gamma <- gamma.cand

      #update Sigma.S also
      Sigma.S <- Sigma.S.cand
    }

    store.rho[i.B] <- rho
    store.nu[i.B] <- nu
    store.gamma[i.B] <- gamma
    store.sigma.S[ , , i.B] <- Sigma.S

    Sinv = solve(Sigma.S)
    
    alisson <- proc.time()[3]
    
    aladeen <-  alisson - burgers
    print(aladeen)
    
    #-----
    #aesthetics
    print(paste0("We are at iteration ", i.B))
    
  }
  
  #-----
  #output
  list.mcmc <- list(  beta = store.beta,
                      mu.beta = store.mu.beta,
                      lambda = store.lambda, 
                      z = store.zt, 
                      sigma.sq = store.sigma.sq, 
                      a = store.a, 
                      Sigma.I = store.sigma.I, 
                      Sigma.B = store.sigma.B, 
                      rho = store.rho, 
                      nu = store.nu, 
                      gamma = store.gamma,
                      Sigma.S = store.sigma.S )
  
  return(list.mcmc)
  
}

#=====================================================================

#-- -- --
#R = 34

# out <- MCMC.chain.for.a.region(34, 200)
# save(out, file = "final.MCMC.R34.Rdata")
# rm(out)

#-- -- --
#R = 4

out <- MCMC.chain.for.a.region(4, 200)
save(out, file = "final.MCMC.R04.Rdata")
rm(out)











