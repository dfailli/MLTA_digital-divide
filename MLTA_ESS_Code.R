## Load the libraries ####

library(readr)
library(igraph)
library(bipartite)
library(ggplot2) 
library(MASS)
library(foreach)
library(mclust)
library(doParallel)
library(plotrix)
library(cds)
library(ggplot2)
library(dplyr)
library(likert)
library(maps)
library(magrittr)


## MLTA functions ####


AitkenAcc <- function(l, lA, iter){
  
  # Perform an Aitken acceleration: given that the EM algorithm converges linearly, 
  # an Aitken acceleration (McLachlan & Peel, 2000) can be used to speed up convergence
  
  a <- (l[3] - l[2]) / (l[2] - l[1])
  lA[1] <- lA[2]
  
  if (a > 0 & a < 1) 
  {
    lA[2] <- l[2] + (l[3] - l[2]) / (1 - a)
  }else{
    lA[2] <- l[3]
  }
  
  diff <- abs(lA[2] - lA[1])
  
  if ((l[3] < l[2]) & (iter > 5))
  {
    stop("Decrease in log-likelihood")
  }	
  
  Out <- c(diff, lA)
  Out
}



ResTable <- function(bicG, restype)
{
  
  # Display the model log-likelihood and BIC (for model selection) 
  
  if (restype == 'll') {
    resBIC <- vector('list', 1)
    names(resBIC) <- 'Table of LL (G-H Quadrature correction)'
    resBIC[[1]] <- bicG
  }
  
  if (restype == 'llva') {
    resBIC <- vector('list', 1)
    names(resBIC) <- 'Table of LL (variational approximation)'
    resBIC[[1]] <- bicG
  }
  
  if (restype == 'BIC') {
    resBIC <- vector('list', 2)
    names(resBIC) <- c('Table of BIC Results', 'Model Selection')
    resBIC[[1]] <- bicG
    resBIC[[2]] <-
      paste(colnames(bicG)[t(bicG == min(bicG)) * seq(1:ncol(bicG))], rownames(bicG)[(bicG ==
                                                                                        min(bicG)) * seq(1:nrow(bicG))], sep = ', ')
    names(resBIC[[2]]) <- 'Model with lower BIC:'
  }
  
  resBIC
}


f_lca <-function(X, DM, G, tol, maxiter, beta0){ 
  
  # Perform a Latent Class Analysis (when D=0, the MLTA reduces to the LCA model)
  # X=incidence matrix, DM=design matrix, G=number of components, tol=tolerance level
  # maxiter=maximum number of iterations, beta0=initial values for beta, if specified
  
  N <- nrow(X)
  M <- ncol(X) 
  J <- ncol(DM)
  
  # Initialize EM Algorithm
  
  if(!is.null(beta0)) beta = beta0 else{
    beta <- rep(0,J*(G-1))
    beta <- matrix(beta,ncol=G-1)
  }
  
  exb <- exp(DM %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  # priors
  
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N)   z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,]))
  
  # Set up
  
  p <- (t(z) %*% X) / colSums(z) # proportion of success in group g
  
  # Initialize Aitken acceleration
  
  l <- numeric(3)
  lA <- numeric(2)
  
  # Iterative process
  
  v <- matrix(0, N, G)
  W <- list()
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  tol <- 0.1 ^ 6
  print(c(beta))
  
  while(diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    beta.old=beta
    
    # M-step 
    
    lk = sum(z*log(eta))
    it = 0; lko = lk
    XXdis = array(0,c(G,(G-1)*ncol(DM),N))
    for(i in 1:N){
      XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(DM[i,]))
    }
    while((lk-lko>10^-6 & it<100) | it==0){
      it = it+1; lko = lk 
      sc = 0; Fi = 0
      for(i in 1:N){
        pdis = eta[i,]
        sc = sc+t(XXdis[,,i])%*%(z[i,]-pdis)
        Fi = Fi+t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
      }
      
      dbe = as.vector(ginv(Fi)%*%sc)
      mdbe = max(abs(dbe))
      if(mdbe>0.5) dbe = dbe/mdbe*0.5
      be0 = c(beta)
      flag = TRUE
      while(flag){
        beta = be0+dbe
        Eta = matrix(0,N,G)
        for(i in 1:N){
          
          if(ncol(DM)==1) Eta[i,] = XXdis[,,i]*beta
          else Eta[i,] = XXdis[,,i]%*%beta
        }	
        if(max(abs(Eta))>100){
          dbe = dbe/2
          flag = TRUE	
        }else{
          flag = FALSE
        }	        	
      }
      if(iter/10 == floor(iter/10))       print(beta)
      
      beta = matrix(beta, J, G-1)    
      exb <- exp(DM %*% beta) # update priors
      eta <- cbind(1,exb)/(rowSums(exb)+1)
      
      lk = sum(z*log(eta))
    }
    
    
    # E-step
    PS = array(apply(p, 1, function(xx) dbinom(c(X), 1, xx)), c(N, M, G))
    v <- eta *  apply(PS, c(1,3), prod)
    vsum <- rowSums(v)
    z <- v / vsum      # posterior
    ll <- rep(1, N) %*% log(vsum) # log-likelihood
    
    
    # Aitken acceleration
    
    l <- c(l[-1], ll)
    Out <- AitkenAcc(l, lA, iter)
    diff <- Out[1]
    lA <- Out[-1]
    
    print(iter)
    
  } # end EM
  
  p <- (t(z) %*% X) / colSums(z)	
  
  BIC <- -2*ll + (J*(G-1)+G*M) * log(N)
  
  expected <- vsum * N
  
  
  colnames(p) <- NULL
  rownames(p) <- rownames(p, do.NULL = FALSE, prefix = "Group ")
  colnames(p) <- colnames(p, do.NULL = FALSE, prefix = "p_g")
  colnames(beta) <- 2:G
  rownames(beta) <- c(colnames(DM))
  names(ll)<- "Log-Likelihood:"
  names(BIC)<-"BIC:"
  
  out <- list(p = p, eta = eta, LL = ll, BIC = BIC)
  list(p = p, eta = eta, LL = ll, BIC = BIC, z = z, expected = expected,
       beta = beta, rrr=exp(beta)) 
}



f_lca_nstarts <- function(X, DM, G, nstarts, tol, maxiter, beta0) {
  
  # Compute the LCA model for different starting values to avoid 
  # the algorithm from remaining trapped in local maxima solutions
  
  out <- f_lca(X, DM, G, tol, maxiter, beta0)
  
  for(i in 2:nstarts){
    out1 <- f_lca(X, DM, G, tol, maxiter, beta0)
    if(out1$LL > out$LL) out <- out1
  }
  out
}


lca <- function(X, DM, G, nstarts = 3, tol = 0.1^2, maxiter = 250, beta0=NULL) {
  
  # Estimation with LCA model
  
  if (any(G < 1)) {
    print("Specify G > 0!")
    return("Specify G > 0!")
  }
  
  if (any(G == 1)) {
    out <- f_lca_nstarts(X, DM, G, nstarts, tol, maxiter, beta0 = beta0)
  } else{
    if (length(G) == 1) {
      out <- f_lca_nstarts(X, DM, G, nstarts, tol, maxiter, beta0 = beta0)
    } else{
      out <- vector("list", length(G) + 1)
      names(out) <- c(paste('G', G, sep = '='), 'BIC')
      i <- 0
      for (g in G) {
        i <- i + 1
        out[[i]] <- f_lca_nstarts(X, DM, g, nstarts, tol, maxiter, beta0 = beta0)
      }
      out[[length(G) + 1]] <- tableBIC(out)
    }
  }
  out
}


buildYYh <- function(YY, mu, dimy, N){
  
  # Compute E[uu'], where u is the D-dimensional continuous latent trait
  # YY=C + mu^2, where C and mu are the covariance matrix and mean vector of u
  # N=number of sending nodes
  # dimy=dimension of the continuous latent trait
  
  YYh <- matrix(NA, (dimy + 1) * (dimy + 1), N)
  
  for(ddi in 1:dimy){
    YYh[((ddi - 1) * dimy + ddi):(ddi - 1 + (ddi * dimy)),] <- YY[((ddi - 1) * dimy + 1):(ddi * dimy),]
    YYh[ddi * dimy + ddi,] <- mu[ddi,]
  }
  
  YYh[(dimy * (dimy + 1) + 1):((dimy + 1) * (dimy + 1) - 1),] <- mu[1:ddi,]
  
  YYh[(dimy + 1) * (dimy + 1),] <- 1
  
  YYh		
}


f_mlta_nstarts <- function(X, DM, G, D, nstarts, tol, maxiter, pdGH, beta0)
{
  
  # Estimate MLTA model for different random starting values to avoid local maxima solutions
  # X=incidence matrix, DM=design matrix, G=number of components, 
  # D=latent trait dimension, nstart=number of random starting values
  # tol=tolerance level, maxiter=maximum number of iterations,
  # pdGH=number of quadrature points, beta0=initial beta if specified
  
  out <- f_mlta(X, DM, G, D, tol, maxiter, pdGH, beta0)
  
  if(nstarts > 1){ 
    for(i in 2:nstarts){
      out1 <- f_mlta(X, DM, G, D, tol, maxiter, pdGH, beta0)
      if(out1$LL > out$LL) out <- out1
    }
  }
  return(out)
}


f_mlta_wfix <- function(X, DM, G, D, tol, maxiter, pdGH, beta0)
{
  
  # Perform the contrained MLTA model (w fixed with respect to g)
  # X=incidence matrix, DM=design matrix, G=number of components, 
  # D=latent trait dimension,
  # tol=tolerance level, maxiter=maximum number of iterations,
  # pdGH=number of quadrature points, beta0=initial beta if specified
  
  N <- nrow(X)
  M <- ncol(X) 
  J <- ncol(DM)
  
  # Initialize EM Algorithm
  
  if(!is.null(beta0)) beta = beta0 else{
    beta <- rep(0,J*(G-1))
    beta <- matrix(beta,ncol=G-1)
  }
  
  exb <- exp(DM %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  # priors
  
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N)   z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,]))
  
  # Set up
  
  p <- (t(z) %*% X) / colSums(z)
  
  
  ###
  
  xi <- array(20, c(N, M, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  w <- matrix(rnorm(M * D), D, M)
  b <- matrix(rnorm(M, G), G, M)
  
  
  wh <- matrix(0, D + G, M)
  
  C <- array(0, c(D * D, N, G))
  
  mu <- array(0, c(N, D, G))
  
  YY <- array(0, c(D * D, N, G))
  
  gam <- matrix(0, D + G, M)
  K <- array(diag(D + G), c(D + G, D + G, M))
  
  lxi <- matrix(0, G, N)
  
  # Iterative process
  
  v <- matrix(0, N, G)
  W <- list()
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  tol <- 0.1 ^ 4
  print(c(beta))
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    beta.old=beta
    ll_old <- ll
    
    if (D == 1) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          1 / (1 - 2 * rowSums(sweep(lambda_xi[, , g], MARGIN = 2, w ^ 2, `*`)))
        mu[, , g] <-
          C[, , g] * rowSums(sweep(
            X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`),
            MARGIN = 2, w, `*`))
        
        YY[, , g] <- matrix(C[, , g] + mu[, , g] ^ 2, ncol = 1)
        
        xi[, , g] <- YY[, , g] %*% w^2 + mu[, , g] %*% (2 * b[g, ] * w) + 
          matrix(b[g, ]^2, nrow = N, ncol = M, byrow = TRUE)
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[1, ] <- colSums(crossprod(z * mu[, 1, ], X - 0.5))
      gam[2:(1 + G), ] <- crossprod(z, X - 0.5)
      
      K[1, 1, ] <- apply(aperm(array(z * YY[1, , ], c(N, G, M)), c(1, 3, 2)) * lambda_xi, 2, sum)
      
      K[2:(G + 1), 1, ] <- t(apply(aperm(array(
        z * mu[, 1, ], c(N, G, M)
      ), c(1, 3, 2)) * lambda_xi, c(2, 3), sum))
      K[1, 2:(G + 1), ] <- K[2:(G + 1), 1, ]
      
      for (m in 1:M)
      {
        for (g in 1:G)
          K[1 + g, 1 + g, m] <- sum(z[, g] * lambda_xi[, m, g])
        
        wh[, m] <- -ginv(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      w <- as.matrix(t(w))
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        lxi[g, ] <-
          0.5 * log(C[1, , g]) + mu[, 1, g] ^ 2 / (2 * C[1, , g]) + rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g] ^ 2 + 
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
              sweep(X - 0.5, MARGIN = 2, b[g, ], `*`))
      } # end for (g in 1:G)
    }
    
    if (D == 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], 
                        tcrossprod(w, X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`))), 
                  2, function(x) matrix(x[1:4], nrow = D) %*% x[-(1:4)]))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums(crossprod(w, matrix(x, ncol = D)) * t(w))) + 
              tcrossprod(2 * b[g, ] * t(w), mu[, , g]) + matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = N,
                byrow = FALSE
              )
          )
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[3:(2 + G), ] <- crossprod(z, X - 0.5)
      
      aa <- aperm(array(z, c(N, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(z, c(N, G, D * D)), c(3, 1, 2)) * YY
      
      K[3:(2 + G), 3:(2 + G), ] <-
        apply(apply(aperm(array(z, c(
          N, G, M
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      for (g in 1:G) {
        K[2 + g, 1:2, ] <- crossprod(aa[, , g], lambda_xi[, , g])
        K[1:2, 2 + g, ] <- K[2 + g, (1:2), ]
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
      }
      
      K[1:2, 1:2, ] <- kk
      
      for (m in 1:M)
      {
        gam[1:D, m] <- apply(aa * (X[, m] - 0.5), 2, sum)
        wh[, m] <- -ginv(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        # Approximation of log(p(y|z))
        
        detC <- C[1, , g] * C[4, , g] - C[3, , g] * C[2, , g]
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C[4, , g] / detC, -C[2, , g] / detC, -C[3, , g] /
                                                detC, C[1, , g] / detC, t(mu[, , g])), 2, function(x)
                                                  t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)])) + 
          rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g]^2 + 
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
              sweep(X - 0.5, MARGIN = 2, b[g, ], `*`))
      } # end for (g in 1:G)
      
    }# end D=2
    
    if (D > 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], w %*% t(
            X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)
          )), 2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% x[-(1:(D * D))]))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums((t(w) %*% matrix(x, ncol = D)) * t(w))) + 
              (2 * b[g, ] * t(w)) %*% t(mu[, , g]) + 
              matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = N,
                byrow = FALSE
              )
          )
        
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[(D + 1):(D + G), ] <- t(z) %*% (X - 0.5)
      
      aa <- aperm(array(z, c(N, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(z, c(N, G, D * D)), c(3, 1, 2)) * YY
      
      K[(D + 1):(D + G), (D + 1):(D + G), ] <-
        apply(apply(aperm(array(z, c(
          N, G, M
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      for (g in 1:G) {
        K[D + g, 1:D, ] <- crossprod(aa[, , g], lambda_xi[, , g])
        K[1:D, D + g, ] <- K[D + g, 1:D, ]
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
      }
      
      K[1:D, 1:D, ] <- kk
      
      for (m in 1:M)	{
        gam[1:D, m] <- apply(aa * (X[, m] - 0.5), 2, sum)
        wh[, m] <- -ginv(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        # Approximation of log(p(y|z))
        
        detC <- apply(C[, , g], 2, function(x)
          det(matrix(x, D, D)))
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C[, , g], t(mu[, , g])), 2, function(x)
            t((x[-(1:(D * D))])) %*% ginv(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(
              1:(D * D))])) + rowSums(
                log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] *
                  xi[, , g]^2 + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
                  sweep(X - 0.5, MARGIN = 2, b[g, ], `*`)
              )
      } # end for (g in 1:G)
      
    } # end if D > 2
    
    # M-step 
    
    lk = sum(z*log(eta))
    it = 0; lko = lk
    XXdis = array(0,c(G,(G-1)*ncol(DM),N))
    for(i in 1:N){
      XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(DM[i,]))
    }
    while((lk-lko>10^-6 & it<100) | it==0){
      it = it+1; lko = lk 
      sc = 0; Fi = 0
      for(i in 1:N){
        pdis = eta[i,]
        sc = sc+t(XXdis[,,i])%*%(z[i,]-pdis)
        Fi = Fi+t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
      }
      
      dbe = as.vector(ginv(Fi)%*%sc)
      mdbe = max(abs(dbe))
      if(mdbe>0.5) dbe = dbe/mdbe*0.5
      be0 = c(beta)
      flag = TRUE
      while(flag){
        beta = be0+dbe
        Eta = matrix(0,N,G)
        for(i in 1:N){
          
          if(ncol(DM)==1) Eta[i,] = XXdis[,,i]*beta
          else Eta[i,] = XXdis[,,i]%*%beta
        }	
        if(max(abs(Eta))>100){
          dbe = dbe/2
          flag = TRUE	
        }else{
          flag = FALSE
        }	        	
      }
      if(iter/10 == floor(iter/10))       print(beta)
      
      beta = matrix(beta, J, G-1)    
      exb <- exp(DM %*% beta) # updfe priors
      eta <- cbind(1,exb)/(rowSums(exb)+1)
      
      lk = sum(z*log(eta))
    }
    
    
    # E-step
    v <- eta * t(exp(lxi))
    if(any(is.nan(v))) browser()
    vsum <- apply(v, 1, sum)
    z <- v / vsum
    ll <- sum(rep(1, N) * log(vsum))
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
  } # end while(diff>tol)
  
  # Correction to the log-likelihood 
  
  # Gauss-Hermite Quadrature
  
  npoints <- round(pdGH ^ (1 / D))
  ny <- npoints ^ D
  GaHer <- glmmML::ghq(npoints, FALSE)
  Ygh <- expand.grid(rep(list(GaHer$zeros), D))
  Ygh <- as.matrix(Ygh)
  Wgh <-
    apply(as.matrix(expand.grid(rep(
      list(GaHer$weights), D
    ))), 1, prod) * apply(exp(Ygh ^ 2), 1, prod)
  
  Hy <- apply(Ygh, 1, mvtnorm::dmvnorm)
  Beta <- Hy * Wgh / sum(Hy * Wgh)
  
  fxy <- array(0, c(N, ny, G))
  
  for (g in 1:G)
  {
    Agh <- t(tcrossprod(t(w), Ygh) + b[g, ])
    pgh <- 1 / (1 + exp(-Agh))
    fxy[, , g] <-
      exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh)))
  }
  
  
  LLva <- ll
  BICva <- -2 * LLva + {
    G * M + M * D - D * {
      D - 1
    } / 2 + J*(G - 1)
  } * log(N)
  
  llGH1 <- apply(rep(Beta, each = N) * fxy, c(1, 3), sum)
  LL <- sum(log(rowSums(eta * (llGH1))))
  
  BIC <- -2 * LL + {
    G * M + M * D - D * {
      D - 1
    } / 2 + J*(G - 1)
  } * log(N)
  
  expected <- colSums(eta * (llGH1)) * N
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- rownames(b, do.NULL = FALSE, prefix = "Group ")
  
  rownames(w) <- rownames(w, do.NULL = FALSE, prefix = "Dim ")
  colnames(w) <- colnames(w, do.NULL = FALSE, prefix = "Item ")
  
  colnames(beta) <- 2:G
  rownames(beta) <- c(colnames(DM))
  
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(BICva) <- c("BIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(BIC) <- c("BIC (G-H Quadrature correction):")
  
  out1 <-
    list(
      b = b,
      w = w,
      eta = eta,
      LL = LL,
      BIC = BIC,
      LLva = LLva,
      BICva = BICva,
      expected = expected,
      mu = mu,
      C = C,
      z = z,
      beta = beta, 
      rrr=exp(beta)
    )
  
  out1
  
}


f_mlta_nstarts_wfix <- function(X, DM, G, D, nstarts, tol, maxiter, pdGH, beta0)
{
  
  # Perform the constrained MLTA model for different starting values to 
  # avoid local maxima solutions
  # X=incidence matrix, DM=design matrix, G=number of components, 
  # D=latent trait dimension, nstart=number of random starting values
  # tol=tolerance level, maxiter=maximum number of iterations,
  # pdGH=number of quadrature points, beta0=initial beta if specified
  
  out <- f_mlta_wfix(X, DM, G, D, tol, maxiter, pdGH, beta0)
  
  if(nstarts > 1){ 
    for(i in 2:nstarts){
      out1 <- f_mlta_wfix(X, DM, G, D, tol, maxiter, pdGH, beta0)
      if(out1$LL > out$LL) out <- out1
    }
  }
  return(out)
}



f_mlta_methods <- function(X, DM, G, D, nstarts, tol, maxiter, pdGH, wfix, beta0=NULL)
{
  
  # Decide which model to estimate according to G (number of components), and
  # D (latent trait dimension)
  
  # X=incidence matrix, DM=design matrix, nstarts=number of random starting values
  # tol=tolerance level, maxiter=maximum number of iterations,
  # pdGH=number of quadrature points, 
  # wfix= indicator for the constrained model, beta0=initial beta if specified
  
  if (D == 0) {
    if (any(G == 1)) {
      out <- f_lca_nstarts(X, DM, G, nstarts, tol, maxiter, beta0=beta0)
    } else{
      if (length(G) == 1) {
        out <- f_lca_nstarts(X, DM, G, nstarts, tol, maxiter, beta0=beta0)
      } else{
        out <- vector("list", length(G) + 1)
        names(out) <- c(paste('G', G, sep = '='), 'BIC')
        i <- 0
        for (g in G) {
          i <- i + 1
          out[[i]] <- f_lca_nstarts(X, DM, g, nstarts, tol, maxiter, beta0=beta0)
        }
        out[[length(G) + 1]] <- tableBIC(out)
      }
    }
  } else {
    if (D > 0 && G == 1){
      if(length(D) == 1){ 
        out <- f_lta_nstarts(X, D, nstarts, tol, maxiter, pdGH)
        out$eta <- 1
      }else{
        out<-vector("list", length(D) + 1)
        names(out) <- c(paste('Dim y', D, sep = '='), 'BIC')
        i<-0
        for(diy in D){
          i <- i + 1
          out[[i]] <- f_lta_nstarts(X, diy, nstarts, tol, maxiter, pdGH)
          out[[i]]$eta <- 1
        }
        
        cat('BIC results',"\n \n")
        out[[length(D) + 1]]<-tableBIC(out)
      }
    }
    if (D > 0 && G > 1) {
      if (wfix == TRUE) {
        
        out <- f_mlta_nstarts_wfix(X, DM, G, D, nstarts, tol, maxiter, pdGH, beta0=beta0)
        class(out) <- c("mlta")
        
      } else{
        out <- f_mlta_nstarts(X, DM, G, D, nstarts, tol, maxiter, pdGH, beta0=NULL)
      }
    }
  }
  class(out) <- c("mlta")
  return(out)
}



ghLLlta <- function(X, w, b, D, pdGH) {
  
  # Gauss Hermite quadrature for Latent Trait Analysis
  
  b <- as.vector(b)
  N <- nrow(X)
  M <- ncol(X)
  npoints <- round(pdGH ^ (1 / D))
  ny <- npoints^D
  GaHer <- glmmML::ghq(npoints, modified = FALSE)
  ygh <- GaHer$zeros
  wgh <- GaHer$weights
  Ygh <- as.matrix(expand.grid(rep(list(ygh), D)))
  Wghm <- as.matrix(expand.grid(rep(list(wgh), D)))
  Wgh <- apply(Wghm, 1, prod) * apply(exp(Ygh ^ 2), 1, prod)
  Hy <- apply(Ygh, 1, mvtnorm::dmvnorm)
  Beta <- Hy * Wgh / sum(Hy * Wgh)
  Agh <- t(t(w) %*% t(Ygh) + b)
  pgh <- 1 / (1 + exp(-Agh))
  
  fxy <- exp(X %*% t(log(pgh)) + (1 - X) %*% t(log(1 - pgh)))
  
  llGH1 <- rowSums(rep(Beta, each = N) * fxy)
  llGH <- sum(log(llGH1))
  expected <- llGH1 * N
  
  list(expected = expected, llGH = llGH)
}


f_lta <- function(X, D, tol, maxiter, pdGH) {
  
  # Latent Trait Analysis (the MLTA coincides with the LTA when G=1)
  
  N <- nrow(X)
  M <- ncol(X) 
  
  # Initialization
  
  xi <- matrix(20, N, M)
  w <- matrix(rnorm(M * D), D, M)
  b <- rnorm(M)
  
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  diff <- 1
  iter <- 0L
  ll <- -Inf
  
  while (diff > tol & iter < maxiter) {
    iter <- iter + 1
    ll_old <- ll
    
    if (D == 1) {
      # # STEP 1: Computing the Latent Posterior Statistics
      
      C <-
        1 / (1 - 2 * rowSums(sweep(lambda_xi, MARGIN = 2, w ^ 2, `*`)))
      mu <-
        C * rowSums(sweep(
          X - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`),
          MARGIN = 2,
          w,
          `*`
        ))
      
      # # STEP 2: Optimising the Variational Parameters (xi)
      
      YY <- matrix(C + mu ^ 2, ncol = 1)
      mu <- matrix(mu, N, D)
      
      xi <-
        YY %*% w ^ 2 + mu %*% (2 * b * w) + matrix(b ^ 2,
                                                   nrow = N,
                                                   ncol = M,
                                                   byrow = TRUE)
      
      xi <- sqrt(xi)
      
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      YYh <- 2 * cbind(c(YY), mu, mu, 1)
      
      den <- t(YYh) %*% lambda_xi
      num <- rbind(c(mu), 1) %*% ((X - 0.5))
      
      wh <-
        apply(rbind(den, num), 2, function(x)
          - crossprod(ginv(matrix(x[1:4], 2, 2)), x[5:6]))
      
      w <- as.matrix(t(wh[1:D, ]))
      b <- wh[D + 1, ]
      
      # Approximation of Likelihood
      
      lxi <-  (0.5 * log(C) + c(mu) ^ 2 / (2 * C) +
                 rowSums(
                   log(sigma_xi) - 0.5 * xi - lambda_xi * (xi ^ 2) +
                     sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`) +
                     sweep(X - 0.5, MARGIN = 2, b, `*`)
                 ))
      
      C <- array(C, c(D, D, N))
      
      ll <- sum(lxi)
      bw <- t(rbind(b, w))
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
    } # end if D=1
    
    if (D == 2) {
      # # STEP 1: Computing the Latent Posterior Statistics
      
      C <- apply(lambda_xi, 1, function(x)
        solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
      
      mu <- apply(rbind(C, w %*% t(
        X - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`)
      )),
      2,
      function(x)
        matrix(x[1:4], nrow = D) %*% (x[-(1:4)]))
      
      # # STEP 2: Optimising the Variational Parameters (xi)
      
      YY <- C + apply(mu, 2, tcrossprod)
      
      xi <-
        apply(YY, 2, function(x)
          rowSums((t(w) %*% matrix(x, ncol = D)) * t(w))) +
        (2 * b * t(w)) %*% mu +
        matrix(b ^ 2,
               nrow = M,
               ncol = N,
               byrow = FALSE)
      
      xi <- t(sqrt(xi))
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      YYh <-
        2 * cbind(YY[1, ], YY[2, ], mu[1, ], YY[3, ], YY[4, ], mu[2, ], mu[1, ], mu[2, ], 1)
      
      den <- t(YYh) %*% lambda_xi
      num <- rbind(mu, 1) %*% ((X - 0.5))
      wh <-
        apply(rbind(den, num), 2, function(x)
          - crossprod(ginv(matrix(x[1:9], 3, 3)), x[10:12]))
      
      w <- wh[1:D, ]
      b <- wh[D + 1, ]
      
      # Approximation of Likelihood
      
      detC <- C[1, ] * C[4, ] - C[3, ] * C[2, ]
      
      lxi <- (
        0.5 * log(detC) +
          0.5 * apply(rbind(C[4, ] / detC,-C[2, ] / detC,-C[3, ] / detC, C[1, ] / detC, mu), 2,
                      function(x)
                        t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)])) +
          rowSums(
            log(sigma_xi) - 0.5 * xi - lambda_xi * (xi ^ 2) +
              sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`) + sweep(X - 0.5, MARGIN = 2, b, `*`)
          )
      )
      
      C <- array(C, c(D, D, N))
      mu <- t(mu)
      
      ll <- sum(lxi)
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    }# end if D=2
    
    if (D > 2) {
      # # STEP 1: Computing the Latent Posterior Statistics
      
      C <-
        apply(lambda_xi, 1, function(x)
          solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
      
      mu <-
        apply(rbind(C, w %*% t(
          X - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`)
        )) ,
        2, function(x)
          matrix(x[1:(D * D)], nrow = D) %*% (x[-(1:(D * D))]))
      
      # # STEP 2: Optimising the Variational Parameters (xi)
      
      YY <- C + apply(mu, 2, tcrossprod)
      
      xi <-
        t(
          apply(YY, 2, function(x)
            rowSums((
              t(w) %*% matrix(x, ncol = D)
            ) * t(w))) +
            (2 * b * t(w)) %*% mu + matrix(
              b ^ 2,
              nrow = M,
              ncol = N,
              byrow = FALSE
            )
        )
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # # STEP 3: Optimising the Model Parameters (w and b)
      
      YYh <- sweep(buildYYh(YY, mu, D, N), MARGIN = 2, 2, `*`)
      
      den <- YYh %*% lambda_xi
      
      num <- rbind(mu, 1) %*% ((X - 0.5))
      
      wh <-
        apply(rbind(den, num), 2, function(x)
          - crossprod(ginv(matrix(x[1:((D + 1) * (D + 1))], D + 1, D + 1)), x[-(1:((D + 1) * (D + 1)))]))
      
      w <- wh[1:D, ]
      b <- wh[D + 1, ]
      
      # Approximation of Likelihood
      
      detC <- apply(C, 2, function(x)
        det(matrix(x, D, D)))
      
      lxi <- (
        0.5 * log(detC) +
          0.5 * apply(rbind(C, mu), 2, function(x)
            t((x[-(1:(D * D))])) %*% ginv(matrix(x[1:(D * D)], nrow = D)) %*%
              (x[-(1:(D * D))])) +
          rowSums(
            log(sigma_xi) - 0.5 * xi - lambda_xi * xi ^ 2 +
              sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`) + sweep(X - 0.5, MARGIN = 2, b, `*`)
          )
      )
      
      C <- array(C, c(D, D, N))
      mu <- t(mu)
      
      ll <- sum(lxi)
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    } # end if D>2
    
  } # End of while statement
  
  # Gauss Hermite correction to the log-likelihood
  
  ghLL <- ghLLlta(X, w, b, D, pdGH)
  expected <- ghLL$expected
  LL <- ghLL$llGH
  
  LLva <- ll
  
  BICva <- -2 * LLva + (M * (D + 1) - D * (D - 1) / 2) * log(N)
  BIC <- -2 * LL + (M * (D + 1) - D * (D - 1) / 2) * log(N)
  
  b <- t(as.matrix(b))
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- ""
  rownames(w) <- NULL
  colnames(w) <- colnames(w, do.NULL = FALSE, prefix = "Item ")
  rownames(w) <- rownames(w, do.NULL = FALSE, prefix = "Dim ")
  names(LLva) <- "Log-Likelihood (variational approximation):"
  names(BICva) <- "BIC (variational approximation):"
  names(LL) <- "Log-Likelihood (G-H Quadrature correction):"
  names(BIC) <- "BIC (G-H Quadrature correction):"
  
  list(
    b = b,
    w = w,
    LLva = LLva,
    BICva = BICva,
    LL = LL,
    BIC = BIC,
    mu = mu,
    C = C,
    expected = expected
  )
}


f_lta_nstarts <- function(X, D, nstarts, tol, maxiter, pdGH)
{
  # Estimate the LTA model for different starting values
  # to avoid local maxima solutions
  
  out <- f_lta(X, D, tol, maxiter, pdGH)
  if(nstarts > 1){ 
    for(i in 2:nstarts){
      out1 <- f_lta(X, D, tol, maxiter, pdGH)
      if(out1$LL > out$LL) out <- out1
    }
  }
  out
}



lta <- function(X, D, nstarts = 3, tol = 0.1^2, maxiter = 250, pdGH = 21) {
  
  # LTA model
  
  if(any(D == 0)) stop("D must be > 0")
  
  if(length(D) == 1){ 
    out <- f_lta_nstarts(X, D, nstarts, tol, maxiter, pdGH)
  }else{
    out<-vector("list", length(D) + 1)
    names(out) <- c(paste('Dim y', D, sep = '='), 'BIC')
    i<-0
    for(diy in D){
      i <- i + 1
      out[[i]] <- f_lta_nstarts(X, diy, nstarts, tol, maxiter, pdGH)
    }
    out[[length(D) + 1]]<-tableBIC(out)
  }
  out
}


mlta <- function(X, DM, G, D, wfix = FALSE, nstarts = 3, tol = 0.1 ^ 2, maxiter = 250, pdGH = 21, 
                 beta0 = NULL){
  
  # MLTA model
  
  if (any(G < 1)) {
    print("Specify G > 0!")
    return("Specify G > 0!")
  }
  
  if (any(D < 0)) {
    print("Specify D >= 0!")
    return("Specify D >= 0!")
  }
  
  
  if (length(D) == 1 && length(G) == 1) {
    out <- f_mlta_methods(X, DM, G, D, nstarts, tol, maxiter, pdGH, wfix, beta0 = beta0)
  } else{
    out <- vector("list", length(D) * length(G) + 3)
    names(out) <- c(t(outer(
      paste('G=', G, sep = ''),
      paste('dim y=', D, sep = ''),
      paste,
      sep = ','
    )),
    'BIC', 'LL', 'LLva')
    
    bictab <- matrix(0, length(D), length(G))
    lltab <- matrix(0, length(D), length(G))
    
    rownames(bictab) <- paste('dim y=', D, sep = '')
    colnames(bictab) <- paste('G=', G, sep = '')
    
    rownames(lltab) <- paste('dim y=', D, sep = '')
    colnames(lltab) <- paste('G=', G, sep = '')
    
    llvatab <- lltab
    
    i <- 0
    for (g in G){
      for (diy in D) {
        i <- i + 1
        
        out[[i]] <- f_mlta_methods(X, DM, g, diy, nstarts, tol, maxiter, pdGH, wfix, beta0=beta0)
        print("ciao")
        bictab[i] <- out[[i]]$BIC
        lltab[i] <- out[[i]]$LL
        
        if (diy == 0) {
          llvatab[i] <- out[[i]]$LL
        } else{
          llvatab[i] <- out[[i]]$LLva
        }
        
        out[[length(G) * length(D) + 1]] <- ResTable(bictab, restype = 'BIC')
        out[[length(G) * length(D) + 2]] <- ResTable(lltab, restype = 'll')
        out[[length(G) * length(D) + 3]] <- ResTable(llvatab, restype = 'llva')
      }
    }
    class(out) <- "mmlta"
  }
  
  out
}



f_mlta <- function(X, DM, G, D, tol, maxiter, pdGH, beta0)
{
  
  # Unconstrained MLTA model
  # X=incidence matrix, DM=design matrix, G=number of components, 
  # D=latent trait dimension,
  # tol=tolerance level, maxiter=maximum number of iterations,
  # pdGH=number of quadrature points, beta0=initial beta if specified
  
  N <- nrow(X)
  M <- ncol(X) 
  J <- ncol(DM)
  
  # Initialize EM Algorithm
  
  if(!is.null(beta0)) beta = beta0 else{
    beta <- rep(0,J*(G-1))
    beta <- matrix(beta,ncol=G-1)
  }
  
  exb <- exp(DM %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  # priors
  
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N)   z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,]))
  
  
  p <- (t(z) %*% X) / colSums(z)
  
  
  # Initialize Variational Approximation
  
  xi <- array(20, c(N, M, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  w <- array(rnorm(M * D * G), c(D, M, G))
  b <- matrix(rnorm(M, G), G, M)
  
  C <- array(0, c(D, D, N, G))
  mu <- array(0, c(N, D, G))
  
  lxi <- matrix(0, G, N)
  
  # Iterative process
  
  v <- matrix(0, N, G)
  W <- list()
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  tol <- 0.1 ^ 6
  print(c(beta))
  
  xin <- xi
  
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    beta.old=beta
    ll_old <- ll
    
    if (D == 1) {
      for (g in 1:G)
      {
        w[, , g] <- as.matrix(t(w[, , g]))
        
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , , g] <-
          1 / (1 - 2 * rowSums(sweep(lambda_xi[, , g], MARGIN = 2, w[, , g] ^ 2, `*`)))
        mu[, , g] <-
          C[, , , g] * rowSums(sweep(
            X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`),
            MARGIN = 2,
            w[, , g],
            `*`
          ))
        
        mu[, , g] <- matrix(mu[, , g], N, D)
        YY <- matrix(C[, , , g] + mu[, , g] ^ 2, ncol = 1)
        
        xi[, , g] <-
          YY %*% (w[, , g] ^ 2) + mu[, , g] %*% t(2 * b[g, ] * w[, , g]) +
          matrix(b[g, ] ^ 2,
                 nrow = N,
                 ncol = M,
                 byrow = TRUE)
        
        xi[, , g] <- sqrt(xi[, , g])
        sigma_xi[, , g] <- 1 / (1 + exp(-xi[, , g]))
        lambda_xi[, , g] <- (0.5 - sigma_xi[, , g]) / (2 * xi[, , g])
        
        # STEP 3: Optimising the Model Parameters (w and b)
        
        YYh <- z[, g] * 2 * cbind(c(YY), mu[, , g], mu[, , g], 1)
        
        den <- t(YYh) %*% lambda_xi[, , g]
        num <- rbind(c(mu[, , g]), 1) %*% (z[, g] * (X - 0.5))
        wh <-
          apply(rbind(den, num), 2, function(x)
            - crossprod(ginv(matrix(x[1:4], 2, 2)), x[5:6]))
        
        w[, , g] <- wh[1:D, ]
        w[, , g] <- as.matrix(t(w[, , g]))
        b[g, ] <- wh[D + 1, ]
        
        # Approximation of log(p(x¦z))
        
        lxi[g, ] <-
          0.5 * log(C[, , , g]) + c(mu[, , g]) ^ 2 / (2 * C[, , , g]) +
          rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g] ^
              2 +
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) +
              sweep(X - 0.5, MARGIN = 2, b[g, ], `*`)
          )
        
      } # end for (g in 1:G)
      
      # M-step 
      
      lk = sum(z*log(eta))
      it = 0; lko = lk
      XXdis = array(0,c(G,(G-1)*ncol(DM),N))
      for(i in 1:N){
        XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(DM[i,]))
      }
      while((lk-lko>10^-6 & it<100) | it==0){
        it = it+1; lko = lk 
        sc = 0; Fi = 0
        for(i in 1:N){
          pdis = eta[i,]
          sc = sc+t(XXdis[,,i])%*%(z[i,]-pdis)
          Fi = Fi+t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
        }
        
        dbe = as.vector(ginv(Fi)%*%sc)
        mdbe = max(abs(dbe))
        if(mdbe>0.5) dbe = dbe/mdbe*0.5
        be0 = c(beta)
        flag = TRUE
        while(flag){
          beta = be0+dbe
          Eta = matrix(0,N,G)
          for(i in 1:N){
            
            if(ncol(DM)==1) Eta[i,] = XXdis[,,i]*beta
            else Eta[i,] = XXdis[,,i]%*%beta
          }	
          if(max(abs(Eta))>100){
            dbe = dbe/2
            flag = TRUE	
          }else{
            flag = FALSE
          }	        	
        }
        if(iter/10 == floor(iter/10))       print(beta)
        
        beta = matrix(beta, J, G-1)    
        exb <- exp(DM %*% beta) # updfe priors
        eta <- cbind(1,exb)/(rowSums(exb)+1)
        
        lk = sum(z*log(eta))
      }
      
      
      # E-step
      v <- eta * t(exp(lxi))
      if(any(is.nan(v))) browser()
      vsum <- apply(v, 1, sum)
      z <- v / vsum
      ll <- sum(rep(1, N) * log(vsum))
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    } # end while(diff>tol)
    
    if (D == 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C2 <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w[, , g]), t(w[, , g]))))
        C[, , , g] <- array(C2, c(D, D, N))
        
        mu[, , g] <-
          t(apply(rbind(C2, w[, , g] %*% t(
            X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)
          )), 2, function(x)
            matrix(x[1:4], nrow = D) %*% (x[-(1:4)])))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY <- C2 + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY, 2, function(x)
              rowSums((
                t(w[, , g]) %*% matrix(x, ncol = D)
              ) * t(w[, , g]))) + (2 * b[g, ] * t(w[, , g])) %*% t(mu[, , g]) + matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = N,
                byrow = FALSE
              )
          )
        
        xi[, , g] <- sqrt(xi[, , g])
        sigma_xi[, , g] <- 1 / (1 + exp(-xi[, , g]))
        lambda_xi[, , g] <- (0.5 - sigma_xi[, , g]) / (2 * xi[, , g])
        
        # # STEP 3: Optimising the Model Parameters (w and b)
        
        YYh <-
          z[, g] * 2 * cbind(YY[1, ], YY[2, ], mu[, 1, g], YY[3, ], YY[4, ], mu[, 2, g], mu[, 1, g], mu[, 2, g], 1)
        
        den <- t(YYh) %*% lambda_xi[, , g]
        num <- rbind(t(mu[, , g]), 1) %*% {
          z[, g] * {
            X - 0.5
          }
        }
        wh <-
          apply(rbind(den, num), 2, function(x)
            - crossprod(ginv(matrix(x[1:9], 3, 3)), x[10:12]))
        
        w[, , g] <- wh[1:D, ]
        w[, , g] <- as.matrix(w[, , g])
        b[g, ] <- wh[D + 1, ]
        
        # Approximation of log(p(y|z))
        
        detC <- C2[1, ] * C2[4, ] - C2[3, ] * C2[2, ]
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C2[4, ] / detC, -C2[2, ] / detC, -C2[3, ] /
                                                detC, C2[1, ] / detC, t(mu[, , g])), 2, function(x)
                                                  t((x[-{
                                                    1:4
                                                  }])) %*% matrix(x[1:4], nrow = D) %*% (x[-{
                                                    1:4
                                                  }])) + rowSums(
                                                    log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * {
                                                      xi[, , g] ^ 2
                                                    } + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + sweep(X - 0.5, MARGIN =
                                                                                                                       2, b[g, ], `*`)
                                                  )
        
      } # end for (g in 1:G)
      
      # M-step 
      
      lk = sum(z*log(eta))
      it = 0; lko = lk
      XXdis = array(0,c(G,(G-1)*ncol(DM),N))
      for(i in 1:N){
        XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(DM[i,]))
      }
      while((lk-lko>10^-6 & it<100) | it==0){
        it = it+1; lko = lk 
        sc = 0; Fi = 0
        for(i in 1:N){
          pdis = eta[i,]
          sc = sc+t(XXdis[,,i])%*%(z[i,]-pdis)
          Fi = Fi+t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
        }
        
        dbe = as.vector(ginv(Fi)%*%sc)
        mdbe = max(abs(dbe))
        if(mdbe>0.5) dbe = dbe/mdbe*0.5
        be0 = c(beta)
        flag = TRUE
        while(flag){
          beta = be0+dbe
          Eta = matrix(0,N,G)
          for(i in 1:N){
            
            if(ncol(DM)==1) Eta[i,] = XXdis[,,i]*beta
            else Eta[i,] = XXdis[,,i]%*%beta
          }	
          if(max(abs(Eta))>100){
            dbe = dbe/2
            flag = TRUE	
          }else{
            flag = FALSE
          }	        	
        }
        if(iter/10 == floor(iter/10))       print(beta)
        
        beta = matrix(beta, J, G-1)    
        exb <- exp(DM %*% beta) # updfe priors
        eta <- cbind(1,exb)/(rowSums(exb)+1)
        
        lk = sum(z*log(eta))
      }
      
      
      # E-step
      v <- eta * t(exp(lxi))
      if(any(is.nan(v))) browser()
      vsum <- apply(v, 1, sum)
      z <- v / vsum
      ll <- sum(rep(1, N) * log(vsum))
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    } # end while(diff>tol)
    
    if (D > 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C2 <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w[, , g]), t(w[, , g]))))
        C[, , , g] <- array(C2, c(D, D, N))
        
        mu[, , g] <-
          t(apply(rbind(C2, w[, , g] %*% t(
            X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)
          )), 2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% (
              x[-(1:(D * D))])))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY <- C2 + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY, 2, function(x)
              rowSums((t(w[, , g]) %*% matrix(x, ncol = D)) * t(w[, , g]))) + (2 * b[g, ] * t(w[, , g])) %*% t(mu[, , g]) +
              matrix(b[g, ] ^ 2,
                     nrow = M,
                     ncol = N,
                     byrow = FALSE)
          )
        
        xi[, , g] <- sqrt(xi[, , g])
        sigma_xi[, , g] <- 1 / (1 + exp(-xi[, , g]))
        lambda_xi[, , g] <- (0.5 - sigma_xi[, , g]) / (2 * xi[, , g])
        
        # # STEP 3: Optimising the Model Parameters (w and b)
        
        YYh <-
          sweep(buildYYh(YY, t(mu[, , g]), D, N), MARGIN = 2, 2 * z[, g], `*`)
        den <- (YYh) %*% lambda_xi[, , g]
        num <- rbind(t(mu[, , g]), 1) %*% (z[, g] * (X - 0.5))
        wh <-
          apply(rbind(den, num), 2, function(x)
            - crossprod(ginv(matrix(x[1:((D + 1) * (D + 1))], D + 1, D + 1)), x[-(1:((D +
                                                                                         1) * (D + 1)))]))
        
        w[, , g] <- wh[1:D, ]
        w[, , g] <- as.matrix(w[, , g])
        b[g, ] <- wh[D + 1, ]
        
        # Approximation of log(p(y|z))
        
        detC <- apply(C2, 2, function(x)
          det(matrix(x, D, D)))
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C2, t(mu[, , g])), 2, function(x)
            t((x[-(1:(D * D))])) %*% ginv(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(1:(D * D))])) + 
          rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * (xi[, , g] ^
                                                                           2) + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + sweep(X - 0.5, MARGIN =
                                                                                                                                               2, b[g, ], `*`)
          )
      } # end for (g in 1:G)
      
      # M-step 
      
      lk = sum(z*log(eta))
      it = 0; lko = lk
      XXdis = array(0,c(G,(G-1)*ncol(DM),N))
      for(i in 1:N){
        XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(DM[i,]))
      }
      while((lk-lko>10^-6 & it<100) | it==0){
        it = it+1; lko = lk 
        sc = 0; Fi = 0
        for(i in 1:N){
          pdis = eta[i,]
          sc = sc+t(XXdis[,,i])%*%(z[i,]-pdis)
          Fi = Fi+t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
        }
        
        dbe = as.vector(ginv(Fi)%*%sc)
        mdbe = max(abs(dbe))
        if(mdbe>0.5) dbe = dbe/mdbe*0.5
        be0 = c(beta)
        flag = TRUE
        while(flag){
          beta = be0+dbe
          Eta = matrix(0,N,G)
          for(i in 1:N){
            
            if(ncol(DM)==1) Eta[i,] = XXdis[,,i]*beta
            else Eta[i,] = XXdis[,,i]%*%beta
          }	
          if(max(abs(Eta))>100){
            dbe = dbe/2
            flag = TRUE	
          }else{
            flag = FALSE
          }	        	
        }
        if(iter/10 == floor(iter/10))       print(beta)
        
        beta = matrix(beta, J, G-1)    
        exb <- exp(DM %*% beta) # update priors
        eta <- cbind(1,exb)/(rowSums(exb)+1)
        
        lk = sum(z*log(eta))
      }
      
      
      # E-step
      v <- eta * t(exp(lxi))
      if(any(is.nan(v))) browser()
      vsum <- apply(v, 1, sum)
      z <- v / vsum
      ll <- sum(rep(1, N) * log(vsum))
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    } # end while(diff>tol)
    
  } # end while(diff>tol)
  
  # Correction to the log-likelihood
  
  # Gauss-Hermite Quadrature
  
  npoints <- round(pdGH ^ (1 / D))
  ny <- npoints ^ D
  GaHer <- glmmML::ghq(npoints, FALSE)
  Ygh <- expand.grid(rep(list(GaHer$zeros), D))
  Ygh <- as.matrix(Ygh)
  Wgh <- apply(as.matrix(expand.grid(rep(
    list(GaHer$weights), D
  ))), 1, prod) *
    apply(exp(Ygh ^ 2), 1, prod)
  
  Hy <- apply(Ygh, 1, mvtnorm::dmvnorm)
  Beta <- Hy * Wgh / sum(Hy * Wgh)
  
  fxy <- array(0, c(N, ny, G))
  
  for (g in 1:G)
  {
    if (D == 1) {
      Agh <- t(tcrossprod(w[, , g], Ygh) + b[g, ])
    } else{
      Agh <- t(tcrossprod(t(w[, , g]), Ygh) + b[g, ])
    }
    pgh <- 1 / (1 + exp(-Agh))
    
    fxy[, , g] <-
      exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh)))
    fxy[is.nan(fxy[, , g])] <- 0
  }
  
  LLva <- ll
  BICva <-
    -2 * LLva + (G * (M * (D + 1) - D * (D - 1) / 2) + J*(G - 1)) * log(N)
  
  llGH1 <- apply(rep(Beta, each = N) * fxy, c(1, 3), sum)
  LL <- sum(log(rowSums(eta * llGH1)))
  
  
  BIC <-
    -2 * LL + (G * (M * (D + 1) - D * (D - 1) / 2) + J*(G - 1)) * log(N)
  
  expected <- colSums(eta * (llGH1)) * N
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- rownames(b, do.NULL = FALSE, prefix = "Group ")
  
  rownames(w) <- rownames(w, do.NULL = FALSE, prefix = "Dim ")
  colnames(w) <- colnames(w, do.NULL = FALSE, prefix = "Item ")
  dimnames(w)[[3]] <- paste("Group ", seq(1, G) , sep = '')
  
  colnames(beta) <- 2:G
  rownames(beta) <- c(colnames(DM))
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(BICva) <- c("BIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(BIC) <- c("BIC (G-H Quadrature correction):")
  
  out <- list(
    b = b,
    w = w,
    eta = eta,
    mu = mu,
    C = C,
    z = z,
    LL = LL,
    BIC = BIC,
    LLva = LLva,
    BICva = BICva,
    expected = expected,
    beta = beta, 
    rrr=exp(beta)
  )
  
  out
}


tableBIC <- function(out){
  
  # BIC table for model selection
  
  lout <- length(out) - 1
  bicG <- numeric(lout)
  names(bicG) <- names(out[- (lout + 1)])
  
  for(i in 1:lout) bicG[i] <- out[[i]]$BIC
  
  resBIC <- vector('list', 2)
  names(resBIC) <- c('Model Selection', 'Table of BIC Results')
  resBIC[[1]] <- names(bicG)[bicG == min(bicG)]
  names(resBIC[[1]]) <- 'Model with lower BIC:'
  resBIC[[2]] <- bicG
  resBIC
}



print.mlta <- function(x){
  stopifnot(inherits(x, 'mlta'))
  cat("b:\n")
  print(x$b) 
  cat("\nw:\n")  
  print(x$w)
  cat("\neta:\n")
  print(x$eta)
  cat("\n")
  print(x$LL)
  print(x$BIC) 
}



print.mmlta <- function(x){
  cat("Log-Likelihood:\n")
  print(x$LL$`Table of LL (G-H Quadrature correction)`)
  cat("BIC:\n")
  print(x$BIC$`Table of BIC Results`)
  print(x$BIC$`Model Selection`)
}


## Likelihood Ratio Test ####
# load("finland_mod.RData")
# load("italy_mod.RData")
# load("bulgaria_mod.RData")

modLL=matrix(NA, 3, 3)
for (l in 1:9) {
  modLL[l]=mod[[l]]$LL
}

lcaLL=matrix(NA, 1, 3)
for (l in 1:3) {
  lcaLL[l]=LCA[[l]]$LL
}

ltaLL=matrix(NA, 1, 3)
for (l in 1:3) {
  ltaLL[l]=LTA[[l]]$LL
}

mod.wfixLL=matrix(NA, 3, 3)
for (l in 1:9) {
  mod.wfixLL[l]=mod.wfix[[l]]$LL
}

Nfin=1446
Nita=1465
Nbulg=2082

bic.fin=8744.80
bic.ita=12858.12
bic.bulg=17972.87

LLmod=cbind(c((((1+1+7) * log(Nfin)) - bic.fin)/2,ltaLL),rbind(lcaLL,modLL)) # complex
LLmod.wfix=cbind(c((((1+1+7) * log(Nfin)) - bic.fin)/2,ltaLL),rbind(lcaLL,mod.wfixLL)) # nested

LLmod=cbind(c((((1+1+7) * log(Nita)) - bic.ita)/2,ltaLL),rbind(lcaLL,modLL)) # complex
LLmod.wfix=cbind(c((((1+1+7) * log(Nita)) - bic.ita)/2,ltaLL),rbind(lcaLL,mod.wfixLL)) # nested

LLmod=cbind(c((((1+1+7) * log(Nbulg)) - bic.bulg)/2,ltaLL),rbind(lcaLL,modLL)) # complex
LLmod.wfix=cbind(c((((1+1+7) * log(Nbulg)) - bic.bulg)/2,ltaLL),rbind(lcaLL,mod.wfixLL)) # nested


teststat <- -2 * (LLmod.wfix-LLmod)
p.val=matrix(NA, 4, 4)
for (g in 1:ncol(LLmod)) {
  for (d in 1:nrow(LLmod)) {
    df = d * 7 * (g-1)
    p.val[d,g] <- pchisq(teststat[d,g], df = df, lower.tail = FALSE)
  }
}
round(p.val, 2)


## ESS Application ####

## Load the data set for the country of interest:

# load("finland.RData") # for Finland
# load("italy.RData") # for Italy
# load("bulgaria.RData") # for Bulgaria 


## Response variables plots ####
# eps size: 1300 650

data.int <- data.frame(
  Indice = c("Never", "Occasionally", "A few times \n a week", "Most days", "Every day"),
  Valore = as.vector(table(ESS10$netusoft)/nrow(ESS10))
)
ordine_desiderato <- c("Never", "Occasionally", "A few times \n a week", "Most days", "Every day")
data.int$Indice <- factor(data.int$Indice, levels = ordine_desiderato)
ggplot(data.int, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Internet use") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.pref <- data.frame(
  Indice = c("Not at all", "Not very", "Somewhat", "Very", "Completely"),
  Valore = as.vector(table(ESS10$fampref)/nrow(ESS10))
)
ordine_desiderato <- c("Not at all", "Not very", "Somewhat", "Very", "Completely")
data.pref$Indice <- factor(data.pref$Indice, levels = ordine_desiderato)
ggplot(data.pref, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Preference settings") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.adv <- data.frame(
  Indice = c("Not at all", "Not very", "Somewhat", "Very", "Completely"),
  Valore = as.vector(table(ESS10$famadvs)/nrow(ESS10))
)
ordine_desiderato <- c("Not at all", "Not very", "Somewhat", "Very", "Completely")
data.adv$Indice <- factor(data.adv$Indice, levels = ordine_desiderato)
ggplot(data.adv, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Advanced search") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.pdf <- data.frame(
  Indice = c("Not at all", "Not very", "Somewhat", "Very", "Completely"),
  Valore = as.vector(table(ESS10$fampdf)/nrow(ESS10))
)
ordine_desiderato <- c("Not at all", "Not very", "Somewhat", "Very", "Completely")
data.pdf$Indice <- factor(data.pdf$Indice, levels = ordine_desiderato)
ggplot(data.pdf, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("PDFs") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.videochild <- data.frame(
  Indice = c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
             "Never", "Not \n applicable"),
  Valore = as.vector(table(ESS10$scrno12)/nrow(ESS10))
)
ordine_desiderato <- c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
                       "Never", "Not \n applicable")
data.videochild$Indice <- factor(data.videochild$Indice, levels = ordine_desiderato)
ggplot(data.videochild, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Video call with child 12+") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))



data.videoparent <- data.frame(
  Indice = c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
             "Never", "Not \n applicable"),
  Valore = as.vector(table(ESS10$scrnpnt)/nrow(ESS10))
)
ordine_desiderato <- c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
                       "Never", "Not \n applicable")
data.videoparent$Indice <- factor(data.videoparent$Indice, levels = ordine_desiderato)
ggplot(data.videoparent, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Video call with parents") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.videomanager <- data.frame(
  Indice = c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
             "Never", "Not \n applicable"),
  Valore = as.vector(table(ESS10$manscrn)/nrow(ESS10))
)
ordine_desiderato <- c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
                       "Never", "Not \n applicable")
data.videomanager$Indice <- factor(data.videomanager$Indice, levels = ordine_desiderato)
ggplot(data.videomanager, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Video call with manager") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))



data.videocoll <- data.frame(
  Indice = c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
             "Never", "Not \n applicable"),
  Valore = as.vector(table(ESS10$colscrn)/nrow(ESS10))
)
ordine_desiderato <- c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
                       "Never", "Not \n applicable")
data.videocoll$Indice <- factor(data.videocoll$Indice, levels = ordine_desiderato)
ggplot(data.videocoll, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Video call with colleagues") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.messchild <- data.frame(
  Indice = c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
             "Never", "Not \n applicable"),
  Valore = as.vector(table(ESS10$como12)/nrow(ESS10))
)
ordine_desiderato <- c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
                       "Never", "Not \n applicable")
data.messchild$Indice <- factor(data.messchild$Indice, levels = ordine_desiderato)
ggplot(data.messchild, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Messages with child 12+") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))



data.messpar <- data.frame(
  Indice = c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
             "Never", "Not \n applicable"),
  Valore = as.vector(table(ESS10$compnt)/nrow(ESS10))
)
ordine_desiderato <- c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
                       "Never", "Not \n applicable")
data.messpar$Indice <- factor(data.messpar$Indice, levels = ordine_desiderato)
ggplot(data.messpar, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Messages with parents") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.messmanager <- data.frame(
  Indice = c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
             "Never", "Not \n applicable"),
  Valore = as.vector(table(ESS10$mancom)/nrow(ESS10))
)
ordine_desiderato <- c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
                       "Never", "Not \n applicable")
data.messmanager$Indice <- factor(data.messmanager$Indice, levels = ordine_desiderato)
ggplot(data.messmanager, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Messages with manager") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.messcoll <- data.frame(
  Indice = c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
             "Never", "Not \n applicable"),
  Valore = as.vector(table(ESS10$colcom)/nrow(ESS10))
)
ordine_desiderato <- c("Often in \n a day", "Once \n a day", "Often in \n a week", "Often in \n a month", "Once \n a month", "Less \n often",
                       "Never", "Not \n applicable")
data.messcoll$Indice <- factor(data.messcoll$Indice, levels = ordine_desiderato)
ggplot(data.messcoll, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Messages with colleagues") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data.post <- data.frame(
  Indice = c("Yes", "No"),
  Valore = as.vector(table(ESS10$pstplonl)/nrow(ESS10))
)
ordine_desiderato <- c("Yes", "No")
data.post$Indice <- factor(data.post$Indice, levels = ordine_desiderato)
ggplot(data.post, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Online politcal posts") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


### Likert plots ####

# load("finland.RData")
ess.finland=ESS10
# load("italy.RData")
ess.italy=ESS10
# load("bulgaria.RData")
ess.bulgaria=ESS10

#load finland.RData
dati.fin=dati
#load italy.RData
dati.ita=dati
#load bulgaria.RData
dati.bulg=dati

levels_order <- c("Bulgaria", "Italy", "Finland")

# eps 965 351

# Internet
ess.finland$netusoft=as.factor(ess.finland$netusoft)
ess.italy$netusoft=as.factor(ess.italy$netusoft)
ess.bulgaria$netusoft=as.factor(ess.bulgaria$netusoft)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$netusoft, ess.finland$country), 
               cbind(ess.italy$netusoft, ess.italy$country),
               cbind(ess.bulgaria$netusoft, ess.bulgaria$country)))
colnames(new.data)=c("Internet", "Country")
str(new.data)
new.data$Internet=as.factor(new.data$Internet)
levels(new.data$Internet)=c("Never","Occasionally","A few times a week","Most days","Every day")
"Internet use"=new.data$Internet
new.data$Country <- factor(new.data$Country, levels = levels_order)
internet <- likert(data.frame(`Internet use`), grouping=new.data$Country)
plot(internet, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


# Preference settings
ess.finland$fampref=as.factor(ess.finland$fampref)
ess.italy$fampref=as.factor(ess.italy$fampref)
ess.bulgaria$fampref=as.factor(ess.bulgaria$fampref)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$fampref, ess.finland$country), 
                          cbind(ess.italy$fampref, ess.italy$country),
                          cbind(ess.bulgaria$fampref, ess.bulgaria$country)))
colnames(new.data)=c("Preference settings", "Country")
str(new.data)
new.data$`Preference settings`=as.factor(new.data$`Preference settings`)
levels(new.data$`Preference settings`)=c("Not familiar", "Not very familiar", "Somewhat familiar", "Very familiar", "Completely familiar")
"Preference settings"=new.data$`Preference settings`
new.data$Country <- factor(new.data$Country, levels = levels_order)
pref <- likert(data.frame(`Preference settings`), grouping=new.data$Country)
plot(pref, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Advanced search
ess.finland$famadvs=as.factor(ess.finland$famadvs)
ess.italy$famadvs=as.factor(ess.italy$famadvs)
ess.bulgaria$famadvs=as.factor(ess.bulgaria$famadvs)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$famadvs, ess.finland$country), 
                          cbind(ess.italy$famadvs, ess.italy$country),
                          cbind(ess.bulgaria$famadvs, ess.bulgaria$country)))
colnames(new.data)=c("Advanced search", "Country")
str(new.data)
new.data$`Advanced search`=as.factor(new.data$`Advanced search`)
levels(new.data$`Advanced search`)=c("Not familiar", "Not very familiar", "Somewhat familiar", "Very familiar", "Completely familiar")
"Advanced search"=new.data$`Advanced search`
new.data$Country <- factor(new.data$Country, levels = levels_order)
adv <- likert(data.frame(`Advanced search`), grouping=new.data$Country)
plot(adv, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# PDFs
ess.finland$fampdf=as.factor(ess.finland$fampdf)
ess.italy$fampdf=as.factor(ess.italy$fampdf)
ess.bulgaria$fampdf=as.factor(ess.bulgaria$fampdf)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$fampdf, ess.finland$country), 
                          cbind(ess.italy$fampdf, ess.italy$country),
                          cbind(ess.bulgaria$fampdf, ess.bulgaria$country)))
colnames(new.data)=c("PDFs", "Country")
str(new.data)
new.data$`PDFs`=as.factor(new.data$`PDFs`)
levels(new.data$`PDFs`)=c("Not familiar", "Not very familiar", "Somewhat familiar", "Very familiar", "Completely familiar")
"PDFs"=new.data$`PDFs`
new.data$Country <- factor(new.data$Country, levels = levels_order)
pdf <- likert(data.frame(`PDFs`), grouping=new.data$Country)
plot(pdf, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


# Video with child 12+
ess.finland$scrno12=as.factor(ess.finland$scrno12)
ess.italy$scrno12=as.factor(ess.italy$scrno12)
ess.bulgaria$scrno12=as.factor(ess.bulgaria$scrno12)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$scrno12, ess.finland$country), 
                          cbind(ess.italy$scrno12, ess.italy$country),
                          cbind(ess.bulgaria$scrno12, ess.bulgaria$country)))
colnames(new.data)=c("Video call with child over 12", "Country")
str(new.data)
new.data$`Video call with child over 12`=as.factor(new.data$`Video call with child over 12`)
levels(new.data$`Video call with child over 12`)=c("Often in a day", "Once a day", "Often in a week", "Often in a month", "Once a month", "Less often",
                                               "Never", "Not applicable")
new.data$`Video call with child over 12` <- factor(new.data$`Video call with child over 12`, 
                                                   levels = rev(levels(new.data$`Video call with child over 12`)), 
                                                   ordered = TRUE)
"Video call with child over 12"=new.data$`Video call with child over 12`
new.data$Country <- factor(new.data$Country, levels = levels_order)
videochild <- likert(data.frame(`Video call with child over 12`), grouping=new.data$Country)
plot(videochild, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Video with parent
ess.finland$scrnpnt=as.factor(ess.finland$scrnpnt)
ess.italy$scrnpnt=as.factor(ess.italy$scrnpnt)
ess.bulgaria$scrnpnt=as.factor(ess.bulgaria$scrnpnt)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$scrnpnt, ess.finland$country), 
                          cbind(ess.italy$scrnpnt, ess.italy$country),
                          cbind(ess.bulgaria$scrnpnt, ess.bulgaria$country)))
colnames(new.data)=c("Video call with parents", "Country")
str(new.data)
new.data$`Video call with parents`=as.factor(new.data$`Video call with parents`)
levels(new.data$`Video call with parents`)=c("Often in a day", "Once a day", "Often in a week", "Often in a month", "Once a month", "Less often",
                                               "Never", "Not applicable")
new.data$`Video call with parents` <- factor(new.data$`Video call with parents`, 
                                                   levels = rev(levels(new.data$`Video call with parents`)), 
                                                   ordered = TRUE)
"Video call with parents"=new.data$`Video call with parents`
new.data$Country <- factor(new.data$Country, levels = levels_order)
videoparent <- likert(data.frame(`Video call with parents`), grouping=new.data$Country)
plot(videoparent, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Video with manager
ess.finland$manscrn=as.factor(ess.finland$manscrn)
ess.italy$manscrn=as.factor(ess.italy$manscrn)
ess.bulgaria$manscrn=as.factor(ess.bulgaria$manscrn)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$manscrn, ess.finland$country), 
                          cbind(ess.italy$manscrn, ess.italy$country),
                          cbind(ess.bulgaria$manscrn, ess.bulgaria$country)))
colnames(new.data)=c("Video call with manager", "Country")
str(new.data)
new.data$`Video call with manager`=as.factor(new.data$`Video call with manager`)
levels(new.data$`Video call with manager`)=c("Often in a day", "Once a day", "Often in a week", "Often in a month", "Once a month", "Less often",
                                             "Never", "Not applicable")
new.data$`Video call with manager` <- factor(new.data$`Video call with manager`, 
                                             levels = rev(levels(new.data$`Video call with manager`)), 
                                             ordered = TRUE)
"Video call with manager"=new.data$`Video call with manager`
new.data$Country <- factor(new.data$Country, levels = levels_order)
videomanager <- likert(data.frame(`Video call with manager`), grouping=new.data$Country)
plot(videomanager, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


# Video with colleagues
ess.finland$colscrn=as.factor(ess.finland$colscrn)
ess.italy$colscrn=as.factor(ess.italy$colscrn)
ess.bulgaria$colscrn=as.factor(ess.bulgaria$colscrn)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$colscrn, ess.finland$country), 
                          cbind(ess.italy$colscrn, ess.italy$country),
                          cbind(ess.bulgaria$colscrn, ess.bulgaria$country)))
colnames(new.data)=c("Video call with colleagues", "Country")
str(new.data)
new.data$`Video call with colleagues`=as.factor(new.data$`Video call with colleagues`)
levels(new.data$`Video call with colleagues`)=c("Often in a day", "Once a day", "Often in a week", "Often in a month", "Once a month", "Less often",
                                             "Never", "Not applicable")
new.data$`Video call with colleagues` <- factor(new.data$`Video call with colleagues`, 
                                             levels = rev(levels(new.data$`Video call with colleagues`)), 
                                             ordered = TRUE)
"Video call with colleagues"=new.data$`Video call with colleagues`
new.data$Country <- factor(new.data$Country, levels = levels_order)
videocolleagues <- likert(data.frame(`Video call with colleagues`), grouping=new.data$Country)
plot(videocolleagues, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


# Messages with child 12+
ess.finland$como12=as.factor(ess.finland$como12)
ess.italy$como12=as.factor(ess.italy$como12)
ess.bulgaria$como12=as.factor(ess.bulgaria$como12)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$como12, ess.finland$country), 
                          cbind(ess.italy$como12, ess.italy$country),
                          cbind(ess.bulgaria$como12, ess.bulgaria$country)))
colnames(new.data)=c("Messages with child over 12", "Country")
str(new.data)
new.data$`Messages with child over 12`=as.factor(new.data$`Messages with child over 12`)
levels(new.data$`Messages with child over 12`)=c("Often in a day", "Once a day", "Often in a week", "Often in a month", "Once a month", "Less often",
                                                "Never", "Not applicable")
new.data$`Messages with child over 12` <- factor(new.data$`Messages with child over 12`, 
                                                levels = rev(levels(new.data$`Messages with child over 12`)), 
                                                ordered = TRUE)
"Messages with child over 12"=new.data$`Messages with child over 12`
new.data$Country <- factor(new.data$Country, levels = levels_order)
messchild <- likert(data.frame(`Messages with child over 12`), grouping=new.data$Country)
plot(messchild, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Messages with parents
ess.finland$compnt=as.factor(ess.finland$compnt)
ess.italy$compnt=as.factor(ess.italy$compnt)
ess.bulgaria$compnt=as.factor(ess.bulgaria$compnt)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$compnt, ess.finland$country), 
                          cbind(ess.italy$compnt, ess.italy$country),
                          cbind(ess.bulgaria$compnt, ess.bulgaria$country)))
colnames(new.data)=c("Messages with parents", "Country")
str(new.data)
new.data$`Messages with parents`=as.factor(new.data$`Messages with parents`)
levels(new.data$`Messages with parents`)=c("Often in a day", "Once a day", "Often in a week", "Often in a month", "Once a month", "Less often",
                                             "Never", "Not applicable")
new.data$`Messages with parents` <- factor(new.data$`Messages with parents`, 
                                                 levels = rev(levels(new.data$`Messages with parents`)), 
                                                 ordered = TRUE)
"Messages with parents"=new.data$`Messages with parents`
new.data$Country <- factor(new.data$Country, levels = levels_order)
messpar <- likert(data.frame(`Messages with parents`), grouping=new.data$Country)
plot(messpar, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


# Messages with manager
ess.finland$mancom=as.factor(ess.finland$mancom)
ess.italy$mancom=as.factor(ess.italy$mancom)
ess.bulgaria$mancom=as.factor(ess.bulgaria$mancom)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$mancom, ess.finland$country), 
                          cbind(ess.italy$mancom, ess.italy$country),
                          cbind(ess.bulgaria$mancom, ess.bulgaria$country)))
colnames(new.data)=c("Messages with manager", "Country")
str(new.data)
new.data$`Messages with manager`=as.factor(new.data$`Messages with manager`)
levels(new.data$`Messages with manager`)=c("Often in a day", "Once a day", "Often in a week", "Often in a month", "Once a month", "Less often",
                                           "Never", "Not applicable")
new.data$`Messages with manager` <- factor(new.data$`Messages with manager`, 
                                           levels = rev(levels(new.data$`Messages with manager`)), 
                                           ordered = TRUE)
"Messages with manager"=new.data$`Messages with manager`
new.data$Country <- factor(new.data$Country, levels = levels_order)
messmanager <- likert(data.frame(`Messages with manager`), grouping=new.data$Country)
plot(messmanager, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Messages with colleagues
ess.finland$colcom=as.factor(ess.finland$colcom)
ess.italy$colcom=as.factor(ess.italy$colcom)
ess.bulgaria$colcom=as.factor(ess.bulgaria$colcom)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$colcom, ess.finland$country), 
                          cbind(ess.italy$colcom, ess.italy$country),
                          cbind(ess.bulgaria$colcom, ess.bulgaria$country)))
colnames(new.data)=c("Messages with colleagues", "Country")
str(new.data)
new.data$`Messages with colleagues`=as.factor(new.data$`Messages with colleagues`)
levels(new.data$`Messages with colleagues`)=c("Often in a day", "Once a day", "Often in a week", "Often in a month", "Once a month", "Less often",
                                           "Never", "Not applicable")
new.data$`Messages with colleagues` <- factor(new.data$`Messages with colleagues`, 
                                           levels = rev(levels(new.data$`Messages with colleagues`)), 
                                           ordered = TRUE)
"Messages with colleagues"=new.data$`Messages with colleagues`
new.data$Country <- factor(new.data$Country, levels = levels_order)
messcoll <- likert(data.frame(`Messages with colleagues`), grouping=new.data$Country)
plot(messcoll, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Online posts
ess.finland$pstplonl=as.factor(ess.finland$pstplonl)
ess.italy$pstplonl=as.factor(ess.italy$pstplonl)
ess.bulgaria$pstplonl=as.factor(ess.bulgaria$pstplonl)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$pstplonl, ess.finland$country), 
                          cbind(ess.italy$pstplonl, ess.italy$country),
                          cbind(ess.bulgaria$pstplonl, ess.bulgaria$country)))
colnames(new.data)=c("Online political posts", "Country")
str(new.data)
new.data$`Online political posts`=as.factor(new.data$`Online political posts`)
levels(new.data$`Online political posts`)=c("Yes","No")
new.data$`Online political posts` <- factor(new.data$`Online political posts`, 
                                              levels = rev(levels(new.data$`Online political posts`)), 
                                              ordered = TRUE)
"Online political posts"=new.data$`Online political posts`
new.data$Country <- factor(new.data$Country, levels = levels_order)
post <- likert(data.frame(`Online political posts`), grouping=new.data$Country)
plot(post, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


## Concomitant variables plots ####

summary(ESS10$agea)
data <- data.frame(
  agea = ESS10$agea
)
data$agea_grouped <- cut(data$agea, breaks = seq(0, 100, by = 5), right = FALSE)
ggplot(data, aes(x = agea_grouped)) +
  geom_bar(fill = "gray70") +
  scale_x_discrete(labels = scales::label_wrap(8)) +
  labs(x = "", y = "Frequency") +
  ggtitle("Age") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5))

# All Countries together: eps 1400 400
# data <- dati %>%
#   filter(!is.na(agea) & !is.na(country)) %>%
#   mutate(agea_grouped = cut(agea, breaks = seq(0, 100, by = 5), right = FALSE))
# data_rel_freq <- data %>%
#   group_by(country, agea_grouped) %>%
#   summarise(count = n()) %>%
#   mutate(freq_relativa = count / sum(count)) %>%
#   ungroup()
# order_of_countries <- c("Finland", "Italy", "Bulgaria") # Sostituisci con i nomi dei paesi nel tuo dataset
# data_rel_freq$country <- factor(data_rel_freq$country, levels = order_of_countries)
# ggplot(data_rel_freq, aes(x = agea_grouped, y = freq_relativa)) +
#   geom_bar(stat = "identity", fill = "gray70") +
#   scale_x_discrete(labels = scales::label_wrap(8)) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Age Group", y = "Relative Frequency") +
#   ggtitle("Age Distribution by Country") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(size = 5),
#     axis.text.y = element_text(size = 28),
#     axis.title = element_text(size = 28),
#     plot.title = element_text(size = 28, hjust = 0.5),
#     strip.text = element_text(size = 28)
#   ) +
#   facet_wrap(~country, scales = "free_x")

data <- data.frame(
  Indice = c("Very good", "Good",	"Fair", "Bad", "Very bad"),
  Valore = as.vector(table(ESS10$health)/nrow(ESS10))
)
ordine_desiderato <- c("Very good", "Good",	"Fair", "Bad", "Very bad")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Health") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Indice = c("Yes a lot", "Yes to some extent", "No"),
  Valore = as.vector(table(ESS10$hlthhmp)/nrow(ESS10))
)
ordine_desiderato <- c("Yes a lot", "Yes to some extent", "No")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Hampered in daily activities") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

data <- data.frame(
  Indice = c("Yes","No"),
  Valore = as.vector(table(ESS10$brncntr)/nrow(ESS10))
)
ordine_desiderato <- c("Yes", "No")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Born in the country") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Indice = c("Male","Female"),
  Valore = as.vector(table(ESS10$gndr)/nrow(ESS10))
)
ordine_desiderato <- c("Male", "Female")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Gender") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Indice = c("< Lower \n secondary", "Lower \n secondary", 
             "Upper \n secondary",
             "Sub-degree", "BA \n level", ">= MA \n level"),
  Valore = as.vector(table(ESS10$eisced)/nrow(ESS10))
)
ordine_desiderato <- c("< Lower \n secondary", "Lower \n secondary", 
                       "Upper \n secondary",
                       "Sub-degree", "BA \n level", ">= MA \n level")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Education") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 22),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Indice = c("< Lower \n secondary", "Lower \n secondary", 
             "Upper \n secondary",
             "Sub-degree", "BA \n level", ">= MA \n level", "Not \n applicable"),
  Valore = as.vector(table(ESS10$eiscedp)/nrow(ESS10))
)
ordine_desiderato <- c("< Lower \n secondary", "Lower \n secondary", 
                       "Upper \n secondary",
                       "Sub-degree", "BA \n level", ">= MA \n level", "Not \n applicable")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Partner education") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 22),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Indice = c("1st", "2nd", "3rd", "4th", "5th",
             "6th", "7th", "8th", "9th", "10th"),
  Valore = as.vector(table(ESS10$hinctnta)/nrow(ESS10))
)
ordine_desiderato <- c("1st", "2nd", "3rd", "4th", "5th",
                       "6th", "7th", "8th", "9th", "10th")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Decile", y = "Frequency") +
  ggtitle("Income") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Indice = c("Over 12", "No", "Under 12", "Not applicable"),
  Valore = as.vector(table(dati$child)/nrow(ESS10))
)
ordine_desiderato <- c("Under 12", "Over 12", "No", "Not applicable")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Children in the household") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 28),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))


data <- data.frame(
  Indice = c("Paid \n work", "Student", "Looking \n for job", "Not \n looking \n for job",
             "Sick/ \n disabled", "Retired", "Community \n military \n service",
             "House \n work", "Other"),
  Valore = as.vector(table(ESS10$mnactic)/nrow(ESS10))
)
ordine_desiderato <- c("Paid \n work", "Student", "Looking \n for job", "Not \n looking \n for job",
                       "Sick/ \n disabled", "Retired", "Community \n military \n service",
                       "House \n work", "Other")
data$Indice <- factor(data$Indice, levels = ordine_desiderato)
ggplot(data, aes(x = factor(Indice), y = Valore)) +
  geom_bar(stat = "identity", fill = "gray70") +
  labs(x = "Category", y = "Frequency") +
  ggtitle("Main activity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 22),  
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28, hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 1))

### Likert plots ####

levels_order <- c("Bulgaria", "Italy", "Finland")

# eps 1005 372

# Age
ess.finland$agea_grouped <- cut(ess.finland$agea, breaks = seq(0, 100, by = 5), right = FALSE)
ess.italy$agea_grouped <- cut(ess.italy$agea, breaks = seq(0, 100, by = 5), right = FALSE)
ess.bulgaria$agea_grouped <- cut(ess.bulgaria$agea, breaks = seq(0, 100, by = 5), right = FALSE)
ess.finland$agea_grouped=as.factor(ess.finland$agea_grouped)
ess.italy$agea_grouped=as.factor(ess.italy$agea_grouped)
ess.bulgaria$agea_grouped=as.factor(ess.bulgaria$agea_grouped)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$agea_grouped, ess.finland$country), 
                          cbind(ess.italy$agea_grouped, ess.italy$country),
                          cbind(ess.bulgaria$agea_grouped, ess.bulgaria$country)))
colnames(new.data)=c("Age", "Country")
str(new.data)
new.data$Age=as.factor(new.data$Age)
levels(new.data$Age)=c("[15-20)","[20-25)","[25-30)","[30-35)","[35-40)","[40-45)","[45-50)","[50-55)",
                       "[55-60)","[60-65)","[65-70)","[70-75)","[75-80)","[80-85)","[85-90)","[90-95)")
"Age"=new.data$Age
age <- likert(data.frame(`Age`), grouping=new.data$Country)
plot(age)


# Health
ess.finland$health=as.factor(ess.finland$health)
ess.italy$health=as.factor(ess.italy$health)
ess.bulgaria$health=as.factor(ess.bulgaria$health)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$health, ess.finland$country), 
                          cbind(ess.italy$health, ess.italy$country),
                          cbind(ess.bulgaria$health, ess.bulgaria$country)))
colnames(new.data)=c("Health", "Country")
str(new.data)
new.data$Health=as.factor(new.data$Health)
levels(new.data$Health)=c("Very good", "Good",	"Fair", "Bad", "Very bad")
new.data$Health <- factor(new.data$Health, 
                          levels = rev(levels(new.data$Health)), 
                          ordered = TRUE)
"Health"=new.data$Health
new.data$Country <- factor(new.data$Country, levels = levels_order)
Health <- likert(data.frame(`Health`), grouping=new.data$Country)
plot(Health, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Hampered
ess.finland$hlthhmp=as.factor(ess.finland$hlthhmp)
ess.italy$hlthhmp=as.factor(ess.italy$hlthhmp)
ess.bulgaria$hlthhmp=as.factor(ess.bulgaria$hlthhmp)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$hlthhmp, ess.finland$country), 
                          cbind(ess.italy$hlthhmp, ess.italy$country),
                          cbind(ess.bulgaria$hlthhmp, ess.bulgaria$country)))
colnames(new.data)=c("Hampered in daily activities", "Country")
str(new.data)
new.data$`Hampered in daily activities`=as.factor(new.data$`Hampered in daily activities`)
levels(new.data$`Hampered in daily activities`)=c("Yes a lot", "Yes to some extent", "No")
"Hampered in daily activities"=new.data$`Hampered in daily activities`
new.data$Country <- factor(new.data$Country, levels = levels_order)
hamp <- likert(data.frame(`Hampered in daily activities`), grouping=new.data$Country)
plot(hamp, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Born in country
ess.finland$brncntr=as.factor(ess.finland$brncntr)
ess.italy$brncntr=as.factor(ess.italy$brncntr)
ess.bulgaria$brncntr=as.factor(ess.bulgaria$brncntr)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$brncntr, ess.finland$country), 
                          cbind(ess.italy$brncntr, ess.italy$country),
                          cbind(ess.bulgaria$brncntr, ess.bulgaria$country)))
colnames(new.data)=c("Born in the country", "Country")
str(new.data)
new.data$`Born in the country`=as.factor(new.data$`Born in the country`)
levels(new.data$`Born in the country`)=c("Yes","No")
new.data$`Born in the country` <- factor(new.data$`Born in the country`, 
                          levels = rev(levels(new.data$`Born in the country`)), 
                          ordered = TRUE)
"Born in the country"=new.data$`Born in the country`
new.data$Country <- factor(new.data$Country, levels = levels_order)
born <- likert(data.frame(`Born in the country`), grouping=new.data$Country)
plot(born, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


# Gender
ess.finland$gndr=as.factor(ess.finland$gndr)
ess.italy$gndr=as.factor(ess.italy$gndr)
ess.bulgaria$gndr=as.factor(ess.bulgaria$gndr)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$gndr, ess.finland$country), 
                          cbind(ess.italy$gndr, ess.italy$country),
                          cbind(ess.bulgaria$gndr, ess.bulgaria$country)))
colnames(new.data)=c("Gender", "Country")
str(new.data)
new.data$`Gender`=as.factor(new.data$`Gender`)
levels(new.data$`Gender`)=c("Male", "Female")
"Gender"=new.data$`Gender`
new.data$Country <- factor(new.data$Country, levels = levels_order)
gender <- likert(data.frame(`Gender`), grouping=new.data$Country)
plot(gender, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Education
table(ess.finland$eisced)
table(ess.italy$eisced)
table(ess.bulgaria$eisced)

educ.fin=ess.finland$eisced[ess.finland$eisced!="55"]
educ.ita=ess.italy$eisced[ess.italy$eisced!="55"&ess.italy$eisced!="3"]

educ.fin=as.factor(educ.fin)
educ.ita=as.factor(educ.ita)
ess.bulgaria$eisced=as.factor(ess.bulgaria$eisced)

country.fin=rep("Finland",length(educ.fin))
country.ita=rep("Italy",length(educ.ita))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(educ.fin, country.fin), 
                          cbind(educ.ita, country.ita),
                          cbind(ess.bulgaria$eisced, ess.bulgaria$country)))
colnames(new.data)=c("Education", "Country")
str(new.data)
new.data$`Education`=as.factor(new.data$`Education`)
levels(new.data$`Education`)=c("< lower secondary","Lower secondary",
                               "Upper tier upper secondary",
                               "Sub-degree", "BA level", ">= MA")
"Education"=new.data$`Education`
new.data$Country <- factor(new.data$Country, levels = levels_order)
educ <- likert(data.frame(`Education`), grouping=new.data$Country)
plot(educ, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


# Partner education
table(ess.finland$eiscedp)
table(ess.italy$eiscedp)
table(ess.bulgaria$eiscedp)

educp.fin=ess.finland$eiscedp[ess.finland$eiscedp!="55"]
educp.ita=ess.italy$eiscedp[ess.italy$eiscedp!="55"&ess.italy$eiscedp!="3"]

educp.fin=as.factor(educp.fin)
educp.ita=as.factor(educp.ita)
ess.bulgaria$eiscedp=as.factor(ess.bulgaria$eiscedp)

country.fin=rep("Finland",length(educp.fin))
country.ita=rep("Italy",length(educp.ita))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(educp.fin, country.fin), 
                          cbind(educp.ita, country.ita),
                          cbind(ess.bulgaria$eiscedp, ess.bulgaria$country)))
colnames(new.data)=c("Partner education", "Country")
str(new.data)
new.data$`Partner education`=as.factor(new.data$`Partner education`)
levels(new.data$`Partner education`)=c("< lower secondary","Lower secondary",
                               "Upper tier upper secondary",
                               "Sub-degree", "BA level", ">= MA", "NA")
`Partner education`=new.data$`Partner education`
new.data$Country <- factor(new.data$Country, levels = levels_order)
peduc <- likert(data.frame(`Partner education`), grouping=new.data$Country)
plot(peduc, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Income
ess.finland$hinctnta=as.factor(ess.finland$hinctnta)
ess.italy$hinctnta=as.factor(ess.italy$hinctnta)
ess.bulgaria$hinctnta=as.factor(ess.bulgaria$hinctnta)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$hinctnta, ess.finland$country), 
                          cbind(ess.italy$hinctnta, ess.italy$country),
                          cbind(ess.bulgaria$hinctnta, ess.bulgaria$country)))
colnames(new.data)=c("Income", "Country")
str(new.data)
new.data$`Income`=as.factor(new.data$`Income`)
levels(new.data$`Income`)=c("1st decile", "2nd decile", "3rd decile", "4th decile", "5th decile",
                            "6th decile", "7th decile", "8th decile", "9th decile", "10th decile")
"Income"=new.data$`Income`
new.data$Country <- factor(new.data$Country, levels = levels_order)
income <- likert(data.frame(`Income`), grouping=new.data$Country)
plot(income, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Children in the household
ess.finland$child=as.factor(dati.fin$child)
ess.italy$child=as.factor(dati.ita$child)
ess.bulgaria$child=as.factor(dati.bulg$child)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(as.factor(ess.finland$child), ess.finland$country), 
                          cbind(as.factor(ess.italy$child), ess.italy$country),
                          cbind(as.factor(ess.bulgaria$child), ess.bulgaria$country)))
colnames(new.data)=c("Children in the household", "Country")
str(new.data)
new.data$`Children in the household`=as.factor(new.data$`Children in the household`)
levels(new.data$`Children in the household`)=c("Over 12", "No", "Under 12", "Not applicable")
"Children in the household"=new.data$`Children in the household`
new.data$Country <- factor(new.data$Country, levels = levels_order)
child <- likert(data.frame(`Children in the household`), grouping=new.data$Country)
plot(child, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))

# Main activity
ess.finland$mnactic=as.factor(ess.finland$mnactic)
ess.italy$mnactic=as.factor(ess.italy$mnactic)
ess.bulgaria$mnactic=as.factor(ess.bulgaria$mnactic)
ess.finland$country=rep("Finland",nrow(ess.finland))
ess.italy$country=rep("Italy",nrow(ess.italy))
ess.bulgaria$country=rep("Bulgaria",nrow(ess.bulgaria))
new.data=data.frame(rbind(cbind(ess.finland$mnactic, ess.finland$country), 
                          cbind(ess.italy$mnactic, ess.italy$country),
                          cbind(ess.bulgaria$mnactic, ess.bulgaria$country)))
colnames(new.data)=c("Main activity", "Country")
str(new.data)
new.data$`Main activity`=as.factor(new.data$`Main activity`)
levels(new.data$`Main activity`)=c("Paid  job", "Student", "Looking for job", "Not looking for job",
                                   "Sick/disabled", "Retired", "Community/military service",
                                   "Housework", "Other","No answer")
"Main activity"=new.data$`Main activity`
new.data$Country <- factor(new.data$Country, levels = levels_order)
empl <- likert(data.frame(`Main activity`), grouping=new.data$Country)
plot(empl, text.size=5) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text=element_text(size=15))


## DESI Map ####

europe <- map_data("world")
head(europe)
table(europe$region)
paesi_da_selezionare <- c("Spain", "Portugal", "France", "Germany", "Italy", "UK", "Ireland", "Netherlands", 
                          "Belgium", "Luxembourg", "Switzerland", "Austria", "Liechtenstein", "Denmark", "Norway", "Sweden",
                          "Finland", "Iceland", "Greece", "Cyprus", "Malta", "Bulgaria", "Romania", "Croatia", 
                          "Slovenia", "Hungary", "Slovakia", "Czech Republic", "Poland", "Lithuania", "Latvia", "Estonia",
                          "Montenegro", "North Macedonia","Albania","Serbia","Bosnia and Herzegovina","Kosovo")
paesi_selezionati <- europe %>%
  filter(region %in% paesi_da_selezionare)
dim(paesi_selezionati)
table(paesi_selezionati$region)
paesi_selezionati <- paesi_selezionati %>%
  arrange(region)
dim(paesi_selezionati)
head(paesi_selezionati)
table(paesi_selezionati$region)
paesi_selezionati$colore <- as.numeric(factor(c(rep(0,length(which(paesi_selezionati$region=="Albania"))),
                                                rep(0,length(which(paesi_selezionati$region=="Austria"))),
                                                rep(0.50,length(which(paesi_selezionati$region=="Belgium"))),
                                                rep(0,length(which(paesi_selezionati$region=="Bosnia and Herzegovina"))),
                                                rep(0.38,length(which(paesi_selezionati$region=="Bulgaria"))),
                                                rep(0.48,length(which(paesi_selezionati$region=="Croatia"))),
                                                rep(0,length(which(paesi_selezionati$region=="Cyprus"))),
                                                rep(0.49,length(which(paesi_selezionati$region=="Czech Republic"))),
                                                rep(0,length(which(paesi_selezionati$region=="Denmark"))),
                                                rep(0.57,length(which(paesi_selezionati$region=="Estonia"))),
                                                rep(0.70,length(which(paesi_selezionati$region=="Finland"))),
                                                rep(0.53,length(which(paesi_selezionati$region=="France"))),
                                                rep(0,length(which(paesi_selezionati$region=="Germany"))),
                                                rep(0,length(which(paesi_selezionati$region=="Greece"))),
                                                rep(0.44,length(which(paesi_selezionati$region=="Hungary"))),
                                                rep("1",length(which(paesi_selezionati$region=="Iceland"))),
                                                rep(0.63,length(which(paesi_selezionati$region=="Ireland"))),
                                                rep(0.48,length(which(paesi_selezionati$region=="Italy"))),
                                                rep(0,length(which(paesi_selezionati$region=="Kosovo"))),
                                                rep(0,length(which(paesi_selezionati$region=="Latvia"))),
                                                rep(0,length(which(paesi_selezionati$region=="Liechtenstein"))),
                                                rep(0.53,length(which(paesi_selezionati$region=="Lithuania"))),
                                                rep(0,length(which(paesi_selezionati$region=="Luxembourg"))),
                                                rep(0,length(which(paesi_selezionati$region=="Malta"))),
                                                rep(0,length(which(paesi_selezionati$region=="Montenegro"))),
                                                rep(0.67,length(which(paesi_selezionati$region=="Netherlands"))),
                                                rep(0,length(which(paesi_selezionati$region=="North Macedonia"))),
                                                rep(0.64,length(which(paesi_selezionati$region=="Norway"))),
                                                rep(0,length(which(paesi_selezionati$region=="Poland"))),
                                                rep(0.48,length(which(paesi_selezionati$region=="Portugal"))),
                                                rep(0,length(which(paesi_selezionati$region=="Romania"))),
                                                rep(0,length(which(paesi_selezionati$region=="Serbia"))),
                                                rep(0.43,length(which(paesi_selezionati$region=="Slovakia"))),
                                                rep(0.53,length(which(paesi_selezionati$region=="Slovenia"))),
                                                rep(0,length(which(paesi_selezionati$region=="Spain"))),
                                                rep(0,length(which(paesi_selezionati$region=="Sweden"))),
                                                rep(0,length(which(paesi_selezionati$region=="Switzerland"))),
                                                rep(0.60,length(which(paesi_selezionati$region=="UK")))
)))

custom_breaks=seq(from = 1, to = 15, length.out = 4)
custom_labels=c("No data","0.40","0.55","0.70")

custom_breaks=seq(from = 1, to = 15, length.out = 4)
custom_labels=c("No data","0.30","0.50","0.70")

ggplot(paesi_selezionati, aes(x = long, y = lat, group = group, fill = colore)) +
  geom_polygon(color = "black", size = 0.5) +
  coord_fixed(1.3) +
  # scale_fill_manual(values = c("1" = "dimgray", "2" = "darkgray","No data" = "lightgray"), name = "Group") +
  scale_fill_gradient(low = "white", high = "dimgray", name = "DESI", 
                      breaks = custom_breaks, labels = custom_labels) +
  theme_minimal() +
  labs(title = "", x = "", y = "") + 
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_blank(),   
        axis.text.y = element_blank(),   
        axis.title.x = element_blank(),  
        axis.title.y = element_blank())  


paesi_da_evidenziare <- c("Italy","Finland","Bulgaria")
paesi_evidenziati <- europe %>%
  filter(region %in% paesi_da_evidenziare)
dim(paesi_evidenziati)
table(paesi_evidenziati$region)
paesi_evidenziati <- paesi_evidenziati %>%
  arrange(region)

paesi_evidenziati$colore <- as.numeric(factor(c(rep(0.38,length(which(paesi_selezionati$region=="Bulgaria"))),
                                                rep(0.70,length(which(paesi_selezionati$region=="Finland"))),
                                                rep(0.48,length(which(paesi_selezionati$region=="Italy")))
                                                ))
)

colors <- c("white", "#FEE090", "#1A9850")

ggplot(paesi_selezionati, aes(x = long, y = lat, group = group, fill = colore)) +
  geom_polygon(color = "black", size = 0.5) +   # Poligoni con confini neri
  geom_path(data = paesi_evidenziati, aes(x = long, y = lat, group = group), 
            color = "#b22222", size = 1) +        # Confini evidenziati in rosso
  coord_fixed(1.3) +
  # scale_fill_gradient(low = "white", high = "dimgray", name = "DESI index", 
                      # breaks = custom_breaks, labels = custom_labels) +
  scale_fill_gradientn(colors = colors, name = "DESI", 
                       breaks = custom_breaks, labels = custom_labels) +
  theme_minimal() +
  labs(title = "", x = "", y = "") + 
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_blank(),   
        axis.text.y = element_blank(),   
        axis.title.x = element_blank(),  
        axis.title.y = element_blank())


## Model estimation ####

G <- 2:4 

mod <- mlta(X=m, DM=DesMat, G = G, D = 1:3, wfix = FALSE, nstarts=10) 
mod$BIC 

mod.wfix <- mlta(X=m, DM=DesMat, G = G, D = 1:3, wfix = TRUE, nstarts=10) 
mod.wfix$BIC 

LTA <- lta(X=m, D=1:3, nstarts=10) # G=1
LTA$BIC

LCA <- lca(X=m, DM=DesMat, G=G, nstarts=10) # D=0
LCA$BIC


# choose the best model
MOD <- mod.wfix[[1]]

MOD$eta
MOD$z
MOD$b
MOD$w
MOD$beta

class <- apply(MOD$z,1,which.max) # clustering
class
table(class)


## Boostrap standard errors ####

N <- nrow(m)
M <- ncol(m)
J <- ncol(DesMat)
S <- 1000
G=3 # choose G for the selected model
D=1

res <- foreach (s=1:S, .packages = c("MASS","igraph")) %dopar% {
  n.new <- sample(V(g)[1:N], size=N, replace=TRUE)
  m.new <- m[n.new, ]
  mod <- mlta(X=m.new, DM=DesMat, G = G, D = D, wfix = FALSE, beta0=MOD$beta, nstarts=1)
}


b1 <- matrix(NA, nrow=S , ncol=M)
b2 <- matrix(NA, nrow=S , ncol=M)
b3 <- matrix(NA, nrow=S , ncol=M)
for (s in 1:S) {
  ord = order(apply(res[[s]]$b, 1, mean))
  bb = res[[s]]$b[ord,]
  b1[s,] <- as.vector(bb[1,])
  b2[s,] <- as.vector(bb[2,])
  b3[s,] <- as.vector(bb[3,])
}


w1 <- matrix(NA, nrow=S , ncol=M)
w2 <- matrix(NA, nrow=S , ncol=M)
w3 <- matrix(NA, nrow=S , ncol=M)
for (s in 1:S) {
  # ord = order(apply(res[[s]]$b, 1, mean)) # for varying wgk across groups
  # ww = res[[s]]$w[,,ord]
  # w1[s,] <- as.vector(ww[,1])
  # w2[s,] <- as.vector(ww[,2])
  # w3[s,] <- as.vector(ww[,3])
  w1[s,] <- as.vector(res[[s]]$w)
}


J=25
beta2 <- matrix(NA, nrow=S, ncol=J)
beta3 <- matrix(NA, nrow=S, ncol=J)
for (s in 1:S) {
  b.sim = res[[s]]$b
  ord = order(apply(res[[s]]$b, 1, mean))
  be = cbind(rep(0,J), res[[s]]$beta)
  be = be[,ord]
  beta2[s,] <- be[,2]
  beta3[s,] <- be[,3]
}


se.b1 <- apply(b1,2,sd)
se.b2 <- apply(b2,2,sd)
se.b3 <- apply(b3,2,sd)


se.w1 <- apply(abs(w1),2,sd)
se.w2 <- apply(abs(w2),2,sd)
se.w3 <- apply(abs(w3),2,sd)


se.beta2 <- apply(beta2, 2, sd)
se.beta3 <- apply(beta3, 2, sd)

z <- qnorm(0.975)

order(apply(MOD$b, 1, mean))


# plot b estimates and confidence intervals                                  
b_1 <- as.vector(MOD$b[1,])
u.b1 <- b_1 + z*se.b1
l.b1 <- b_1 - z*se.b1
y <- 1:7
par(mar=c(5.1,15,4.1,0), xpd=T)
# par(mar=c(5.1,2.1,4.1,0), xpd=T) # comment
plotCI(y=y, x=b_1, ui=u.b1, li=l.b1, xlab = "", ylab="", pch = 19, col="gray40", main = "Group 1", yaxt="n",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, cex=1, xlim=c(-11,11), lwd=2, err="x", frame.plot=F)
par(mar=c(5.1,15,4.1,0), xpd=F)
# par(mar=c(5.1,2.1,4.1,0), xpd=F) # comment
grid(lty="solid")
plotCI(y=y, x=b_1, ui=u.b1, li=l.b1, xlab = "", ylab="", pch = 19, col="gray40", main = "Group 1", yaxt="n",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, cex=1, xlim=c(-11,11), lwd=2, err="x", add=TRUE)
axis(2, at=1:7, labels=c("Internet", "Preference setting", "Advanced search", "PDF",
                         "Video calls", "Messages", "Posts"), las=1, cex=7, cex.axis=2)
abline(v=0, lty=2)
mtext(expression("b"),side=1,las=1,line=3, cex=2)


# plot w estimates and confidence intervals     
# w_1 <- as.vector(MOD$w[,,3]) # for varying wgk across groups
w_1 <- as.vector(MOD$w)
u.w1 <- w_1 + z*se.w1
l.w1 <- w_1 - z*se.w1
y <- 1:7
par(mar=c(5.1,15,4.1,0), xpd=T)
# par(mar=c(5.1,2.1,4.1,0), xpd=T) # comment
plotCI(y=y, x=w_1, ui=u.w1, li=l.w1, xlab = "", ylab="", pch = 19, col="gray40", main = "", yaxt="n",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, cex=1, xlim=c(-3,3), lwd=2, err="x", frame.plot=F)
par(mar=c(5.1,15,4.1,0), xpd=F)
# par(mar=c(5.1,2.1,4.1,0), xpd=F) # comment
grid(lty="solid")
plotCI(y=y, x=w_1, ui=u.w1, li=l.w1, xlab = "", ylab="", pch = 19, col="gray40", main = "", yaxt="n",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, cex=1, xlim=c(-3,3), lwd=2, err="x", add=TRUE)
axis(2, at=1:7, labels=c("Internet", "Preference setting", "Advanced search", "PDF",
                         "Video calls", "Messages", "Posts"), las=1, cex=7, cex.axis=2)
abline(v=0, lty=2)
mtext(expression("w"),side=1,las=1,line=3, cex=2)


# plot beta estimates and confidence intervals     

ord = order(apply(MOD$b, 1, mean))
be = cbind(0, MOD$beta)
be = be[,ord]

beta_2 <- as.vector(be[,2])
u.beta2 <- beta_2 + z*(se.beta2)
l.beta2 <- beta_2 - z*(se.beta2)
y <- 1:25
par(mfrow=c(1,1))
par(mar=c(4.3,17,2.1,0), xpd = T)
# par(mar=c(5.1,2.1,4.1,0), xpd=T) # comment
plotCI(y=y, x=beta_2, ui=u.beta2, li=l.beta2, xlab = "", ylab = "", pch = 19, col="gray40", main = "Group 2 vs 1",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, cex=1, xlim=c(-15,15), lwd=2, yaxt="n", err="x", frame.plot=F)
par(mar=c(4.3,17,2.1,0), xpd = F)
# par(mar=c(5.1,2.1,4.1,0), xpd=F) # comment
grid(lty="solid")
plotCI(y=y, x=beta_2, ui=u.beta2, li=l.beta2, xlab = "", ylab = "", pch = 19, col="gray40", main = "Group 2 vs 1",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, cex=1, xlim=c(-15,15), lwd=2, add=TRUE, err="x")
axis(2, at=1:25, labels=c("Intercept", "Age 30-50 (<30)", "Age 50-65 (<30)", "Age 65-75 (<30)",
                          "Age 75+ (<30)", "Fair health (bad)", "Good health (bad)",
                          "Hampered", "Born in country", "Male", "Low educ. (high)", "Medium educ. (high)",
                          "In a relationship", "Partner educ low (high)", "Partner educ med. (high)",
                          "Partner educ not appl. (high)", "Income low (high)", "Income medium (high)",
                          "No children (children 12+)", "Children <12 (children 12+)", "Not appl. (children 12+)",
                          "Retired (employed)", "Unemployed (employed)", "Age gap >2 years (same)",
                          "Age gap <2 years (same)"), las=1, cex.axis=1.5)
abline(v=0, lty=2)
mtext(expression(beta),side=1,las=1,line=3.5, cex=1.5)


## eta ####

jj=c(2,3,4,5,6,7,11,12,17,18,19,20,21)
gg=c(1,2)
gg=1
# MOD=mod.wfix
# MOD=mod[[1]]
# MOD$beta=cbind(rep(0,nrow(MOD$beta)),MOD$beta)
ris=matrix(NA,max(jj),length(gg))
for(j in jj){
  for(g in gg){
    print(j)
    ris[j,g]=exp(MOD$beta[j,g])/(1+sum(exp(MOD$beta[j,])))
  }
}
round(ris,2)


jj=c(2,3,4,5,6,7,11,12,17,18,19,20,21)
gg=c(1,2)
# gg=1
# MOD=mod.wfix
# MOD=mod[[1]]
# MOD$beta=cbind(rep(0,nrow(MOD$beta)),MOD$beta)
ris=matrix(NA,max(jj),length(gg))
for(j in jj){
  for(g in gg){
    print(j)
    ris[j,g]=(1)/(1+sum(exp(MOD$beta[j,])))
  }
}
round(ris,2)


jj=c(2,3,4,5)
gg=c(1,2)
for(g in gg){
  print(round(exp(1-sum(MOD$beta[jj,g]))/(1+sum(exp(1-colSums(MOD$beta[jj,])))),2))
}
jj=c(6,7)
gg=c(1,2)
for(g in gg){
  print(round(exp(1-sum(MOD$beta[jj,g]))/(1+sum(exp(1-colSums(MOD$beta[jj,])))),2))
}
jj=c(11,12)
gg=c(1,2)
for(g in gg){
  print(round(exp(1-sum(MOD$beta[jj,g]))/(1+sum(exp(1-colSums(MOD$beta[jj,])))),2))
}
jj=c(17,18)
gg=c(1,2)
for(g in gg){
  print(round(exp(1-sum(MOD$beta[jj,g]))/(1+sum(exp(1-colSums(MOD$beta[jj,])))),2))
}
jj=c(20)
for(g in gg){
  print(round(exp(1-(MOD$beta[jj,g]))/(1+sum(exp(1-(MOD$beta[jj,])))),2))
}


## Plot predicted probabilities ####

N=dim(MOD$mu)[1]
# G=2
# M=7
PI=array(0,c(G,M,N))
for (g in 1:G) {
  for(m in 1:M) {
    # pi=1/(1+exp(-(MOD$b[g,m]+MOD$w[1,m,g]*(MOD$mu[,1,g]))))
    # pi=1/(1+exp(-(MOD$b[g,m]+MOD$w[1,m]*(MOD$mu[,1,g])))) # bulgaria
    PI[g,m,]=pi
  }
}
order(apply(MOD$b, 1, mean))
PI11=PI[1,1,]
PI12=PI[1,2,]
PI13=PI[1,3,]
PI14=PI[1,4,]
PI15=PI[1,5,]
PI16=PI[1,6,]
PI17=PI[1,7,]
PI21=PI[2,1,]
PI22=PI[2,2,]
PI23=PI[2,3,]
PI24=PI[2,4,]
PI25=PI[2,5,]
PI26=PI[2,6,]
PI27=PI[2,7,]
PI31=PI[3,1,]
PI32=PI[3,2,]
PI33=PI[3,3,]
PI34=PI[3,4,]
PI35=PI[3,5,]
PI36=PI[3,6,]
PI37=PI[3,7,]

data <- data.frame(
  Value = c(PI11, PI12, PI13, PI14, PI15, PI16, PI17,
            PI21, PI22, PI23, PI24, PI25, PI26, PI27,
            PI31, PI32, PI33, PI34, PI35, PI36, PI37),
  Group = rep(c("g=1, Internet", "g=1, Pref. settings", "g=1, Adv. search", "g=1, PDFs", "g=1, Video calls", "g=1, Messages", "g=1, Online posts",
                "g=2, Internet", "g=2, Pref. settings", "g=2, Adv. search", "g=2, PDFs", "g=2, Video calls", "g=2, Messages", "g=2, Online posts",
                "g=3, Internet", "g=3, Pref. settings", "g=3, Adv. search", "g=3, PDFs", "g=3, Video calls", "g=3, Messages", "g=3, Online posts"), 
              each = length(PI11))
)

means <- aggregate(data$Value, list(Group = data$Group), mean)
data$ord <- seq_len(nrow(data))

data <- merge(data, means, by = "Group")
data <- data[order(data$ord), ]
colnames(data)[ncol(data)] <- "Mean"
head(data)

data$Group <- factor(data$Group, levels = unique(data$Group))

ggplot(data, aes(x = Value, fill = Group)) +
  geom_histogram(binwidth = 0.08, color = "gray40") +
  geom_vline(data = data, aes(xintercept = Mean), color = "firebrick", linetype = "dashed", size = 1) +
  facet_wrap(~Group, nrow = 3, scales = "free") +
  scale_fill_manual(values = rep("gray85", 21)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(strip.text = element_text(face = "italic", size = 15),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") + 
  coord_cartesian(xlim = c(0, 1))