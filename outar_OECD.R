library(BDgraph)
# library(geoR)
# library(expm)
library(MASS)
# library(mvtnorm)
# library(statmod)
# library(expm) 
library(pracma)
# library(fda)
# library(rstiefel)
library(parallel)
# library(ssgraph)
# library(CVglasso)
# library(spectralGP)
library(igraph)
library(combinat)
library(INLA)
library(Matrix)
# library(foreach)
# library(doParallel)
# library(doRNG)
library(astsa)
library(tictoc)


# dict = 2:8
# 
# args = commandArgs(trailingOnly = TRUE)
# argslen = length(args)
# if(argslen > 1) stop('Error: Too Many Arguments')
# if(argslen < 1) stop('Error: Takes One Argument Only')
# args1 = dict[as.numeric(args)]
# 
# 
# cat("\n Process order: ", args1, "\n")

OUTAR = function(Y, K, n1, n2, n.iter = 1, verbose = T, adaptA = T, adaptAR = T){
  # Y: p-variate time series data
  # K: K of the AR models fit to the latent components
  # verbose: If TRUE, model displays the number of iterations, acceptance probability and RMSE of the estimator of Omega
  # adapt: If TRUE, the adaptive covariance function of Haario et al. (2001) is used to evaluate the covariance matrix of the proposal.
  # n.iter: Number of Monte Carlo iterations
  
  Y = as.matrix(Y)
  # p: dimension of the time series data
  p = nrow(Y)
  # N: number of observations in each time series
  N = ncol(Y)
  # print(dim(Y))
  cat("\n Processing OUTAR model with ", p, " time series with ", n1+n2, " observations in each...\n")
  # Store the parameters. Note that the matrices are stored using their vectorized representations.
  A.samples = array(NA, dim=c(n.iter, p*(p-1)/2))
  AZ.samples = array(NA, dim=c(n.iter, p*(p-1)/2))
  D.samples = array(NA, dim=c(n.iter, p))
  L1.samples = L2.samples = array(NA, dim=c(n.iter, p*(p-1)/2))
  LZ1.samples = LZ2.samples = array(NA, dim=c(n.iter, p*(p-1)/2))
  PHI.samples = array(NA, dim=c(n.iter, p, K))
  SIGMA.samples = array(NA, dim=c(n.iter, p))
  Omega1.samples = array(NA, dim=c(n.iter, p, p))
  Omega2.samples = array(NA, dim=c(n.iter, p, p))
  PACF.samples = array(NA, dim=c(n.iter, p, K))
  Y_onestep_exp = Y_onestep_raw = matrix(NA, nrow=n.iter, ncol=p)
  # LASSO = CVglasso(X = t(Y), lam.min.ratio = 1e-7)[["Omega"]]
  # cat("\n RMSE by LASSO: ", RMSE_LASSO, "\n")
  # success = FALSE
  # while (!success) {
  #   # print(1)
  #   tryCatch({
  #     GGM = bdgraph( data = t(Y), method="ggm", iter = 10000 )
  #     success = TRUE  # If no error, set success to TRUE
  #   }, error = function(e) {
  #     # message("Error in AR simulation: ", e$message, ". Retrying...")
  #     success = FALSE  # Continue the loop to retry
  #   })
  # }
  GGM1 = bdgraph( data = t(Y[, 1:n1]), method="ggm", iter = 10000 )
  GGM2 = bdgraph( data = t(Y[, n1+1:n2]), method="ggm", iter = 10000 )
  # RMSE_GGM = sqrt(mean(c(Omega0 - GGM$K_hat)^2))
  # cat("\n RMSE by GGM: ", RMSE_GGM, "\n")
  # GCGM = bdgraph( data = t(Y), method="gcgm", iter = 10000 )
  # RMSE_GCGM = sqrt(mean(c(Omega0 - GCGM$K_hat)^2))
  # cat("\n RMSE by GCGM: ", RMSE_GCGM, "\n")
  
  # Initialize L, D and A
  Omega1 = GGM1$K_hat
  Omega2 = GGM2$K_hat
  cat("\n Running process with order = ", K, "\n")
  L1.chol = t(chol(Omega1))
  L2.chol = t(chol(Omega2))
  # print(L.chol)
  # D.samples[1, ] = diag(IL)
  IL1 = L1.chol %*% diag( 1 / diag(L1.chol) )
  IL2 = L2.chol %*% diag( 1 / diag(L2.chol) )
  # L.samples[1, ] = L[lower.tri(L)]
  # A.samples[1, ] = rep(0, p*(p-1)/2)
  # print(L.samples[1,])
  # print(D.samples[1,])
  # L_test = L[1:5, 1:5]
  # print(L_test)
  # print(L_test[lower.tri(L_test)])
  # LL = diag(5)
  # LL[lower.tri(LL)] = L_test[lower.tri(L_test)]
  # print(LL)
  
  # Initialize lambda parameters
  lambda_A_U = lambda_L_U = 1;
  lambda_A = lambda_L1 = lambda_L2 = 0
  # eps_L = 1e-3 / (p*(p-1)/2); eps_D = 1e-3 / p # initial step size
  eps_L = 1e-6; eps_D = 1e-6;
  sigma_D = 10
  sigma_A = 1e-2; sigma_L = 1e-2
  print(dim(Y))
  # compute Z from U, D L and Y
  # mat<> stores <> in matrix format
  matA = matrix(0, p, p)
  matA[lower.tri(matA)] = A = AZ = rep(0, p*(p-1)/2)
  matA = matA - t(matA)
  matU = (diag(p) - matA) %*% solve(diag(p) + matA)
  matD = D = diag(L1.chol) + diag(L2.chol)
  matD = diag(matD)
  matL1 = matL2 = matrix(0, p, p)
  LZ1 = - IL1[lower.tri(IL1)]
  LZ2 = - IL2[lower.tri(IL2)]
  L1 = LZ1 * (abs(LZ1) > lambda_L1)
  L2 = LZ2 * (abs(LZ2) > lambda_L2)
  matL1[lower.tri(matL1)] = L1
  matL2[lower.tri(matL2)] = L2
  Z = t(matU) %*% matD %*% cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                                  (diag(p) - matL2) %*% Y[, n1+1:n2])
  # print("check")
  # print(dim(Z))
  PHI = PACF = matrix(0, p, K)
  SIGMA = rep(0, p)
  # initialize phis an sigmas based on Z
  for(i in 1:p){
    # init_model = arima(Z[i, ], order = c(K, 1, 0), transform.pars = TRUE, include.mean = FALSE, method = "ML")
    init_model = ar(Z[i, ], aic = FALSE, order.max = K)
    # PHI[i, ] = phi = init_model[["coef"]]
    PHI[i, ] = phi = init_model[["ar"]]
    PACF[i, ] = INLA::inla.ar.phi2pacf(phi)
    acf = INLA::inla.ar.phi2acf(phi)[-1]
    SIGMA[i] = sigma = sqrt(1 - sum(phi * acf))
    # print(c(phi, acf, sigma))
  }
  # print(PHI)
  # cat("\n", SIGMA)
  # PHI = PHI.samples[1,,]
  # SIGMA = SIGMA.samples[1,]
  # cat("\n check 1")
  cat("\n")
  
  ae_A = 0.1 + p * (p-1)/2
  be_A = 0.1 + sum(AZ ^ 2)/2
  
  sigma_T_A = 1/rgamma(1, ae_A, be_A)
  
  ae_L1 = 0.1 + p * (p-1)/2
  be_L1 = 0.1 + sum(LZ1 ^ 2)/2
  
  sigma_T_L1 = 1/rgamma(1, ae_L1, be_L1)
  
  ae_L2 = 0.1 + p * (p-1)/2
  be_L2 = 0.1 + sum(LZ2 ^ 2)/2
  
  sigma_T_L2 = 1/rgamma(1, ae_L2, be_L2)
  
  # Joint log-likelihood of a single AR(K) processes ----- (1)
  # Note: Y has L, D and U
  ll_Zi = function(tms, phi, sigma){
    tsvals = embed(ts(tms), dimension = K + 1)
    return( - 0.5 * ( 2 * N * log(sigma) + sum((tsvals[, 1] - tsvals[, -1] %*% phi)^2) /  sigma^2 ) )
  }
  # print(ll_Zi(Z[1,], PHI[1,], SIGMA[1]))
  # Wrapper for the joint log-likelihood of the p AR(K) processes
  ll_Z = function(data, phi.values, sigma.values, ID = NULL){
    # numCores = detectCores() - 1
    # cl = makeCluster(numCores)
    # registerDoParallel(cl)
    # result = rep(NA, p)
    # results = foreach(i = 1:p, .combine = c, .export = c("N", "K", "ll_Zi")) %dopar% {
    # for(i in 1:p){
    #   result[i] = ll_Zi(tms=data[i, ], phi=phi.values[i,], sigma=sigma.values[i])
    # }
    
    if(is.null(ID))
      result = lapply(1:p, function(i) ll_Zi(tms=data[i, 1:n1], phi=phi.values[i,], sigma=sigma.values[i]) +
                  ll_Zi(tms=data[i, n1+1:n2], phi=phi.values[i,], sigma=sigma.values[i]))
    else
      result = lapply(1:p, function(i) ll_Zi(tms=data[i, ], phi=phi.values[i,], sigma=sigma.values[i]))
    # stopCluster(cl)
    # print(result)
    return( sum(unlist(result)) )
  }
  
  # print(ll_Z(Z, PHI, SIGMA))
  # Gradient of the negative log-likelihood of Z wrt Z
  gradient_ll_Zi = function(tms, phi, sigma) {
    # print(length(tms))
    N = length(tms)
    # Create the lagged matrix using the embed function
    tsvals = embed(tms, K + 1)
    residuals = tsvals[, 1] - tsvals[, -1] %*% phi
    
    # Initialize the derivative vector
    grad_Zi = numeric(N)
    
    # Calculate the derivative for time points from (p+1) to T
    for (t in (K+1):N) {
      grad_Zi[t] = -1 / sigma^2 * residuals[t - K]
      for (k in 1:K) {
        if ((t + k) <= N) {
          grad_Zi[t] = grad_Zi[t] + 1 / sigma^2 * residuals[t - K + k] * phi[k]
        }
      }
    }
    # print(grad_Zi)
    # Handle the first p components separately
    for (t in 1:K) {
      for (k in 1:t) {
        if ((t + k) <= N) {
          grad_Zi[t] = grad_Zi[t] + 1 / sigma^2 * residuals[k] * phi[K-t+k]
        }
      }
    }
    
    return(grad_Zi)
  }
  
  # print(gradient_ll_Zi(Z[1,], PHI[1,], SIGMA[1]))
  # print(pracma::grad(ll_Zi, Z[1,], phi=PHI[1,], sigma=SIGMA[1]))
  
  
  gradient_ll_Z = function(Z, phi.values, sigma.values, ID = NULL){
    # Number of cores to use
    # numCores = detectCores() - 1
    # cl = makeCluster(numCores)
    # registerDoParallel(cl)
    grad_Z = NULL
    # print(dim(Z))
    # Parallel computation of the gradient
    # grad_Z = foreach(i = 1:p, .combine = rbind, .packages = c("foreach"), .export = c("N", "K")) %dopar% {
    # for(i in 1:p){
    #   # cat("\n check 3")
    #   # }
    #   # print(grad_i)
    #   grad_Z = rbind(grad_Z, gradient_ll_Zi(tms=Z[i, ], phi=PHI[i, ], sigma=SIGMA[i]))
    # 
    #   # return(grad_i / sigma.values[i]^2)
    # }
    if(is.null(ID))
      grad_Z = lapply(1:p, function(i) c(gradient_ll_Zi(tms=Z[i, 1:n1], phi=PHI[i, ], sigma=SIGMA[i]),
                                         gradient_ll_Zi(tms=Z[i, n1 + 1:n2], phi=PHI[i, ], sigma=SIGMA[i])))
    else
      grad_Z = lapply(1:p, function(i) gradient_ll_Zi(tms=Z[i, ], phi=PHI[i, ], sigma=SIGMA[i]))
    
    grad_Z = do.call(rbind, grad_Z)
    # print(dim(grad_Z))
    # stopCluster(cl)
    # stopifnot(dim(grad_Z)==c(p, N))
    return(grad_Z)
  }
  # cat("\n check 4")
  # print(pracma::grad(f=ll_Z, x0=Z, phi.values=PHI, sigma.values=SIGMA))
  # print(c(gradient_ll_Z(Z, PHI, SIGMA)))
  # gradient_ll_Z(Z, PHI, SIGMA)
  # Joint log-likelihood of the prior put on A ----- (2)
  # log_prior_A = function()
  
  # Gradient of (2)
  
  # Joint log-likelihood of the prior put on L ----- (3)
  # ll_L = function(L.values, lambda, h0 = 1e-8){
  #   return( sum(L.values + log((1 + 2 * (pi)^(-1) * atan((L.values^2 - lambda^2)/h0))/2)) )
  # }
  
  # Gradient of (1) + (3) wrt L
  gradient_L = function(L.values, LZ.values, Y, Z, lambda, sigma_T, ID, h0=1e-8){
    # Return the gradient as a flattened vector
    # Contribution from (1)
    grad_ll = - kronecker(Y, matD %*% matU)
    # print(dim(grad_nll))
    # print(dim(grad_ll))
    # print(dim(Z))
    # print(dim(gradient_ll_Z(Z=Z, PHI, SIGMA)))
    grad_ll = matrix(grad_ll %*% c(gradient_ll_Z(Z=Z, PHI, SIGMA, ID)), p, p)
    grad_ll = grad_ll[lower.tri(grad_ll)]
    
    # Contribution from threshholding function
    grad_thrshld = 0.5 * (1 + 2 * (pi)^(-1) * atan((LZ.values^2 - lambda^2)/h0)) + 2 * LZ.values^2 / ( pi*h0 * ( 1 + (LZ.values^2 - lambda^2)^2 / h0^2 ) )
    
    return(grad_ll*grad_thrshld - LZ.values / sigma_T^2)
  }
  # return(gradient_L(L[lower.tri(L)], 0))
  
  # Joint log-likelihood of the prior put on log(d) ----- (4)
  # d has an inverse-gaussian distribution with mean mu and scale = 1
  # For the transformation d -> log(d), add the log-Jacobian term: log(d)
  ll_logD = function(D.values, mu){
    ll = sum( -(3/2 * log(D.values) + (D.values - mu)^2 / (2 * D.values)) ) # log-densities of d
    return( ll + sum(log(D.values)) ) # add the Jacobian adjustment
  }
  
  
  # Gradient of (1) + (4) wrt log(d)
  # The effective difference between gradients with current and proposed values is only due to the log-densities of d
  gradient_logD = function(D.values, mu, ID = NULL){
    # Contribution from (1)
    grad_ll = kronecker(cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                               (diag(p) - matL2) %*% Y[, n1+1:n2]), matU)
    # print(dim(grad_ll))
    # print(dim(gradient_ll_Z(Z, PHI, SIGMA))) 
    # print(dim(grad_nll))
    grad_ll = matrix(grad_ll %*% c(gradient_ll_Z(Z, PHI, SIGMA, ID)), p, p)
    grad_ll = diag(grad_ll) * D.values + c(-(D.values^2 - mu^2) / (2 * D.values))
    # grad = diag( grad_ll + c(-(D.values^2 - mu^2) / (2 * D.values)) )
    return(grad_ll)
  }
  # return(gradient_logD(diag(D), rnorm(1)))
  # Joint log-likelihood of the prior put on the pacfs ----- (5)
  # Let U be a pacf. We put a prior Uniform(-1, 1) on U
  # Conveniently, the prior is put on V = log((U+1)/2) ~ Logistic(0,1)
  # ll_PACF = function(V.values){
  #   return(sum(dlogis(V.values, mean=0, scale=1, log=TRUE)))
  # }
  
  # Gradient of (5)
  
  # Joint log-likelihood of the prior put on lambda ----- (6)
  
  # Gradient of (6)
  
  
  # Metropolis-adjusted Langevin Monte Carlo (MALA) implementations
  update_L = function(init_L, init_LZ, matL, Y, step_size, lambda, sigma_T, ar_counter, index=NULL){
    # print(dim(Y))
    # print(init_LZ)
    # vecL = init_L
    Z = t(matU) %*% matD %*% (diag(p) - matL) %*% Y
    log_den = ll_Z(Z, PHI, SIGMA, index) - 0.5 * sum(init_LZ^2) / sigma_T^2
    gradient = gradient_L(L.values=init_L, LZ.values=init_LZ, Y=Y, Z=Z, lambda=lambda, sigma_T=sigma_T, ID=index)
    proposal_mean = init_LZ + 0.5 * step_size * gradient
    proposal_LZ = proposal_mean + sqrt(step_size) * rnorm(p*(p-1)/2)
    # print(proposal_LZ)
    proposal_L = proposal_LZ * (abs(proposal_LZ) > lambda) 
    log_den = log_den - 0.5 * sum((proposal_LZ - proposal_mean)^2) / step_size
    
    matL_prop = matL
    matL_prop[lower.tri(matL_prop)] = proposal_L
    Z = t(matU) %*% matD %*% (diag(p) - matL_prop) %*% Y
    log_num = ll_Z(Z, PHI, SIGMA, index) - 0.5 * sum(proposal_LZ^2) / sigma_T^2
    proposal_gradient = gradient_L(L.values=proposal_L, LZ.values=proposal_LZ, Y=Y, Z=Z, lambda=lambda, sigma_T=sigma_T, ID=index)
    reverse_proposal_mean = proposal_LZ + 0.5 * step_size * proposal_gradient
    log_num = log_num - 0.5 * sum((init_LZ - reverse_proposal_mean)^2) / step_size
    log_accept_ratio = log_num - log_den
    
    if(is.na(log_accept_ratio) || is.nan(log_accept_ratio))
      log_accept_ratio = -Inf
    
    if (log(runif(1)) < log_accept_ratio) {
      init_L = proposal_L
      init_LZ = proposal_LZ
      matL = matL_prop
      ar_counter = ar_counter + 1
    }
    
    return(list(L = init_L, LZ = init_LZ, matL = matL, arc = ar_counter))
  }
  
  update_D = function(init_D, xi, step_size, ar_counter, index=NULL){
    
    init_logD = log(init_D)
    Z = t(matU) %*% diag(init_D) %*% cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                                            (diag(p) - matL2) %*% Y[, n1+1:n2])
    log_den = ll_Z(Z, PHI, SIGMA, index) + ll_logD(D.values=init_D, mu=xi)
    gradient = gradient_logD(D.values=init_D, mu=xi, ID = index)
    proposal_mean = init_logD + 0.5 * step_size * gradient
    proposal_logD = proposal_mean + sqrt(step_size) * rnorm(p)
    log_den = log_den - 0.5 * sum((proposal_logD - proposal_mean)^2) / step_size
    # print(log_den)
    proposal_D = exp(proposal_logD)
    Z = t(matU) %*% diag(proposal_D) %*% cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                                                (diag(p) - matL2) %*% Y[, n1+1:n2])
    log_num = ll_Z(Z, PHI, SIGMA, index) + ll_logD(D.values=proposal_D, mu=xi)
    # print(log_num)
    proposal_gradient = gradient_logD(D.values=proposal_D, mu=xi, ID = index)
    reverse_proposal_mean = proposal_logD + 0.5 * step_size * proposal_gradient
    log_num = log_num - 0.5 * sum((init_logD - reverse_proposal_mean)^2) / step_size
    # print(log_num)
    log_accept_ratio = log_num - log_den
    
    if(is.na(log_accept_ratio) || is.nan(log_accept_ratio))
      log_accept_ratio = -Inf
    
    if (log(runif(1)) < log_accept_ratio) {
      init_D = proposal_D
      matD = diag(proposal_D)
      ar_counter = ar_counter + 1
    }
    
    return(list(D = init_D, matD = matD, arc = ar_counter))
  }
  
  # MH Sampler for A
  update_A = function(init_A, init_AZ, lambda, sigma_T, proposal_sigma, ar_counter, index=NULL){
    
    # A0 = matrix(0, p, p)
    # A0[lower.tri(A0)] = init_A
    # init_U = (diag(p) - A0) %*% solve(diag(p) + A0)
    Z = t(matU) %*% matD %*% cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                                    (diag(p) - matL2) %*% Y[, n1+1:n2])
    log_posterior = ll_Z(Z, PHI, SIGMA, index) - sum(init_AZ^2) / (2 * sigma_T^2)
    proposal_AZ = init_AZ + MASS::mvrnorm(1, mu=rep(0, p*(p-1)/2), Sigma = proposal_sigma)
    proposal_A = (abs(proposal_AZ) > lambda) * (abs(proposal_AZ) - lambda) * sign(proposal_AZ)
    # cat("proposal_A", proposal_A[1:10], "\n")
    A_prop = matrix(0, p, p)
    A_prop[lower.tri(A_prop)] = proposal_A
    matU_prop = (diag(p) - A_prop) %*% solve(diag(p) + A_prop)
    Z = t(matU_prop) %*% matD %*% cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                                         (diag(p) - matL2) %*% Y[, n1+1:n2])
    log_proposed_posterior = ll_Z(Z, PHI, SIGMA, index) - sum(proposal_AZ^2) / (2 * sigma_T^2)
    log_accept_prob = log_proposed_posterior - log_posterior
    # print(acc_prob)
    # cat("\n", log_accept_prob)
    if(is.na(log_accept_prob) || is.nan(log_accept_prob))
      log_accept_prob = -Inf
    
    if(log(runif(1)) < log_accept_prob)
    {
      init_A = proposal_A
      init_AZ = proposal_AZ
      matU = matU_prop
      ar_counter = ar_counter + 1
    }
    return(list(A = init_A, AZ = init_AZ, matU = matU, arc = ar_counter))
  }
  
  # MH Sampler for AR processes
  update_AR = function(data, current_params, proposal_sigma, ar_counter){
    current_phi = current_params[1:K]
    # current_pacf = INLA::inla.ar.phi2pacf(current_phi)
    current_pacf = current_params[(K+1):(2*K)]
    # print(current_pacf)
    current_sigma = current_params[2*K+1]
    # flag = FALSE # indicate acceptance
    
    log_posterior = ll_Zi(tms=data[1:n1], phi=current_phi, sigma=current_sigma) + 
      ll_Zi(tms=data[n1+1:n2], phi=current_phi, sigma=current_sigma) + 
      sum(dlogis(qlogis(((current_pacf + 1) / 2)), location = 0, scale = 1, log = TRUE))
    
    logit_proposal = MASS::mvrnorm(1, mu = qlogis(((current_pacf + 1) / 2)), Sigma = proposal_sigma)   # V = logit((U + 1) / 2) ~ Logistic(0, 1)
    proposed_pacf = plogis(logit_proposal) * 2 - 1  # inverse_logit(V) * 2 - 1
    # print(proposed_pacf)
    proposed_phi = INLA::inla.ar.pacf2phi(proposed_pacf)
    proposed_acf = INLA::inla.ar.pacf2acf(proposed_pacf, lag.max = K)[-1]
    # if(iter %% 100 == 0){
    #   print(proposed_phi)
    #   print(proposed_acf)
    # }
    if(length(proposed_acf) == 0 || 0 %in% proposed_phi || 1 %in% proposed_acf ||
       1 %in% proposed_phi || -1 %in% proposed_acf)
      log_proposed_posterior = NaN
    else{
      # print(current_acf)
      stopifnot(length(proposed_phi) == length(proposed_acf))
      proposed_sigma = sqrt(1 - sum(proposed_phi * proposed_acf))
      log_proposed_posterior = ll_Zi(tms=data[1:n1], phi=proposed_phi, sigma=proposed_sigma) + 
        ll_Zi(tms=data[n1+1:n2], phi=proposed_phi, sigma=proposed_sigma) +
        sum(dlogis(logit_proposal, location = 0, scale = 1, log = TRUE))
    }
    log_accept_prob = log_proposed_posterior - log_posterior
    # print(acc_prob)
    if(is.na(log_accept_prob) || is.nan(log_accept_prob))
      log_accept_prob = -Inf
    
    if(log(runif(1)) < log_accept_prob)
    {
      current_phi = proposed_phi
      current_pacf = proposed_pacf
      ar_counter = ar_counter + 1
      current_sigma = proposed_sigma
    }
    return(c(current_phi, current_pacf, current_sigma, ar_counter))
  }
  
  
  # parameter initializations for the MCMC
  a_rate = a_rate_L1 = a_rate_L2 = a_rate_A = a_rate_D = 0
  sdA = 1e-10; sdAR = rep(1e-8, p)
  arcA = arcL1 = arcL2 = arcD = arcLU = arcLL1 = arcLL2 = 0; arcAR = rep(0, p)
  Cov_prop_AR = array(rep(diag(K), p), dim=c(K, K, p))
  Cov_prop_A = diag(p*(p-1)/2)
  
  # MCMC sampling. 10000 samples with 5000 burnin
  for(iter in 1:n.iter){
    # print(iter)
    # Draw sample for A (U)
    if(iter > 1500){
      result = update_A(init_A=A, init_AZ=AZ, lambda=lambda_A, sigma_T = sigma_T_A, 
                        proposal_sigma=sdA*Cov_prop_A, ar_counter=arcA, index = NULL)
      A.samples[iter,] = A = result$A
      AZ.samples[iter,] = AZ = result$AZ
      matU = result$matU
      arcA = result$arc
      # cat("\n check 1")  # ************** Checkpoint
      # cat("\n", A[30:50])
    }
    
    # if(iter > 1500 && iter %% 100 == 0)
    #   print(A)
    
    # Draw sample for L1
    if(iter > 1500){
      result = update_L(init_L=L1, init_LZ=LZ1, matL = matL1, Y = Y[, 1:n1], 
                        step_size=eps_L, lambda=lambda_L1, sigma_T = sigma_T_L1, 
                        ar_counter=arcL1, index = 1)
      L1.samples[iter,] = L1 = result$L
      LZ1.samples[iter,] = LZ1 = result$LZ
      matL1 = result$matL
      # L = matL[lower.tri(matL)]
      arcL1 = result$arc
      # cat("\n check 2")  # ************** Checkpoint
      # cat("\n", sum(L>10), "\t", eps_L)
    }
    
    # Draw sample for L2
    if(iter > 1500){
      result = update_L(init_L=L2, init_LZ=LZ2, matL = matL2, Y = Y[, n1+1:n2], 
                        step_size=eps_L, lambda=lambda_L2, sigma_T = sigma_T_L2,
                        ar_counter=arcL2, index = 2)
      L2.samples[iter,] = L2 = result$L
      LZ2.samples[iter,] = LZ2 = result$LZ
      matL2 = result$matL
      # L = matL[lower.tri(matL)]
      arcL2 = result$arc
      # cat("\n check 2")  # ************** Checkpoint
      # cat("\n", sum(L>10), "\t", eps_L)
    }
    
    # if(iter > 1500 && iter %% 100 == 0)
    #   print(LZ[1:20])
    
    # Draw sample for D
    if(iter > 1500){
      xi = rnorm(1, 0, sigma_D) # mean of the inverse gaussian prior
      result = update_D(init_D=D, step_size=eps_D, xi=xi, ar_counter=arcD, index = NULL)
      D.samples[iter,] = D = result$D
      matD = result$matD
      # D = diag(matD)
      arcD = result$arc
      # cat("\n check 3")  # ************** Checkpoint
      # cat("\n", sum(D > 15), "\t", eps_D)
    }
    
    # print(matD[1:5, 1:5])
    
    if(iter > 1500){
      Omega1.samples[iter, , ] = Omega_pred1 = (diag(p) - matL1) %*% (matD ^ 2) %*% t(diag(p) - matL1)
      Omega2.samples[iter, , ] = Omega_pred2 = (diag(p) - matL2) %*% (matD ^ 2) %*% t(diag(p) - matL2)
    }
    
    # if(iter > 1500 && iter%%100 == 0){ 
    #   print(c(Omega_pred1[1:5,1:5]))
    #   print(c(Omega_pred2[1:5,1:5]))
    # }
    
    # print(mean((Omega_pred-Omega0)^2))
    
    # Draw samples for AR process parameters
    # numCores = detectCores() - 1
    # cl = makeCluster(numCores)
    # registerDoParallel(cl)
    # registerDoRNG()
    Z = t(matU) %*% matD %*% cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                                    (diag(p) - matL2) %*% Y[, n1+1:n2])
    result = NULL
    # Use foreach with doRNG to ensure reproducibility
    # result = foreach(i = 1:p, .packages = 'doRNG') %dopar% {
    # for(i in 1:p){
    #   result = rbind(result, update_AR(data=Z[i,], current_params=c(PHI[i,], PACF[i, ], SIGMA[i]), proposal_sigma=sdAR[i]*Cov_prop_AR[,,i], ar_counter=arcAR[i]))
    # }
    result = mclapply(1:p, function(i) update_AR(data=Z[i, ], current_params=c(PHI[i,], PACF[i, ], SIGMA[i]), 
                                                 proposal_sigma=sdAR[i]*Cov_prop_AR[,,i], ar_counter=arcAR[i]), mc.cores = 1L)
    result = do.call(rbind, result)
    # stopCluster(cl)
    
    # extract results
    # print(result)
    # result = do.call(rbind, result)
    PHI.samples[iter, , ] = PHI = result[, 1:K]   # p x K matrix
    PACF.samples[iter, , ] = PACF = result[, (K+1):(2*K)]
    SIGMA.samples[iter, ] = SIGMA = result[, 2*K+1]
    arcAR = result[, 2*K+2]
    
    # One-step ahead predictions
    Y_onestep_exp[iter, ] = solve(diag(p) - matL2) %*% diag(1/diag(matD)) %*% matU %*% 
      sapply( 1:p, function(i) sum(Z[i, N:(N-K+1)] * PHI[i, ]), simplify = "vector" )
    Y_onestep_raw[iter, ] = solve(diag(p) - matL2) %*% diag(1/diag(matD)) %*% matU %*% 
      sapply( 1:p, function(i) sum(Z[i, N:(N-K+1)] * PHI[i, ]) + SIGMA[i], simplify = "vector" )
    
    # cat("\n", ll_Z(Z, PHI, SIGMA))
    # 
    # if(iter %% 100 == 0)
    #   cat("\n", sum(SIGMA>1))
    
    # Draw lambda_L1
    if(iter == 2500){
      LZ1_avg = colMeans(LZ1.samples[iter:(iter-499), ])
      lambda_L1 = max(2e-3, quantile(abs(LZ1_avg), probs = 0.8))
      L1 = (abs(LZ1_avg) > lambda_L1) * LZ1_avg
      LZ1 = LZ1_avg
    }
    # Note: Update in log scale
    if(iter > 2500 && iter %% 20 == 0){
      matL1 = matL1_prop = matrix(0, p, p)
      matL1[lower.tri(matL1)] = L1.samples[iter, ]
      Z1 = t(matU) %*% matD %*% (diag(p) - matL1) %*% Y[, 1:n1]
      log_accept_prob = - ll_Z(Z1, PHI, SIGMA, ID=1) - log(lambda_L1)
      llambda_L1_proposal = log(lambda_L1) + rnorm(1, 0, sigma_L) # sigma_L
      llambda_L1_proposal = log(lambda_L_U) * (llambda_L1_proposal > log(lambda_L_U)) + llambda_L1_proposal * (llambda_L1_proposal <= log(lambda_L_U))
      lambda_L1_proposal =  max(2e-3, exp(llambda_L1_proposal))
      proposal_L1 = (abs(LZ1.samples[iter, ])>lambda_L1_proposal) * LZ1.samples[iter, ]
      matL1_prop[lower.tri(matL1_prop)] = proposal_L1
      Z1 = t(matU) %*% matD %*% (diag(p) - matL1_prop) %*% Y[, 1:n1]
      log_accept_prob = log_accept_prob + ll_Z(Z1, PHI, SIGMA, ID=1) + llambda_L1_proposal
      
      if(is.na(log_accept_prob) || is.nan(log_accept_prob))
        log_accept_prob = -Inf
      
      if(log(runif(1)) < log_accept_prob){
        L1.samples[iter, ] = L1 = proposal_L1
        arcLL1 = arcLL1 + 1
        lambda_L1 = lambda_L1_proposal
        matL1 = matL1_prop
      }
      # cat("\n lambda_L", lambda_L)
      # Draw sigma_T
      ae_L1 = 0.01 + p * (p-1)/2
      be_L1 = 0.01 + sum(LZ1 ^ 2)/2
      
      sigma_T_L1 = sqrt(1/rgamma(1, ae_L1, be_L1))
    }
    # print(lambda_L1)
    
    # Draw lambda_L2
    if(iter == 2500){
      LZ2_avg = colMeans(LZ2.samples[iter:(iter-499), ])
      lambda_L2 = max(2e-3, quantile(abs(LZ2_avg), probs = 0.8))
      L2 = (abs(LZ2_avg) > lambda_L2) * LZ2_avg
      LZ2 = LZ2_avg
    }
    # Note: Update in log scale
    if(iter > 2500 && iter %% 20 == 0){
      matL2 = matL2_prop = matrix(0, p, p)
      matL2[lower.tri(matL2)] = L2.samples[iter, ]
      Z2 = t(matU) %*% matD %*% (diag(p) - matL2) %*% Y[, n1+1:n2]
      log_accept_prob = - ll_Z(Z2, PHI, SIGMA, ID=2) - log(lambda_L2)
      llambda_L2_proposal = log(lambda_L2) + rnorm(1, 0, sigma_L) # sigma_L
      llambda_L2_proposal = log(lambda_L_U) * (llambda_L2_proposal > log(lambda_L_U)) + llambda_L2_proposal * (llambda_L2_proposal <= log(lambda_L_U))
      lambda_L2_proposal =  max(2e-3, exp(llambda_L2_proposal))
      proposal_L2 = (abs(LZ2.samples[iter, ])>lambda_L2_proposal) * LZ2.samples[iter, ]
      matL2_prop[lower.tri(matL2_prop)] = proposal_L2
      Z2 = t(matU) %*% matD %*% (diag(p) - matL2_prop) %*% Y[, n1+1:n2]
      log_accept_prob = log_accept_prob + ll_Z(Z2, PHI, SIGMA, ID=2) + llambda_L2_proposal
      
      if(is.na(log_accept_prob) || is.nan(log_accept_prob))
        log_accept_prob = -Inf
      
      if(log(runif(1)) < log_accept_prob){
        L2.samples[iter, ] = L2 = proposal_L2
        arcLL2 = arcLL2 + 1
        lambda_L2 = lambda_L2_proposal
        matL2 = matL2_prop
      }
      # cat("\n lambda_L", lambda_L)
      # Draw sigma_T
      ae_L2 = 0.01 + p * (p-1)/2
      be_L2 = 0.01 + sum(LZ2 ^ 2)/2
      
      sigma_T_L2 = sqrt(1/rgamma(1, ae_L2, be_L2))
    }
    
    # Draw lambda_A
    if(iter == 2500){
      AZ = AZ_avg = colMeans(AZ.samples[iter:(iter-499), ])
      lambda_A = max(2e-4, quantile(abs(AZ_avg), probs = 0.8))
      A = (abs(AZ_avg) > lambda_A) * (abs(AZ_avg) - lambda_A) * sign(AZ_avg)
    }
    # Note: Update in log scale
    if(iter > 2500 && iter %% 30 == 0){
      matA = matA_prop = matrix(0, p, p)
      matA[lower.tri(matA)] = A.samples[iter, ]
      matA = matA - t(matA)
      matU = (diag(p) - matA) %*% solve(diag(p) + matA)
      Z = t(matU) %*% matD %*% cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                                      (diag(p) - matL2) %*% Y[, n1+1:n2])
      log_accept_prob = - ll_Z(Z, PHI, SIGMA, ID=NULL) - log(lambda_A)
      llambda_A_proposal = log(lambda_A) + rnorm(1, 0, sigma_A)
      llambda_A_proposal = log(lambda_A_U) * (llambda_A_proposal > log(lambda_A_U)) + llambda_A_proposal * (llambda_A_proposal <= log(lambda_A_U))
      lambda_A_proposal =  max(2e-4, exp(llambda_A_proposal))
      proposal_A = (abs(AZ.samples[iter, ])>lambda_A_proposal) * (abs(AZ.samples[iter, ]) - lambda_A_proposal) * sign(AZ.samples[iter, ])
      matA_prop[lower.tri(matA_prop)] = proposal_A
      matA_prop = matA_prop - t(matA_prop)
      matU_prop = (diag(p) - matA_prop) %*% solve(diag(p) + matA_prop)
      Z = t(matU_prop) %*% matD %*% cbind( (diag(p) - matL1) %*% Y[, 1:n1], 
                                           (diag(p) - matL2) %*% Y[, n1+1:n2])
      log_accept_prob = log_accept_prob + ll_Z(Z, PHI, SIGMA, ID=NULL) + llambda_A_proposal
      
      if(is.na(log_accept_prob) || is.nan(log_accept_prob))
        log_accept_prob = -Inf
      
      if(log(runif(1)) < log_accept_prob){
        A.samples[iter, ] = A = proposal_A
        arcLU = arcLU + 1
        lambda_A = lambda_A_proposal
        matU = matU_prop
      }
      # cat("\n lambda_A", lambda_A)
      # Draw sigma_T
      ae_A = 0.1 + p * (p-1)/2
      be_A = 0.1 + sum(AZ ^ 2)/2
      
      sigma_T_A = sqrt(1/rgamma(1, ae_A, be_A))
    }
    
    # # Adapt the Covariance matrix
    # # For A
    # if(adaptA){
    #   # if(iter > 5000 && iter <= 7500){
    #   #    Cov_prop_A = cor(as.matrix(A.samples[(iter-499):iter, ])) + 0.01 * diag((p*(p-1)/2))
    #   # }
    #   if(iter > 3500 && iter %% 100 == 0){
    #     Cov_prop_A = cor(as.matrix(A.samples[(iter-499):iter, ])) + 0.01 * diag((p*(p-1)/2))
    #   }
    # }
    
    # Check acceptance ratio
    if(iter%%100 == 0 && iter > 1500){
      a_rate_A = arcA / ( iter - 1500 )
      # if(a_rate_A < 0.05)
      #   sdA = sdA * 0.5
      if(a_rate_A < 0.15)
        sdA = sdA / 5
      else if(a_rate_A > 0.40)
        sdA = sdA * 5
    }
    # For AR processes
    if(adaptAR){
      for (i in 1:p){
        # if(iter > 5000 && iter <= 7500){
        #    Cov_prop_AR[,,i] = cor(as.matrix(PHI.samples[(iter-499):iter, i, ])) + 0.01 * diag(K)
        # }
        if(iter > 7500 && iter %% 100 == 0){
          Cov_prop_AR[,,i] = cor(as.matrix(PHI.samples[(iter-499):iter, i, ])) + 0.01 * diag(K)
        }
      }
    }
    
    # Check acceptance ratio
    if(iter%%100 == 0){
      for(i in 1:p){
        a_rate = arcAR[i] / iter
        if(a_rate < 0.15)
          sdAR[i] = sdAR[i] / 10
        else if(a_rate > 0.40)
          sdAR[i] = sdAR[i] * 10
      }
    }
    
    # Adapt step size for LMC
    # For L1
    if(iter%%100 == 0 && iter > 1500){
      # a_rate_L = arcL / iter
      a_rate_L1 = arcL1 / (iter - 1500)
      if(a_rate_L1 < 0.45)
        eps_L = eps_L / 10
      else if(a_rate_L1 > 0.70)
        eps_L = eps_L * 10
    }
    # For L2
    if(iter%%100 == 0 && iter > 1500){
      # a_rate_L = arcL / iter
      a_rate_L2 = arcL2 / (iter - 1500)
      if(a_rate_L2 < 0.45)
        eps_L = eps_L / 10
      else if(a_rate_L2 > 0.70)
        eps_L = eps_L * 10
    }
    # For D
    if(iter%%100 == 0 && iter > 1500){
      # a_rate_D = arcD / iter
      a_rate_D = arcD / (iter - 1500)
      if(a_rate_D < 0.45)
        eps_D = eps_D / 10
      else if(a_rate_D > 0.70)
        eps_D = eps_D * 10
    }
    
    # Print acceptance ratios and Monitor RMSE
    if(verbose){
      if(iter%%50 == 0 && iter > 1500){
        # RMSE = sqrt(mean(c(Omega_pred - Omega0)^2))
        cat("\n Iterations: ", iter)
      }
      if(iter%%200 == 0 && iter > 1500){
        # cat("> ")
        # Omega_hat = (diag(p) - matL) %*% matD^2 %*% t(diag(p) - matL)
        # RMSE = sqrt(mean(c(Omega_pred - Omega0)^2))
        # print(SIGMA)
        # print(lambda_L)
        cat("\n\n Iterations: ", iter, "| Acceptance Rates: L:", a_rate_L1, a_rate_L2, " D:", a_rate_D, " A:", a_rate_A, "\n")
        cat("\n Average acceptance rate for AR processes: ", mean(arcAR)/iter)
        cat("\n lambda values for L and A:", lambda_L1, lambda_L2, lambda_A, "\n")
      }
    }
  }
  return(list(Omega1.samples = Omega1.samples, 
              Omega2.samples = Omega2.samples, 
              A.samples = A.samples, 
              L1.samples = L1.samples, 
              L2.samples = L2.samples, 
              D.samples = D.samples, 
              PHI.samples = PHI.samples, 
              SIGMA.samples = SIGMA.samples, 
              Y_onestep_exp = Y_onestep_exp,
              Y_onestep_raw = Y_onestep_raw
              # Omega.GGM = GGM$K_hat
              # RMSE_GGM = RMSE_GGM,
              # RMSE_GCGM = RMSE_GCGM
  ))
} # End of OUTAR


results = rep(0, 2)
Y = t(rbind(oecd1, oecd2))
cat("\n", dim(Y))
tm1 = tic()
model = OUTAR(Y = 1000*Y[, 1:78], K = 2, n1 = 38, n2 = 40, n.iter = 10000, verbose = TRUE, adaptA = FALSE, adaptAR = FALSE)
tm2 = toc()

omegas1 = model$Omega1.samples[5001:10000,,] # before
omegas2 = model$Omega2.samples[5001:10000,,] # after

# thresholded omegas
threshold = 0.0015
omegas1[abs(omegas1) < threshold] = 0
omegas2[abs(omegas2) < threshold] = 0

omegas_diff = omegas2 - omegas1

Omega1.pred = apply(model$Omega1.samples[5001:10000,,], c(2, 3), mean) # before
Omega2.pred = apply(model$Omega2.samples[5001:10000,,], c(2, 3), mean) # after

Omega_diff.pred = Omega2.pred - Omega1.pred  # posterior mean

lower = apply(omegas_diff, c(2,3), quantile, probs=0.05)
upper = apply(omegas_diff, c(2,3), quantile, probs=0.975)

adjmat = apply((lower > 0 | upper < 0), c(1,2), as.integer)  # flag 1 if 95% credible interval doesn't contain 0.
W.pred = graph_from_adjacency_matrix(adjmat, mode = "undirected", diag=FALSE)
E(W.pred)$width = 1.5
E(W.pred)$length = 4
V(W.pred)$name = countries$`Reference.area
# layout <- layout_with_fr(W.pred)     # Fruchterman-Reingold (force-directed)
# layout <- layout_with_kk(W.pred)  # Kamada-Kawai (spring layout)
# layout <- layout_in_circle(W.pred)  # Circular layout
# layout <- layout_with_drl(W.pred)   # Large graphs (DRL layout)
# plot(W.pred,
#      layout = layout,
#      vertex.size = 15,          # Increase node size
#      vertex.label.cex = 0.7,    # Decrease font size of labels
#      vertex.label.color = "black",  # Optional: set label color
#      vertex.color = "skyblue")

# Determine the sign of significant change
edge_signs = matrix(0, nrow = nrow(adjmat), ncol = ncol(adjmat))
edge_signs[lower > 0] = 1   # Positive change
edge_signs[upper < 0] = -1  # Negative change

# Extract edge list to match color with edge
edge_list = as_edgelist(W.pred, names = FALSE)

# Initialize color vector
edge_colors = character(ecount(W.pred))

# Assign color based on sign of change
for (i in seq_len(nrow(edge_list))) {
  v1 <- edge_list[i, 1]
  v2 <- edge_list[i, 2]
  sign_val <- edge_signs[v1, v2]
  
  edge_colors[i] <- if (sign_val == 1) {
    "blue"     # positive effect
  } else if (sign_val == -1) {
    "red"      # negative effect
  } else {
    "black"    # fallback (should not happen with filtered adjmat)
  }
}

# Assign to graph
E(W.pred)$color <- edge_colors


# Total number of nodes
n <- vcount(W.pred)

# Example: split nodes into inner (1/3) and outer (2/3) circles
n_inner <- ceiling(n / 3)
n_outer <- n - n_inner

# Radii for inner and outer circles
r_inner <- 1
r_outer <- 2

# Create empty layout matrix
layout <- matrix(NA, nrow = n, ncol = 2)

# Inner circle coordinates
angles_inner <- seq(0, 2 * pi, length.out = n_inner + 1)[- (n_inner + 1)]
layout[1:n_inner, ] <- cbind(
  r_inner * cos(angles_inner),
  r_inner * sin(angles_inner)
)

# Outer circle coordinates
angles_outer <- seq(0, 2 * pi, length.out = n_outer + 1)[- (n_outer + 1)]
layout[(n_inner + 1):n, ] <- cbind(
  r_outer * cos(angles_outer),
  r_outer * sin(angles_outer)
)

# Plot with custom layout
plot(
  W.pred,
  layout = layout,
  vertex.size = 17,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",  # Optional: set label color
  vertex.color = "oldlace", # ivory, ghostwhite, whiteSmoke, oldlace
  vertex.frame.width = 2
)

# Add legend
legend(
  "topright",
  legend = c("Increase after Recession", "Decrease after Recession"),
  col = c("blue", "red"),
  lty = 1,
  lwd = 2,
  title = "Shifts in Conditional Economic \nDependencies After the 2009 Recession",
  box.lwd = 0,
  inset = 0.001
)

Y_hat_raw = colMeans(model$Y_onestep_raw[5001:10000, ]) # one step ahead predictions
Y_hat_exp = colMeans(model$Y_onestep_exp[5001:10000, ])
print(Y_hat_raw)
print(Y_hat_exp)
print(as.numeric(tm2$toc - tm2$tic))
results[1] = sqrt(mean((10000*Y[, 79] - Y_hat_raw)^2)) / 10000
results[2] = sqrt(mean((10000*Y[, 79] - Y_hat_exp)^2)) / 10000
# 
# ci_raw = apply(model$Y_onestep_raw[5001:10000, ], 2, quantile, probs=c(0.05, 0.975))
# ci_exp = apply(model$Y_onestep_exp[5001:10000, ], 2, quantile, probs=c(0.05, 0.975))
# 
# D = 1 / sqrt(outer(diag(Omega.pred), diag(Omega.pred)))
# O.pred = Omega.pred * D
# adj.pred = abs(O.pred) > 0.03 # quantile(c(abs(O.pred)), 0.9)
# W.pred = graph_from_adjacency_matrix(adj.pred, mode = "undirected", diag=FALSE)
# E(W.pred)$width = 2.5
# plot(W.pred)
# cat("\n")
# print(round(results, 4))
# # write.csv(x=round(results, 4), file=paste0("/share/statistics/sghosh27/BayesTS/Final/realdata/ospvals/result_l_", args1, ".csv"))
# # saveRDS(Omega.pred, file = paste0("/share/statistics/sghosh27/BayesTS/Final/realdata/omegavals/Omega_l_", args1, ".rds"))
# 
# # from GGM
# Omega.pred1 = model$Omega.GGM
# D1 = 1 / sqrt(outer(diag(Omega.pred1), diag(Omega.pred1)))
# O.pred1 = Omega.pred1 * D1
# adj.pred1 = abs(O.pred1) > 0.1
# W.pred1 = graph_from_adjacency_matrix(adj.pred1, mode = "undirected", diag=FALSE)
# E(W.pred1)$width = 2.5
# plot(W.pred1)
# 
# 
# omegas = model$Omega.samples[5001:10000,,]
# 
# # Create a list to store the graph objects
# graphs <- list()
# 
# # Loop over each precision matrix and create the corresponding graph
# for (i in 1:5000) {
#   D = 1 / sqrt(outer(diag(omegas[i,,]), diag(omegas[i,,])))
#   o.pred = omegas[i,,] * D
#   adj_matrix <- abs(o.pred) > 0.1  # Adjacency matrix from omega # 9e-2
#   graphs[[i]] <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
#   E(graphs[[i]])$width = 2.5
# }
# 
# prob = rep(NA, 5000)
# for (i in 1:5000) {
#   D = 1 / sqrt(outer(diag(omegas[i,,]), diag(omegas[i,,])))
#   o.pred = omegas[i,,] * D
#   # adj_matrix <- abs(o.pred) > 0.1  # Adjacency matrix from omega # 9e-2
#   prob[i] = mean(abs(o.pred) > 0.1)
#   # graphs[[i]] <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
#   # E(graphs[[i]])$width = 2.5
# }
# 
# # Combine all the graphs into one unified graph (overlay)
# g_combined <- graphs[[1]]  # Start with the first graph
# 
# for (i in 2:5000) {
#   g_combined <- union(g_combined, graphs[[i]])  # Overlay all graphs using union
# }
# 
# # Plot the combined graph with distinct edge colors and line types for each graph
# plot(g_combined, vertex.label = 1:vcount(g_combined),
#      main = "Overlay of Graphs from Multiple Precision Matrices")
