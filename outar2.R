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

# orders 2, 5, 8, 10
# p = 3 * 10, 3 * 20, 3 * 30
# T = 50, 100 and 150 for 2 and 5 orders and
# T = 150, 200 and 250 for orders 8 and 10

dict = list(list(50, 2, 0.9, 3, 10, 2), 
            list(100, 2, 0.9, 3, 10, 2),
            list(150, 2, 0.9, 3, 10, 2),
            list(50, 2, 0.9, 3, 20, 2), 
            list(100, 2, 0.9, 3, 20, 2),
            list(150, 2, 0.9, 3, 20, 2), 
            list(50, 2, 0.9, 3, 30, 2), 
            list(100, 2, 0.9, 3, 30, 2),
            list(150, 2, 0.9, 3, 30, 2), 
            list(50, 5, 0.9, 3, 10, 5), 
            list(100, 5, 0.9, 3, 10, 5),
            list(150, 5, 0.9, 3, 10, 5),
            list(50, 5, 0.9, 3, 20, 5), 
            list(100, 5, 0.9, 3, 20, 5),
            list(150, 5, 0.9, 3, 20, 5), 
            list(50, 5, 0.9, 3, 30, 5), 
            list(100, 5, 0.9, 3, 30, 5),
            list(150, 5, 0.9, 3, 30, 5),
            list(150, 8, 0.9, 3, 10, 8), 
            list(200, 8, 0.9, 3, 10, 8),
            list(250, 8, 0.9, 3, 10, 8),
            list(150, 8, 0.9, 3, 20, 8), 
            list(200, 8, 0.9, 3, 20, 8),
            list(250, 8, 0.9, 3, 20, 8), 
            list(150, 8, 0.9, 3, 30, 8), 
            list(200, 8, 0.9, 3, 30, 8),
            list(250, 8, 0.9, 3, 30, 8),
            list(150, 10, 0.9, 3, 10, 10), 
            list(200, 10, 0.9, 3, 10, 10),
            list(250, 10, 0.9, 3, 10, 10),
            list(150, 10, 0.9, 3, 20, 10), 
            list(200, 10, 0.9, 3, 20, 10),
            list(250, 10, 0.9, 3, 20, 10), 
            list(150, 10, 0.9, 3, 30, 10), 
            list(200, 10, 0.9, 3, 30, 10),
            list(250, 10, 0.9, 3, 30, 10)
            # list(40, c(rep(2, 50), rep(5, 10)), c(rep(0.9, 50), rep(0.9, 10)), 3, 20, 2),
            # list(40, c(rep(2, 60), rep(5, 30)), c(rep(0.9, 60), rep(0.9, 30)), 3, 30, 2),
            # list(40, c(rep(2, 90), rep(5, 30)), c(rep(0.9, 90), rep(0.9, 30)), 4, 30, 2),
            # list(40, c(rep(2, 10), rep(5, 50)), c(rep(0.9, 10), rep(0.9, 50)), 3, 20, 3),
            # list(40, c(rep(2, 30), rep(5, 60)), c(rep(0.9, 30), rep(0.9, 60)), 3, 30, 3),
            # list(40, c(rep(2, 50), rep(5, 70)), c(rep(0.9, 50), rep(0.9, 70)), 4, 30, 3),
            # list(40, c(rep(2, 10), rep(3, 30), rep(5, 20)), c(rep(0.9, 10), rep(0.9, 30), rep(0.9, 20)), 3, 20, 3),
            # list(40, c(rep(2, 20), rep(3, 50), rep(5, 20)), c(rep(0.9, 20), rep(0.9, 50), rep(0.9, 20)), 3, 30, 3),
            # list(40, c(rep(2, 60), rep(3, 20), rep(5, 40)), c(rep(0.9, 60), rep(0.9, 20), rep(0.9, 40)), 4, 30, 3),
            # list(40, c(rep(2, 10), rep(3, 10), rep(5, 40)), c(rep(0.9, 10), rep(0.9, 10), rep(0.9, 40)), 3, 20, 5),
            # list(40, c(rep(2, 20), rep(3, 30), rep(5, 40)), c(rep(0.9, 20), rep(0.9, 30), rep(0.9, 40)), 3, 30, 5),
            # list(40, c(rep(2, 40), rep(3, 20), rep(5, 60)), c(rep(0.9, 40), rep(0.9, 20), rep(0.9, 60)), 4, 30, 5),
            # list(40, c(rep(2, 40), rep(3, 10), rep(5, 10)), c(rep(0.9, 40), rep(0.9, 10), rep(0.9, 10)), 3, 20, 2),
            # list(40, c(rep(2, 60), rep(3, 10), rep(5, 20)), c(rep(0.9, 60), rep(0.9, 10), rep(0.9, 20)), 3, 30, 3),
            # list(40, c(rep(2, 70), rep(3, 20), rep(5, 30)), c(rep(0.9, 70), rep(0.9, 20), rep(0.9, 30)), 4, 30, 5),
            # list(200, c(rep(8, 40), rep(10, 20)), c(rep(0.9, 40), rep(0.9, 20)), 3, 20, 8),
            # list(200, c(rep(8, 60), rep(10, 30)), c(rep(0.9, 60), rep(0.9, 30)), 3, 30, 8),
            # list(200, c(rep(8, 80), rep(10, 40)), c(rep(0.9, 80), rep(0.9, 40)), 4, 30, 8),
            # list(200, c(rep(2, 20), rep(7, 20), rep(10, 20)), c(rep(0.9, 20), rep(0.9, 20), rep(0.9, 20)), 3, 20, 8),
            # list(200, c(rep(2, 20), rep(7, 30), rep(10, 40)), c(rep(0.9, 20), rep(0.9, 30), rep(0.9, 40)), 3, 30, 8),
            # list(200, c(rep(2, 30), rep(7, 30), rep(10, 60)), c(rep(0.9, 30), rep(0.9, 30), rep(0.9, 60)), 4, 30, 8)
            )

args = commandArgs(trailingOnly = TRUE)
argslen = length(args)
if(argslen > 1) stop('Error: Too Many Arguments')
if(argslen < 1) stop('Error: Takes One Argument Only')
args1 = dict[[as.numeric(args)]]

# args1 = dict[[30]]   # check with 28-30

cat("\n Process parameters: ", unlist(args1), "\n")

# simulate multivariate time-series data
generate_data = function(Tm, order, corr = 0.9, n_graphs=3, n_nodes=20, prob=0.1){
  # Tm: the number of time-points in each time series
  # order: orders of the underlying AR processes. If a scalar, all AR processes have the same order
  # prob (q): Sparsity level in each small-world network
  # n_graphs: Number of small world networks
  # n_nodes: Number of nodes in each small world network
  # returns: (1) Y: A p-multivariate time series with Tm time points each
  #          (2) Omega0: The conditional dependence structure in Y
  
  G = W = list()
  if(length(n_nodes)==1)
    n_nodes = rep(n_nodes, n_graphs)
  p = sum(n_nodes)
  if(length(order)==1)
    order = rep(order, p)
  if(length(corr)==1)
    corr = rep(corr, p)
  index <- as.matrix(combinat::combn(1:p, 2))
  for(graph in 1:n_graphs){
    G[[graph]] = sample_smallworld(1, size=n_nodes[graph], nei=10, p=prob)
    W[[graph]] = as.matrix(as_adjacency_matrix(G[[graph]], sparse = F))
  }
  W = bdiag(W)
  W[1:(p/2), (p/2+1):p] = rbinom((p/2)^2, 1, 0.02)
  W[(p/2+1):p, 1:(p/2)] = t(W[1:(p/2), (p/2+1):p])
  OmegaA = as.matrix(W)
  # image(W)
  W.g = graph_from_adjacency_matrix(W, mode = "undirected")
  E(W.g)$width = 2.5
  # plot(W.g)
  
  mean(OmegaA[upper.tri(OmegaA)]==1)
  
  Omega0 = rgwish( n = 1, adj = OmegaA, b = 10, D = diag(p), threshold = 1e-8 )
  Omega0 = OmegaA * Omega0 + diag(diag(Omega0))
  Omega0[which(abs(Omega0)<1)] = 0
  
  vec = abs(Omega0[t(index)])
  # cat("\n", mean(vec!=0)*100, "% \n")
  
  dis   = dist(1:Tm)
  phi.true = matrix(NA, p, max(order))
  sigma.true = rep(NA, p)
  Z = matrix(0, p, Tm)
  adj.v = matrix(0, Tm, Tm)
  adj.v = (abs(row(adj.v)-col(adj.v))<10)
  diag(adj.v) = rep(0, Tm)
  
  # for(i in 1:p){
  #   pacfs_sign = runif(order[i], -1, 1)
  #   pacfs = runif(order[i], corr[i], 1) * sign(pacfs_sign)
  #   # pacfs = runif(order[i], -1, 1)
  #   phi.true[i, 1:order[i]] = INLA::inla.ar.pacf2phi(pacfs)
  #   acfs = INLA::inla.ar.pacf2acf(pacfs, lag.max = order[i])[-1]
  #   # print(sum(acfs*phi.true[i, ]))
  #   sigma.true[i] = sqrt(abs(1 - sum(acfs * phi.true[i, 1:order[i]])))
  #   Z[i, ] = astsa::sarima.sim(ar = phi.true[i, 1:order[i]], n = Tm, sd = sigma.true[i])
  #   # Z[i, ] = arima.sim(model = list(ar = phi.true[i, 1:order[i]]), n = Tm, sd = sigma.true[i])
  # }
  
  for (i in 1:p) {
    success = FALSE
    while (!success) {
      tryCatch({
        pacfs_sign = runif(order[i], -1, 1)
        pacfs = runif(order[i], corr[i], 1) * sign(pacfs_sign)
        phi.true[i, 1:order[i]] = INLA::inla.ar.pacf2phi(pacfs)
        acfs = INLA::inla.ar.pacf2acf(pacfs, lag.max = order[i])[-1]
        sigma.true[i] = sqrt(abs(1 - sum(acfs * phi.true[i, 1:order[i]])))
        
        # Attempt to generate the AR process simulation
        Z[i, ] = astsa::sarima.sim(ar = phi.true[i, 1:order[i]], n = Tm, sd = sigma.true[i])
        success = TRUE  # If no error, set success to TRUE
      }, error = function(e) {
        # message("Error in AR simulation: ", e$message, ". Retrying...")
        success = FALSE  # Continue the loop to retry
      })
      if(diff(range(Z[i,]))>5){
        success = FALSE
      }
    }
  }
  
  # print(phi.true)
  # print(sigma.true)
  # A.values = rnorm(p*(p-1)/2) 
  # skvec = rep(0, p*(p-1)/2)
  A = matrix(0, p, p)
  # A[upper.tri(A)] = A.values
  # A = A - t(A)
  U = (diag(p) + A) %*% solve(diag(p) - A)
  svpd = svd(Omega0)
  sqrtroot = svpd$u %*% diag(1/sqrt(svpd$d)) %*% t(svpd$v)
  Y = sqrtroot %*% U %*% Z
  
  return(list(Y=Y, Omega0 = Omega0))
}

# trial run
# data = generate_data(Tm=200, order=2)
# data = generate_data(Tm=50, order=c(rep(2, 50), rep(5, 10)))


OUTAR = function(Y, K, n.iter = 1, verbose = T, adaptA = T, adaptAR = T){
  # Y: p-variate time series data
  # K: K of the AR models fit to the latent components
  # verbose: If TRUE, model displays the number of iterations, acceptance probability and RMSE of the estimator of Omega
  # adapt: If TRUE, the adaptive covariance function of Haario et al. (2001) is used to evaluate the covariance matrix of the proposal.
  # n.iter: Number of Monte Carlo iterations
  
  cat("\n")
  # p: dimension of the time series data
  p = nrow(Y)
  # N: number of observations in each time series
  N = ncol(Y)
  # Store the parameters. Note that the matrices are stored using their vectorized representations.
  A.samples = array(NA, dim=c(n.iter, p*(p-1)/2))
  AZ.samples = array(NA, dim=c(n.iter, p*(p-1)/2))
  D.samples = array(NA, dim=c(n.iter, p))
  L.samples = array(NA, dim=c(n.iter, p*(p-1)/2))
  LZ.samples = array(NA, dim=c(n.iter, p*(p-1)/2))
  PHI.samples = array(NA, dim=c(n.iter, p, K))
  SIGMA.samples = array(NA, dim=c(n.iter, p))
  Omega.samples = array(NA, dim=c(n.iter, p, p))
  PACF.samples = array(NA, dim=c(n.iter, p, K))
  
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
  GGM = bdgraph( data = t(Y), method="ggm", iter = 10000, df.prior = 10 )
  RMSE_GGM = sqrt(mean(c(Omega0 - GGM$K_hat)^2))
  cat("\n RMSE by GGM: ", RMSE_GGM, "\n")
  GCGM = bdgraph( data = t(Y), method="gcgm", iter = 10000 )
  RMSE_GCGM = sqrt(mean(c(Omega0 - GCGM$K_hat)^2))
  cat("\n RMSE by GCGM: ", RMSE_GCGM, "\n")
  
  # Initialize L, D and A
  Omega = GGM$K_hat
  cat("\n Running process with order = ", K, "\n")
  L.chol = t(chol(Omega))
  # D.samples[1, ] = diag(IL)
  IL = L.chol %*% diag( 1 / diag(L.chol) )
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
  lambda_A = lambda_L = 0
  # eps_L = 1e-3 / (p*(p-1)/2); eps_D = 1e-3 / p # initial step size
  eps_L = 1e-6; eps_D = 1e-6;
  sigma_D = 10
  sigma_A = sigma_L = 1e-2
  
  # compute Z from U, D L and Y
  # mat<> stores <> in matrix format
  matA = matrix(0, p, p)
  matA[lower.tri(matA)] = A = AZ = rep(0, p*(p-1)/2)
  matA = matA - t(matA)
  matU = (diag(p) - matA) %*% solve(diag(p) + matA)
  matD = D = diag(L.chol)
  matD = diag(matD)
  matL = matrix(0, p, p)
  LZ = - IL[lower.tri(IL)]
  L = LZ * (abs(LZ) > lambda_L)
  matL[lower.tri(matL)] = L
  Z = t(matU) %*% matD %*% (diag(p) - matL) %*% Y
  
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
  
  ae_L = 0.1 + p * (p-1)/2
  be_L = 0.1 + sum(LZ ^ 2)/2
  
  sigma_T_L = 1/rgamma(1, ae_L, be_L)
  
  # Joint log-likelihood of a single AR(K) processes ----- (1)
  # Note: Y has L, D and U
  ll_Zi = function(tms, phi, sigma){
    tsvals = embed(ts(tms), dimension = K + 1)
    return( - 0.5 * ( 2 * N * log(sigma) + sum((tsvals[, 1] - tsvals[, -1] %*% phi)^2) /  sigma^2 ) )
  }
  # print(ll_Zi(Z[1,], PHI[1,], SIGMA[1]))
  # Wrapper for the joint log-likelihood of the p AR(K) processes
  ll_Z = function(data, phi.values, sigma.values){
    # numCores = detectCores() - 1
    # cl = makeCluster(numCores)
    # registerDoParallel(cl)
    # result = rep(NA, p)
    # results = foreach(i = 1:p, .combine = c, .export = c("N", "K", "ll_Zi")) %dopar% {
    # for(i in 1:p){
    #   result[i] = ll_Zi(tms=data[i, ], phi=phi.values[i,], sigma=sigma.values[i])
    # }
    result = mclapply(1:p, function(i) ll_Zi(tms=data[i, ], phi=phi.values[i,], sigma=sigma.values[i]), mc.cores = 1L)
    # stopCluster(cl)
    # print(result)
    return( sum(unlist(result)) )
  }
  
  # print(ll_Z(Z, PHI, SIGMA))
  # Gradient of the negative log-likelihood of Z wrt Z
  gradient_ll_Zi = function(tms, phi, sigma) {
    
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
  
  
  gradient_ll_Z = function(Z, phi.values, sigma.values){
    # Number of cores to use
    # numCores = detectCores() - 1
    # cl = makeCluster(numCores)
    # registerDoParallel(cl)
    grad_Z = NULL
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
    grad_Z = mclapply(1:p, function(i) gradient_ll_Zi(tms=Z[i, ], phi=PHI[i, ], sigma=SIGMA[i]), mc.cores = 1L)
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
  gradient_L = function(L.values, LZ.values, lambda, sigma_T, h0=1e-8){
    # Return the gradient as a flattened vector
    # Contribution from (1)
    grad_ll = - kronecker(Y, matD %*% matU)
    # print(dim(grad_nll))
    grad_ll = matrix(grad_ll %*% c(gradient_ll_Z(Z, PHI, SIGMA)), p, p)
    grad_ll = grad_ll[lower.tri(grad_ll)]
    
    # Contribution from threshholding function
    grad_thrshld = 0.5 * (1 + 2 * (pi)^(-1) * atan((LZ.values^2 - lambda^2)/h0)) + 2 * LZ.values^2 / ( pi*h0 * ( 1 + (LZ.values^2 - lambda^2)^2 / h0^2 ) )
    
    return(grad_ll*grad_thrshld - LZ.values/sigma_T^2)
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
  gradient_logD = function(D.values, mu){
    # Contribution from (1)
    grad_ll = kronecker((diag(p) - matL) %*% Y, matU)
    # print(dim(grad_nll))
    grad_ll = matrix(grad_ll %*% c(gradient_ll_Z(Z, PHI, SIGMA)), p, p)
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
  update_L = function(init_L, init_LZ, step_size, lambda, sigma_T, ar_counter){
    
    # print(init_LZ)
    # vecL = init_L
    Z = t(matU) %*% matD %*% (diag(p) - matL) %*% Y
    log_den = ll_Z(Z, PHI, SIGMA) - 0.5 * sum(init_LZ^2) / sigma_T^2
    gradient = gradient_L(L.values=init_L, LZ.values=init_LZ, lambda=lambda, sigma_T=sigma_T)
    proposal_mean = init_LZ + 0.5 * step_size * gradient
    proposal_LZ = proposal_mean + sqrt(step_size) * rnorm(p*(p-1)/2)
    # print(proposal_LZ)
    proposal_L = proposal_LZ * (abs(proposal_LZ) > lambda) 
    log_den = log_den - 0.5 * sum((proposal_LZ - proposal_mean)^2) / step_size
    
    matL_prop = matL
    matL_prop[lower.tri(matL_prop)] = proposal_L
    Z = t(matU) %*% matD %*% (diag(p) - matL_prop) %*% Y
    log_num = ll_Z(Z, PHI, SIGMA) - 0.5 * sum(proposal_LZ^2) / sigma_T^2
    proposal_gradient = gradient_L(L.values=proposal_L, LZ.values=proposal_LZ, lambda=lambda, sigma_T=sigma_T)
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
  
  update_D = function(init_D, xi, step_size, ar_counter){
    
    init_logD = log(init_D)
    Z = t(matU) %*% diag(init_D) %*% (diag(p) - matL) %*% Y
    log_den = ll_Z(Z, PHI, SIGMA) + ll_logD(D.values=init_D, mu=xi)
    gradient = gradient_logD(D.values=init_D, mu=xi)
    proposal_mean = init_logD + 0.5 * step_size * gradient
    proposal_logD = proposal_mean + sqrt(step_size) * rnorm(p)
    log_den = log_den - 0.5 * sum((proposal_logD - proposal_mean)^2) / step_size
    # print(log_den)
    proposal_D = exp(proposal_logD)
    Z = t(matU) %*% diag(proposal_D) %*% (diag(p) - matL) %*% Y
    log_num = ll_Z(Z, PHI, SIGMA) + ll_logD(D.values=proposal_D, mu=xi)
    # print(log_num)
    proposal_gradient = gradient_logD(D.values=proposal_D, mu=xi)
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
  update_A = function(init_A, init_AZ, lambda, sigma_T, proposal_sigma, ar_counter){
    
    # A0 = matrix(0, p, p)
    # A0[lower.tri(A0)] = init_A
    # init_U = (diag(p) - A0) %*% solve(diag(p) + A0)
    Z = t(matU) %*% matD %*% (diag(p) - matL) %*% Y
    log_posterior = ll_Z(Z, PHI, SIGMA) - sum(init_AZ^2) / (2 * sigma_T^2)
    proposal_AZ = init_AZ + MASS::mvrnorm(1, mu=rep(0, p*(p-1)/2), Sigma = proposal_sigma)
    proposal_A = (abs(proposal_AZ) > lambda) * (abs(proposal_AZ) - lambda) * sign(proposal_AZ)
    # cat("proposal_A", proposal_A[1:10], "\n")
    A_prop = matrix(0, p, p)
    A_prop[lower.tri(A_prop)] = proposal_A
    matU_prop = (diag(p) - A_prop) %*% solve(diag(p) + A_prop)
    Z = t(matU_prop) %*% matD %*% (diag(p) - matL) %*% Y
    log_proposed_posterior = ll_Z(Z, PHI, SIGMA) - sum(proposal_AZ^2) / (2 * sigma_T^2)
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
    
    log_posterior = ll_Zi(tms=data, phi=current_phi, sigma=current_sigma) + sum(dlogis(qlogis(((current_pacf + 1) / 2)), location = 0, scale = 1, log = TRUE))
    
    logit_proposal = MASS::mvrnorm(1, mu = qlogis(((current_pacf + 1) / 2)), Sigma = proposal_sigma)   # V = logit((U + 1) / 2) ~ Logistic(0, 1)
    proposed_pacf = plogis(logit_proposal) * 2 - 1  # inverse_logit(V) * 2 - 1
    # print(proposed_pacf)
    proposed_phi = INLA::inla.ar.pacf2phi(proposed_pacf)
    proposed_acf = INLA::inla.ar.pacf2acf(proposed_pacf, lag.max = K)[-1]
    if(length(proposed_acf) == 0)
      log_proposed_posterior = NaN
    else{
        # print(current_acf)
        stopifnot(length(proposed_phi) == length(proposed_acf))
        proposed_sigma = sqrt(1 - sum(proposed_phi * proposed_acf))
        log_proposed_posterior = ll_Zi(tms=data, phi=proposed_phi, sigma=proposed_sigma) + sum(dlogis(logit_proposal, location = 0, scale = 1, log = TRUE))
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
  a_rate = a_rate_L = a_rate_A = a_rate_D = 0
  sdA = 1e-10; sdAR = rep(1e-8, p)
  arcA = arcL = arcD = arcLU = arcLL = 0; arcAR = rep(0, p)
  Cov_prop_AR = array(rep(diag(K), p), dim=c(K, K, p))
  Cov_prop_A = diag(p*(p-1)/2)
  
  # MCMC sampling. 10000 samples with 2000 burnin
  for(iter in 1:n.iter){
    
    # Draw sample for A (U)
    if(iter > 1500){
      result = update_A(init_A=A, init_AZ=AZ, lambda=lambda_A, sigma_T = sigma_T_A, proposal_sigma=sdA*Cov_prop_A, ar_counter=arcA)
      A.samples[iter,] = A = result$A
      AZ.samples[iter,] = AZ = result$AZ
      matU = result$matU
      arcA = result$arc
      # cat("\n check 1")  # ************** Checkpoint
      # cat("\n", A[30:50])
    }
    
    # Draw sample for L
    if(iter > 1500){
      result = update_L(init_L=L, init_LZ=LZ, step_size=eps_L, lambda=lambda_L, sigma_T = sigma_T_L, ar_counter=arcL)
      L.samples[iter,] = L = result$L
      LZ.samples[iter,] = LZ = result$LZ
      matL = result$matL
      # L = matL[lower.tri(matL)]
      arcL = result$arc
      # cat("\n check 2")  # ************** Checkpoint
      # cat("\n", sum(L>10), "\t", eps_L)
    }
    
    # Draw sample for D
    if(iter > 1500){
      xi = rnorm(1, 0, sigma_D) # mean of the inverse gaussian prior
      result = update_D(init_D=D, step_size=eps_D, xi=xi, ar_counter=arcD)
      D.samples[iter,] = D = result$D
      matD = result$matD
      # D = diag(matD)
      arcD = result$arc
      # cat("\n check 3")  # ************** Checkpoint
      # cat("\n", sum(D > 15), "\t", eps_D)
    }
    
    # print(matD[1:5, 1:5])
    
    if(iter >= 1500)
      Omega.samples[iter, , ] = Omega_pred = (diag(p) - matL) %*% (matD ^ 2) %*% t(diag(p) - matL)
    
    # print(mean((Omega_pred-Omega0)^2))
    
    # Draw samples for AR process parameters
    # numCores = detectCores() - 1
    # cl = makeCluster(numCores)
    # registerDoParallel(cl)
    # registerDoRNG()
    Z = t(matU) %*% matD %*% (diag(p) - matL) %*% Y
    result = NULL
    # Use foreach with doRNG to ensure reproducibility
    # result = foreach(i = 1:p, .packages = 'doRNG') %dopar% {
    # for(i in 1:p){
    #   result = rbind(result, update_AR(data=Z[i,], current_params=c(PHI[i,], PACF[i, ], SIGMA[i]), proposal_sigma=sdAR[i]*Cov_prop_AR[,,i], ar_counter=arcAR[i]))
    # }
    result = mclapply(1:p, function(i) update_AR(data=Z[i,], current_params=c(PHI[i,], PACF[i, ], SIGMA[i]), proposal_sigma=sdAR[i]*Cov_prop_AR[,,i], ar_counter=arcAR[i]), mc.cores = 1L)
    result = do.call(rbind, result)
    # stopCluster(cl)
    
    # extract results
    # print(result)
    # result = do.call(rbind, result)
    PHI.samples[iter, , ] = PHI = result[, 1:K]
    PACF.samples[iter, , ] = PACF = result[, (K+1):(2*K)]
    SIGMA.samples[iter, ] = SIGMA = result[, 2*K+1]
    arcAR = result[, 2*K+2]
    
    # cat("\n", ll_Z(Z, PHI, SIGMA))
    # 
    # if(iter %% 100 == 0)
    #   cat("\n", sum(SIGMA>1))
    
    # Draw lambda_L
    if(iter == 2500){
      LZ_avg = colMeans(LZ.samples[iter:(iter-499), ])
      lambda_L = max(1e-3, quantile(abs(LZ_avg), probs = 0.8))
      L = (abs(LZ_avg) > lambda_L) * LZ_avg
    }
    # Note: Update in log scale
    if(iter > 2500 && iter %% 20 == 0){
      matL = matL_prop = matrix(0, p, p)
      matL[lower.tri(matL)] = L.samples[iter, ]
      Z = t(matU) %*% matD %*% (diag(p) - matL) %*% Y
      log_accept_prob = - ll_Z(Z, PHI, SIGMA) - log(lambda_L)
      llambda_L_proposal = log(lambda_L) + rnorm(1, 0, sigma_L)
      llambda_L_proposal = log(lambda_L_U) * (llambda_L_proposal > log(lambda_L_U)) + llambda_L_proposal * (llambda_L_proposal <= log(lambda_L_U))
      lambda_L_proposal =  max(1e-3, exp(llambda_L_proposal))
      proposal_L = (abs(LZ.samples[iter, ])>lambda_L_proposal) * LZ.samples[iter, ]
      matL_prop[lower.tri(matL_prop)] = proposal_L
      Z = t(matU) %*% matD %*% (diag(p) - matL_prop) %*% Y
      log_accept_prob = log_accept_prob + ll_Z(Z, PHI, SIGMA) + llambda_L_proposal
      
      if(is.na(log_accept_prob) || is.nan(log_accept_prob))
        log_accept_prob = -Inf
      
      if(log(runif(1)) < log_accept_prob){
        L.samples[iter, ] = L = proposal_L
        arcLL = arcLL + 1
        lambda_L = lambda_L_proposal
        matL = matL_prop
      }
      # cat("\n lambda_L", lambda_L)
      # Draw sigma_T
      ae_L = 0.1 + p * (p-1)/2
      be_L = 0.1 + sum(LZ ^ 2)/2
      
      sigma_T_L = 1/rgamma(1, ae_L, be_L)
    }
    
    # Draw lambda_A
    if(iter == 2500){
      AZ = AZ_avg = colMeans(AZ.samples[iter:(iter-499), ])
      lambda_A = max(1e-6, quantile(abs(AZ_avg), probs = 0.8))
      A = (abs(AZ_avg) > lambda_A) * (abs(AZ_avg) - lambda_A) * sign(AZ_avg)
    }
    # Note: Update in log scale
    if(iter > 2500 && iter %% 30 == 0){
      matA = matA_prop = matrix(0, p, p)
      matA[lower.tri(matA)] = A.samples[iter, ]
      matA = matA - t(matA)
      matU = (diag(p) - matA) %*% solve(diag(p) + matA)
      Z = t(matU) %*% matD %*% (diag(p) - matL) %*% Y
      log_accept_prob = - ll_Z(Z, PHI, SIGMA) - log(lambda_A)
      llambda_A_proposal = log(lambda_A) + rnorm(1, 0, sigma_A)
      llambda_A_proposal = log(lambda_A_U) * (llambda_A_proposal > log(lambda_A_U)) + llambda_A_proposal * (llambda_A_proposal <= log(lambda_A_U))
      lambda_A_proposal =  max(1e-6, exp(llambda_A_proposal))
      proposal_A = (abs(AZ.samples[iter, ])>lambda_A_proposal) * (abs(AZ.samples[iter, ]) - lambda_A_proposal) * sign(AZ.samples[iter, ])
      matA_prop[lower.tri(matA_prop)] = proposal_A
      matA_prop = matA_prop - t(matA_prop)
      matU_prop = (diag(p) - matA_prop) %*% solve(diag(p) + matA_prop)
      Z = t(matU_prop) %*% matD %*% (diag(p) - matL) %*% Y
      log_accept_prob = log_accept_prob + ll_Z(Z, PHI, SIGMA) + llambda_A_proposal
      
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
      
      sigma_T_A = 1/rgamma(1, ae_A, be_A)
    }
    
    # Adapt the Covariance matrix
    # For A
    if(adaptA){
      # if(iter > 5000 && iter <= 7500){
      #    Cov_prop_A = cor(as.matrix(A.samples[(iter-499):iter, ])) + 0.01 * diag((p*(p-1)/2))
      # }
      if(iter > 3500 && iter %% 100 == 0){
        Cov_prop_A = cor(as.matrix(A.samples[(iter-499):iter, ])) + 0.01 * diag((p*(p-1)/2))
      }
    }
    
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
    # For L
    if(iter%%100 == 0 && iter > 1500){
      # a_rate_L = arcL / iter
      a_rate_L = arcL / (iter - 1500)
      if(a_rate_L < 0.45)
        eps_L = eps_L / 10
      else if(a_rate_L > 0.70)
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
      if(iter%%5 == 0 && iter > 1500){
        RMSE = sqrt(mean(c(Omega_pred - Omega0)^2))
        cat("\n ", RMSE)
      }
      if(iter%%200 == 0 && iter > 1500){
        # cat("> ")
        # Omega_hat = (diag(p) - matL) %*% matD^2 %*% t(diag(p) - matL)
        RMSE = sqrt(mean(c(Omega_pred - Omega0)^2))
        cat("\n\n Iterations: ", iter, "| RMSE: ", RMSE, "Acceptance Rates: L:", a_rate_L, " D:", a_rate_D, " A:", a_rate_A, "\n")
        cat("\n Average acceptance rate for AR processes: ", mean(arcAR)/iter)
        cat("\n lambda values for L and A:", lambda_L, lambda_A, "\n")
      }
    }
  }
  return(list(Omega.samples = Omega.samples, 
              A.samples = A.samples, 
              L.samples = L.samples, 
              D.samples = D.samples, 
              PHI.samples = PHI.samples, 
              SIGMA.samples = SIGMA.samples, 
              RMSE_GGM = RMSE_GGM,
              RMSE_GCGM = RMSE_GCGM))
} # End of OUTAR

# try catch implementation
data = generate_data(Tm=args1[[1]], order=args1[[2]], corr = args1[[3]], n_graphs=args1[[4]], n_nodes=args1[[5]], prob=0.1)
# generate_datasets = function(n, args1) {
#   successful_data = list()  # Initialize an empty list to store successful datasets
#   attempts = 0              # Initialize counter for attempts
#   
#   while (length(successful_data) < n) {
#     attempts = attempts + 1
#     print(attempts)
#     print(length(successful_data))
#     # Attempt to generate a dataset
#     result = tryCatch({
#       data = generate_data(Tm=500, order=args1[[2]], corr = args1[[3]], n_graphs=args1[[4]], n_nodes=args1[[5]], prob=0.1)
#       # Return data if successful
#       return(data)
#     }, error = function(e) {
#       # Handle the error: return NULL or log the error
#       message("Error in data generation: ", e$message)
#       return(NULL)
#     })
#     
#     # Check if result is not NULL, meaning data generation was successful
#     if (!is.null(result)) {
#       successful_data = c(successful_data, result)
#     } else {
#       message("Retrying data generation...")
#     }
#   }
#   
#   message("Generated ", length(successful_data), " successful datasets after ", attempts, " attempts.")
#   return(successful_data)
# }
# 
# # Usage
# # args1 = list(Tm = 100, order = 2, corr = 0.5, n_graphs = 10, n_nodes = 50)
# n = 5
# set.seed(999)
# datasets = generate_datasets(5, args1)


# trial run
results = matrix(0, 30, 4)
for(sim in 1:30){
  cat("\n Simulation:= ", sim, "\n")
  data = generate_data(Tm=args1[[1]], order=args1[[2]], corr = args1[[3]], n_graphs=args1[[4]], n_nodes=args1[[5]], prob=0.1)
  Y = data$Y
  # par(mfrow = c(2, 2))
  # for(j in 1:4){
  #   plot(Y[j, ], type="l")
  # }
  # par(mfrow = c(1, 1))
  # print(Y[1:4, 1:20])
  Omega0 = data$Omega0
  tm1 = tic()
  model = OUTAR(Y = Y, K = args1[[6]], n.iter = 10000, verbose = TRUE, adaptA = FALSE, adaptAR = FALSE)
  tm2 = toc()
  Omega.pred = apply(model$Omega.samples[5001:10000,,], c(2, 3), mean)
  results[sim, 1] = model$RMSE_GGM
  results[sim, 2] = model$RMSE_GCGM
  results[sim, 3] = sqrt(mean(c(Omega0 - Omega.pred)^2))
  results[sim, 4] = as.numeric(tm2$toc - tm2$tic)
}

print(round(results, 4))
cat("\n")
print(round(colMeans(results), 4))
write.csv(x=round(colMeans(results), 4), file=paste0("/share/statistics/sghosh27/BayesTS/Final/csvs/result_", as.numeric(args), ".csv"))

# Omega.pred = apply(model$Omega.samples[7501:10000,,], c(2, 3), mean)
# sqrt(mean(c(Omega0 - Omega.pred)^2))
# Omega.pred = apply(model$Omega.samples[10001:15000,,], c(2, 3), mean)
# sqrt(mean(c(Omega0 - Omega.pred)^2))
