# FUnções
# 
# 


#### COntroi a matriz X, aqui vou definir se vai ter splines ou Não


## ============================================================
## BLOCO A - MATRIZ X (Regressao fixa: linear ou spline)
## ============================================================

build_design_matrix <- function(
    X,
    X_imp = NULL, # Identificação de Outlier
    nT,
    use_splines = FALSE,
    df_spline = 4
) {
  
  if (!use_splines) {
    ## Caso GLARMA linear
    Matriz_X <- cbind(
      Intercept = 1,
      X_imp[1:nT, , drop = FALSE]
    )
    return(Matriz_X)
  }
  
  ## Caso GAM-GLARMA
  library(splines)
  
  B_list <- lapply(seq_len(ncol(X_cont)), function(j) {
    bs(X[, j], df = df_spline)
  })
  
  B <- do.call(cbind, B_list)
  
  Matriz_X <- cbind(
    Intercept = 1,
    X_imp[1:nT, , drop = FALSE],
    B[1:nT, , drop = FALSE]
  )
  
  return(Matriz_X)
}




## ============================================================
## BLOCO B - Nucleo Temporal
## ============================================================

# (re)particiona delta atual
r <- ncol(X)
beta  <- delta[1:r]
idx   <- r
phi   <- if (p > 0) delta[(idx+1):(idx+p)] else numeric(0); idx <- idx + p
theta <- if (q > 0) delta[(idx+1):(idx+q)] else numeric(0)

## zera acumuladores desta iteração
Z[] <- 0; e[] <- 0; W[] <- 0; mu[] <- 0
Z.d[,] <- 0; W.d[,] <- 0; e.d[,] <- 0

ll    <- 0
ll.d  <- matrix(0, nrow = s, ncol = 1)
ll.dd <- matrix(0, nrow = s, ncol = s)
## -------------------------------------------------------------------------------- ##
## ---------------- laço temporal no TREINO: t = 1..nT ---------------
## -------------------------------------------------------------------------------- ##
for (t in 1:nT) {
  tt <- t + mpq
  
  ## d eta / d beta começa em X[t,]
  W.d[, tt] <- 0
  W.d[1:r, tt] <- X[t, ]
  
  ## AR(4): contribuições em Z e em derivadas
  if (p > 0) {
    Z.d[(r+1):(r+p), tt] <- Z[tt - phiLags] + e[tt - phiLags]
    for (i in seq_len(p)) {
      lag_i    <- phiLags[i]
      Z[tt]     <- Z[tt] + phi[i] * (Z[tt - lag_i] + e[tt - lag_i])
      Z.d[, tt] <- Z.d[, tt] + phi[i] * (Z.d[, tt - lag_i] + e.d[, tt - lag_i])
    }
  }
  
  ## MA(0):
  
  if (q > 0) {
    Z.d[(r+p+1):(r+p+q), tt] <- e[tt - thetaLags]
    for (j in seq_len(q)) {
      lag_j <- thetaLags[j]
      Z[tt] <- Z[tt] + theta[j] * e[tt - lag_j]
      Z.d[,tt] <- Z.d[,tt] + theta[j] * e.d[,tt - lag_j]
    }
  }
  ## -------------------------------------------------------------------------------- ##
  ## ----------- parte dependente de distribuição -----------------
  ## -------------------------------------------------------------------------------- ##
  
  ## eta, mu e resíduo score (Poisson)
  eta_t  <- as.numeric( drop(X[t, ] %*% beta) + Z[tt] + offset_vec[t] )
  W[tt]  <- eta_t
  
  ####Ligação exponencial
  mu[tt] <- exp(eta_t)
  ####Residuo score
  e[tt]  <- (yT[t] - mu[tt]) / mu[tt]
  
  ## derivadas finais
  W.d[, tt] <- W.d[, tt] + Z.d[, tt]
  e.d[, tt] <- -(1 + e[tt]) * W.d[, tt]
  
  ## loglik, escore e info esperada, mudar a verossimilhança de acordo com a distribuição
  ll    <- ll + (yT[t] * eta_t - mu[tt] - lfactorial(yT[t]))
  ll.d  <- ll.d + (yT[t] - mu[tt]) * matrix(W.d[, tt], ncol = 1)
  ll.dd <- ll.dd - mu[tt] * (W.d[, tt, drop = FALSE] %*% t(W.d[, tt, drop = FALSE]))
}


###  
###  
###  
###  
###  
###  
###  
###  
###  
###  
###  
glarma_temporal_step <- function(
    t, tt,
    X, beta,
    Z, e, Z_d, e_d, W_d,
    phi, phiLags,
    theta, thetaLags,
    p, q, r
) {
  
  ## Derivada inicial de eta wrt beta
  W_d[, tt] <- 0
  W_d[1:r, tt] <- X[t, ]
  
  ## --------- AR(p) ----------
  if (p > 0) {
    Z_d[(r+1):(r+p), tt] <- Z[tt - phiLags] + e[tt - phiLags]
    
    for (i in seq_len(p)) {
      lag_i <- phiLags[i]
      Z[tt] <- Z[tt] + phi[i] * (Z[tt - lag_i] + e[tt - lag_i])
      Z_d[, tt] <- Z_d[, tt] +
        phi[i] * (Z_d[, tt - lag_i] + e_d[, tt - lag_i])
    }
  }
  
  ## --------- MA(q) ----------
  if (q > 0) {
    Z_d[(r+p+1):(r+p+q), tt] <- e[tt - thetaLags]
    
    for (j in seq_len(q)) {
      lag_j <- thetaLags[j]
      Z[tt] <- Z[tt] + theta[j] * e[tt - lag_j]
      Z_d[, tt] <- Z_d[, tt] +
        theta[j] * e_d[, tt - lag_j]
    }
  }
  
  list(
    Z   = Z,
    Z_d = Z_d,
    W_d = W_d
  )
}


dist_poisson_score <- function(y_t, mu_t, eta_t, W_dt) {
  
  e_t <- (y_t - mu_t) / mu_t
  e_dt <- -(1 + e_t) * W_dt
  
  ll <- y_t * eta_t - mu_t - lfactorial(y_t)
  ll_d <- (y_t - mu_t) * W_dt
  ll_dd <- -mu_t * (W_dt %*% t(W_dt))
  
  list(
    e     = e_t,
    e_d   = e_dt,
    ll    = ll,
    ll_d  = ll_d,
    ll_dd = ll_dd
  )
}


################# -Estrutura exemplo de uso- #################3
 

for (t in 1:nT) {
  tt <- t + mpq
  
  ## --- BLOCO TEMPORAL ---
  tmp <- glarma_temporal_step(
    t, tt,
    X, beta,
    Z, e, Z.d, e.d, W.d,
    phi, phiLags,
    theta, thetaLags,
    p, q, r
  )
  
  Z   <- tmp$Z
  Z.d <- tmp$Z_d
  W.d <- tmp$W_d
  
  ## Previsor
  eta_t <- drop(X[t, ] %*% beta) + Z[tt] + offset_vec[t]
  mu_t  <- exp(eta_t)
  
  W.d[, tt] <- W.d[, tt] + Z.d[, tt]
  
  ## --- BLOCO DISTRIBUICAO ---
  dist <- dist_poisson_score(
    y_t  = yT[t],
    mu_t = mu_t,
    eta_t = eta_t,
    W_dt = W.d[, tt]
  )
  
  e[tt]     <- dist$e
  e.d[, tt] <- dist$e_d
  
  ll    <- ll + dist$ll
  ll.d  <- ll.d + matrix(dist$ll_d, ncol = 1)
  ll.dd <- ll.dd + dist$ll_dd
}