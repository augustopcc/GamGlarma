Iteracao_IRLS<- function(delta, Xf, y, offset, type, residuals, alpha,link = NULL,
                         phiLags, thetaLags, p, q, mpq,...) {

  r <- ncol(Xf)
  s <- length(delta)
  # particiona delta
  beta <- delta[1:r]
  idx  <- r
  phi  <- if (p > 0) delta[(idx + 1):(idx + p)] else numeric(0)
  idx  <- idx + p
  theta <- if (q > 0) delta[(idx + 1):(idx + q)] else numeric(0)

  n <- length(y)

  # memoria
  Z  <- numeric(n + mpq)
  e  <- numeric(n + mpq)
  W  <- numeric(n + mpq)
  mu <- numeric(n + mpq)

  Z_d <- matrix(0, nrow = s, ncol = n + mpq)
  W_d <- matrix(0, nrow = s, ncol = n + mpq)
  e_d <- matrix(0, nrow = s, ncol = n + mpq)

  ll    <- 0
  ll_d  <- matrix(0, nrow = s, ncol = 1)
  ll_dd <- matrix(0, nrow = s, ncol = s)

  # loop temporal
  for (t in seq_len(n)) {
    tt <- t + mpq

    # derivada wrt beta
    W_d[, tt] <- 0
    W_d[1:r, tt] <- Xf[t, ]

    # ---------------- AR(p) ----------------
    if (p > 0) {
      Z_d[(r + 1):(r + p), tt] <- Z[tt - phiLags] + e[tt - phiLags]

      for (i in seq_len(p)) {
        lag_i <- phiLags[i]
        Z[tt] <- Z[tt] + phi[i] * (Z[tt - lag_i] + e[tt - lag_i])
        Z_d[, tt] <- Z_d[, tt] +
          phi[i] * (Z_d[, tt - lag_i] + e_d[, tt - lag_i])
      }
    }

    # ---------------- MA(q) ----------------
    if (q > 0) {
      Z_d[(r + p + 1):(r + p + q), tt] <- e[tt - thetaLags]

      for (j in seq_len(q)) {
        lag_j <- thetaLags[j]
        Z[tt] <- Z[tt] + theta[j] * e[tt - lag_j]
        Z_d[, tt] <- Z_d[, tt] +
          theta[j] * e_d[, tt - lag_j]
      }
    }

    # eta e derivada final de eta
    eta_t <- as.numeric(drop(Xf[t, ] %*% beta) + Z[tt] + offset[t])
    W[tt] <- eta_t
    W_d[, tt] <- W_d[, tt] + Z_d[, tt]

    # media
    if (link == "log") {
      mu[tt] <- exp(eta_t)
    } else if (link == "identity"){
      mu[tt] <- eta_t
    } else if (link == "inverse") {
      if (type == "Gamma") {
        if (eta_t < 1e-6) eta_t <- 1e-6 # Trava de segurança para link inverso
        mu[tt] <- 1/eta_t
      } else {
        stop("Link inverso nao suportado para esta distribuicao nesta funcao.")
      }
    }

    # componente da distribuicao
    dp <- parte_distribuicao(
      y_t = y[t],
      mu_t = mu[tt],link = link,
      eta_t = eta_t,
      W_dt = matrix(W_d[, tt], ncol = 1),
      type = type,
      residuals = residuals,
      alpha = alpha
    )

    e[tt] <- dp$e
    e_d[, tt] <- dp$e_d

    ll    <- ll + dp$ll
    ll_d  <- ll_d + dp$ll_d
    ll_dd <- ll_dd + dp$ll_dd
  }

  list(
    ll = ll,
    ll_d = ll_d,
    ll_dd = ll_dd,
    Z = Z,
    W = W,
    mu = mu,
    e = e,
    Z_d = Z_d,
    W_d = W_d,
    e_d = e_d,
    r = r
  )
}
