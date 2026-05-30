estimacao_distribuicao <- function(y_t,
                       mu_t,
                       eta_t,
                       W_dt,
                       type = "Poi",
                       residuals = "Score",
                       alpha = 1) {
  
  # ----- parte 1: resíduo e derivada do resíduo -----
  if (type == "Poi") {
    if (residuals == "Score") {
      e_t  <- (y_t - mu_t) / mu_t
      de_deta <- -(1 + e_t)
    } else if (residuals == "Pearson") {
      e_t <- (y_t - mu_t) / sqrt(mu_t)
      # de/deta = -0.5 * (y + mu)/sqrt(mu)
      de_deta <- -0.5 * (y_t + mu_t) / sqrt(mu_t)
    }
    
    e_dt <- de_deta * W_dt
    
    # ----- parte 2: ll, score e Fisher (Poisson) -----
    ll_contrib   <- y_t * eta_t - mu_t - lfactorial(y_t)
    ll_d_contrib <- (y_t - mu_t) * W_dt
    
    # Fisher esperada em eta para Poisson = mu_t
    fisher_weight <- mu_t
    ll_dd_contrib <- -fisher_weight * (W_dt %*% t(W_dt))
    
  } else if (type == "Gamma") {
    # Gamma com link log e shape alpha fixo
    # media = mu_t, var = mu_t^2 / alpha
    
    if (residuals == "Score") {
      e_t <- (y_t / mu_t) - 1
      de_deta <- -(1 + e_t)
    } else if (residuals == "Pearson") {
      e_t <- (y_t - mu_t) / (mu_t / sqrt(alpha))
      # e_t = sqrt(alpha) * (y/mu - 1)
      de_deta <- -(sqrt(alpha) + e_t)
    }
    
    e_dt <- de_deta * W_dt
    
    # ll completa via dgamma (mais segura)
    # shape = alpha, scale = mu/alpha  => mean = alpha*(mu/alpha)=mu
    ll_contrib <- dgamma(y_t,
                         shape = alpha,
                         scale = mu_t / alpha,
                         log = TRUE)
    
    # score wrt eta = alpha*(y/mu - 1)
    ll_eta <- alpha * ((y_t / mu_t) - 1)
    ll_d_contrib <- ll_eta * W_dt
    
    # Fisher wrt eta = alpha
    fisher_weight <- alpha
    ll_dd_contrib <- -fisher_weight * (W_dt %*% t(W_dt))
  }
  
  list(
    e = e_t,
    e_d = e_dt,
    ll = ll_contrib,
    ll_d = ll_d_contrib,
    ll_dd = ll_dd_contrib
  )
}