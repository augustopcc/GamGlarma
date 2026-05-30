parte_distribuicao <- function(y_t,
                               mu_t,
                               eta_t,
                               W_dt,
                               type = "Poisson",
                               link,
                               residuals = "Score",
                               alpha = 1,
                               ...) {
  # ----- parte 1: resíduo e derivada do resíduo -----
  if (type == "Poisson") {
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

  } else if (type == "Normal"){
    if (link == "identity") {
      ll_eta <- alpha * (y_t - mu_t)
      fisher_weight <- alpha # A Fisher da Normal com identidade é 1/sigma^2 constante

      if (residuals == "Score") {
        e_t <- alpha * (y_t - mu_t)
        de_deta <- -alpha
      } else { # Pearson
        # variancia da normal é 1/alpha, então desvio padrão é 1/sqrt(alpha)
        e_t <- (y_t - mu_t) * sqrt(alpha)
        de_deta <- -sqrt(alpha)
      }
    } else if (link == "log") {
      # --- Lógica para LINK LOG (Normal com resposta estritamente positiva) ---
      # mu = exp(eta) -> dmu/deta = mu

      ll_eta <- alpha * (y_t - mu_t) * mu_t
      fisher_weight <- alpha * (mu_t^2)

      if (residuals == "Score") {
        e_t <- alpha * (y_t - mu_t)
        de_deta <- -alpha * mu_t # (Derivada simplificada aproximada para o score)
      } else { # Pearson
        e_t <- (y_t - mu_t) * sqrt(alpha)
        de_deta <- -sqrt(alpha) * mu_t
      }
    }

    # Log-Verossimilhança da Normal
    # Como alpha = 1/variancia, o desvio padrão (sd) é 1/sqrt(alpha)
    ll_contrib <- dnorm(y_t, mean = mu_t, sd = 1/sqrt(alpha), log = TRUE)

    ll_d_contrib <- ll_eta * W_dt
    ll_dd_contrib <- -fisher_weight * (W_dt %*% t(W_dt))
    e_dt <- de_deta * W_dt


  } else if (type == "Gamma") {
    if (link == "log") {
      # --- Lógica para LINK LOG ---
      ll_eta <- alpha * ((y_t / mu_t) - 1)
      fisher_weight <- alpha * (y_t / mu_t)

      if (residuals == "Score") {
        e_t <- (y_t / mu_t) - 1
        de_deta <- -(y_t / mu_t)
      } else {# Pearson

        e_t <- (y_t - mu_t) / (mu_t / sqrt(alpha))
        de_deta <- -sqrt(alpha) * (y_t / mu_t)
      }

    } else if (link == "inverse") {
      # --- Lógica para LINK INVERSO (Canônico) ---
      ll_eta <- alpha * (mu_t - y_t)
      fisher_weight <- alpha * (mu_t^2)

      if (residuals == "Score") {
        e_t <- alpha * (mu_t - y_t)
        de_deta <- -alpha * (mu_t^2)
      } else {
        # Pearson
        e_t <- (y_t - mu_t) / (mu_t / sqrt(alpha))
        de_deta <- -sqrt(alpha) * mu_t
      }
    }

    ll_contrib <- dgamma(y_t,
                         shape = alpha,
                         scale = mu_t / alpha,
                         log = TRUE)
    ll_d_contrib <- ll_eta * W_dt
    ll_dd_contrib <- -fisher_weight * (W_dt %*% t(W_dt))
    e_dt <- de_deta * W_dt
  }

  list(
    e = e_t,
    e_d = e_dt,
    ll = ll_contrib,
    ll_d = ll_d_contrib,
    ll_dd = ll_dd_contrib
  )
}
