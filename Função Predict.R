
predict.gamglarma <- function(object,n.ahead = 0, novos_dados = NULL, ...){


  #caso nao haja novos dados e nem range de novas previsão, realiza a previsão em cima dos valores usados pra treino, ou seja retorna o fittedvalues
  if (is.null(novos_dados) && n.ahead == 0) {
    message("Nenhum horizonte futuro ou dados informados. Retornando fitted.values.")
    return(object$fitted.values)
  }

  n_hist <- length(object$y)

  if (is.null(novos_dados) && n.ahead > 0) {
    if (!requireNamespace("forecast", quietly = TRUE)) {
      stop("O pacote 'forecast' é necessário para previsão automática das preditoras.")
    }

    dados_orig <- as.data.frame(object$X_original)
    novos_dados <- data.frame(matrix(ncol = ncol(dados_orig), nrow = n.ahead))
    colnames(novos_dados) <- colnames(dados_orig)

    for (col in colnames(dados_orig)) {
      if (is.numeric(dados_orig[[col]])) {
        modelo_arima <- forecast::auto.arima(dados_orig[[col]])
        novos_dados[[col]] <- as.numeric(forecast::forecast(modelo_arima, h = n.ahead)$mean)
      }
    }
    message(sprintf("Previsão gerada automaticamente por auto.arima para %d passos.", n.ahead))

  } else {
    n.ahead <- nrow(novos_dados)
  }


  X_new <- matrix(nrow = n.ahead, ncol = 0)

  # a) Intercepto
  if ("(Intercept)" %in% colnames(object$X)) {
    X_new <- cbind(X_new, `(Intercept)` = rep(1, n.ahead))
  }

  # b) Variáveis Lineares Básicas
  cols_lineares <- setdiff(object$variaveis_lineares, "(Intercept)")
  if (length(cols_lineares) > 0) {
    X_new <- cbind(X_new, as.matrix(novos_dados[, cols_lineares, drop = FALSE]))
  }

  #Utiliza a estrutura dos splines antigos para ajustar os novos splines nos dados gerados pelo autoarima.
  if (isTRUE(object$ind_spline)) {
    for (var_nome in names(object$Formula_spline)) {
      base_historica <- object$Formula_spline[[var_nome]]
      matriz_spline_futura <- predict(base_historica, newx = novos_dados[[var_nome]])
      colnames(matriz_spline_futura) <- paste0(var_nome, "_s", seq_len(ncol(matriz_spline_futura)))
      X_new <- cbind(X_new, matriz_spline_futura)
    }
  }

  # d) Componentes Determinísticos (Controlados pelos novos indicadores)

  # - Tendência
  if (isTRUE(object$ind_tendencia)) {
    tend_futura <- (n_hist + 1:n.ahead) * (1 / n_hist)
    X_new <- cbind(X_new, Tendencia = tend_futura)
  }

  # - Seno e Cosseno
  if (isTRUE(object$ind_sen_cos)) {
    for (i in object$sen_cos) {

      t_futuro_ciclico <- ((n_hist + 1:n.ahead - 1) %% i) + 1

      sen_futuro <- sin((2 * pi * t_futuro_ciclico) / i)
      cos_futuro <- cos((2 * pi * t_futuro_ciclico) / i)

      matriz_onda <- cbind(sen_futuro, cos_futuro)
      colnames(matriz_onda) <- c(paste0("sen", i), paste0("cos", i))

      X_new <- cbind(X_new, matriz_onda)
    }
  }

  # - Impulso (Futuro padrão é zero para todos os impulsos passados)
  if (isTRUE(object$ind_impulso)) {
    cols_impulso <- grep("Impulso_", colnames(object$X), value = TRUE)
    matriz_impulso_futura <- matrix(0, nrow = n.ahead, ncol = length(cols_impulso))
    colnames(matriz_impulso_futura) <- cols_impulso
    X_new <- cbind(X_new, matriz_impulso_futura)
  }

  # - Passo (Futuro mantém o último platô conhecido)
  if (isTRUE(object$ind_passo)) {
    cols_passo <- grep("Passo_", colnames(object$X), value = TRUE)
    ultimo_passo <- object$X[n_hist, cols_passo, drop = FALSE]
    matriz_passo_futura <- matrix(rep(ultimo_passo, each = n.ahead),
                                  nrow = n.ahead, ncol = length(cols_passo))
    colnames(matriz_passo_futura) <- cols_passo
    X_new <- cbind(X_new, matriz_passo_futura)
  }


  X_new <- X_new[, colnames(object$X), drop = FALSE]

  # =========================================================================
  # 5. CÁLCULO DO PREDITOR LINEAR FIXO (X * Beta)

  idx_fixed <- 1:ncol(object$X)
  betas_fixos <- object$delta[idx_fixed]

  eta_futuro <- as.numeric(X_new %*% betas_fixos)

  # =========================================================================
  # 6. INCORPORAÇÃO DO GLARMA E RETORNO

  idx_phi   <- grep("^phi", names(object$delta))
  idx_theta <- grep("^theta", names(object$delta))

  phi_vals   <- object$delta[idx_phi]
  theta_vals <- object$delta[idx_theta]

  p <- length(object$phiLags)
  q <- length(object$thetaLags)

  # 6.2 Prepara os vetores com o histórico e o futuro zerado
  # Regras 3 e 4 do PDF: Erros passados são os resíduos; erros futuros são 0.
  e_total <- c(object$residuals, rep(0, n.ahead))

  # Regras 1 e 2 do PDF: Estados passados vêm do modelo; futuros começam em 0 e serão iterados.
  Z_total <- c(object$Z, rep(0, n.ahead))

  # 6.3 O Loop de Decaimento Temporal GLARMA
  if (p > 0 || q > 0) {
    for (k in 1:n.ahead) {

      t_atual <- n_hist + k

      # Parte Autoregressiva (AR): phi * (Z + e)
      soma_ar <- 0
      if (p > 0) {
        for (i in 1:p) {
          lag <- object$phiLags[i]

          # AQUI ESTÁ A CORREÇÃO DO PDF:
          # phi_i * ( E[Z_{t-i}] + E[\epsilon_{t-i}] )
          soma_ar <- soma_ar + phi_vals[i] * (Z_total[t_atual - lag] + e_total[t_atual - lag])
        }
      }

      # Parte de Médias Móveis (MA): theta * e
      soma_ma <- 0
      if (q > 0) {
        for (j in 1:q) {
          lag <- object$thetaLags[j]

          # theta_j * E[\epsilon_{t-j}]
          soma_ma <- soma_ma + theta_vals[j] * e_total[t_atual - lag]
        }
      }

      # Atualiza a previsão iterativa de Z_t (Regra 2 do PDF)
      Z_total[t_atual] <- soma_ar + soma_ma
    }
  }

  Z_futuro <- Z_total[(n_hist + 1):(n_hist + n.ahead)]
  eta_futuro <- eta_futuro + Z_futuro

  # Reverte a função de ligação
  mu_futuro <- if(object$type == "Gamma" && object$link == "inverse") {
    1 / eta_futuro
  } else if (object$link == "log") {
    exp(eta_futuro)
  } else {
    eta_futuro # identity
  }

  # Restaura calendário temporal se existir
  if (isTRUE(object$is_ts_y)) {
    fim_treino <- object$ts_end
    freq <- object$ts_freq

    # O início do futuro é exatamente 1 período após o fim do treino
    start_futuro <- fim_treino[1] + (fim_treino[2]) / freq

    # Converte o vetor numérico de volta para o formato de calendário
    mu_futuro <- ts(mu_futuro, start = start_futuro, frequency = freq)
  }

  return(mu_futuro)

}
