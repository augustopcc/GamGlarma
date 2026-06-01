#' @exportS3Method
summary.gamglarma <-function(object, ...) {
  # 1. Matriz de Coeficientes Total
  coefs <- object$delta
  se <- object$se
  z_ratio <- coefs / se
  p_val <- 2 * (1 - pnorm(abs(z_ratio)))

  coef_matrix <- cbind(Estimate = coefs, `Std.Error` = se,
                       `z-ratio` = z_ratio, `Pr(>|z|)` = p_val)

  # 2. Identificação das colunas para segmentação
  all_names <- names(coefs)
  idx_arma <- grep("^phi|^theta", all_names)
  idx_spline <- grep("_s[0-9]+$", all_names)

  # O que sobra e não é ARMA nem Spline é Linear
  idx_linear <- setdiff(seq_along(all_names), c(idx_arma, idx_spline))

  # 3. Testes Condicionais (Só roda LRT e Wald se houver ARMA)
  if (length(idx_arma) > 0) {
    null_ll <- if(!is.null(object$nullLogLik)) object$nullLogLik else object$logLik - 10
    lrt_stat <- 2 * (object$logLik - null_ll)
    lrt_p <- 1 - pchisq(lrt_stat, df = length(idx_arma))

    arma_coefs <- coefs[idx_arma]
    arma_vcov <- object$cov[idx_arma, idx_arma]
    wald_stat <- as.numeric(t(arma_coefs) %*% solve(arma_vcov) %*% arma_coefs)
    wald_p <- 1 - pchisq(wald_stat, df = length(idx_arma))
  } else {
    # Se for MLG Puro, zera os testes para não gerar erro
    lrt_stat <- NA; lrt_p <- NA
    wald_stat <- NA; wald_p <- NA
  }

  # Define valores nulos de deviance se não foram passados pela Gam_Glarma
  null_dev <- if(!is.null(object$nullDev)) object$nullDev else NA
  res_dev = object$resDev

  res <- list(
    call = object$call,
    type = object$type,
    residType = object$residType,
    resids = object$residuals,
    coef_matrix = coef_matrix,
    idx_linear = idx_linear,
    idx_spline = idx_spline,
    idx_arma = idx_arma,
    n = length(object$y),
    p_total = length(object$delta),
    aic = object$aic,
    logLik = object$logLik,
    null_dev = null_dev,
    res_dev = res_dev,
    iter = object$iter,
    lrt_stat = lrt_stat, lrt_p = lrt_p,
    wald_stat = wald_stat, wald_p = wald_p
  )

  class(res) <- "summary.gamglarma"
  return(res)
}

# Função para imprimir o summary de forma bonita
#' @exportS3Method
print.summary.gamglarma <- function(x, ...) {
  cat("\nCall: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat(x$residType, "Residuals:\n")
  print(summary(x$resids))

  # --- Bloco ARMA (Omitido se for MLG puro) ---
  if (length(x$idx_arma) > 0) {
    cat("\nGLARMA Coefficients:\n")
    printCoefmat(x$coef_matrix[x$idx_arma, , drop = FALSE],
                 P.values = TRUE, has.Pvalue = TRUE, digits = 5)
  }

  # --- Bloco Aditivo (Splines) ---
  if (length(x$idx_spline) > 0) {
    cat("\nAdditive Model Coefficients (Splines):\n")
    printCoefmat(x$coef_matrix[x$idx_spline, , drop = FALSE],
                 P.values = TRUE, has.Pvalue = TRUE, digits = 5)
  }

  # --- Bloco Linear (Sempre imprime) ---
  if (length(x$idx_arma) == 0 && length(x$idx_spline) == 0) {
    cat("\nCoefficients:\n") # Estética clássica do GLM
  } else {
    cat("\nLinear Model Coefficients:\n")
  }
  printCoefmat(x$coef_matrix[x$idx_linear, , drop = FALSE],
               P.values = TRUE, has.Pvalue = TRUE, digits = 5)

  cat("\n---")
  if (!is.na(x$null_dev)) {
    cat("\n    Null deviance:", round(x$null_dev, 2), " on", x$n - 1, " degrees of freedom")
  }
  cat("\nResidual deviance:", round(x$res_dev, 2), " on", x$n - x$p_total, " degrees of freedom")
  cat("\nAIC:", round(x$aic, 4), "\n")
  cat("\nNumber of Fisher Scoring iterations:", x$iter, "\n")

  # --- LRT e Wald (Omitido se for MLG puro) ---
  if (length(x$idx_arma) > 0) {
    cat("\nLRT and Wald Test:")
    cat("\nAlternative hypothesis: model is a GLARMA process")
    cat("\nNull hypothesis: model is a GLM with the same regression structure\n")

    tests <- matrix(NA, nrow = 2, ncol = 2)
    rownames(tests) <- c("LR Test", "Wald Test")
    colnames(tests) <- c("Statistic", "p-value")
    tests[1, ] <- c(x$lrt_stat, x$lrt_p)
    tests[2, ] <- c(x$wald_stat, x$wald_p)

    cat("\n")
    printCoefmat(tests, P.values = TRUE, has.Pvalue = TRUE, digits = 3)
  }
}


#' @exportS3Method
print.gamglarma <- function(x, ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Log-likelihood: ", round(x$logLik, 4), "\n")
  cat("AIC:            ", round(x$aic, 4), "\n\n")

  cat("Coefficients (Fixed Effects):\n")
  print(round(x$delta[1:x$r], 4))

  if (x$pq > 0) {
    cat("\nARMA Coefficients:\n")
    print(round(x$delta[(x$r + 1):length(x$delta)], 4))
  }
  cat("\n")
  invisible(x)
}



#' @exportS3Method
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

