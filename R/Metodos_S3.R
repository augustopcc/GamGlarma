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


plotAjuste <- function(modelo = NULL){
  if (is.null(modelo) || is.null(modelo$residuals)) {
    stop("Forneça um modelo válido que contenha 'fitted.values'.")
  }

  real   <- modelo$y
  ajuste <- modelo$fitted.values
  n      <- length(real)

  if (inherits(real, "ts")) {
    x <- as.numeric(time(real))
    label_x <- "Tempo"
  } else {
    x <- 1:n
    label_x <- "Índice"
  }

  oldpar <- par(no.readonly = TRUE)

  # 3. Define a disposição dos gráficos com margens otimizadas (evita erros de tamanho)

  par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))

  ylim_range <- range(c(real, ajuste), na.rm = TRUE)

  grafico <- plot(x, real,
                  type = "l",
                  col = "black",
                  lwd = 1.5,
                  ylim = ylim_range,
                  main = "Ajuste do Modelo no Tempo",
                  xlab = label_x,
                  ylab = "Variável Resposta (Y)")

  # Adiciona a linha do modelo estimado
  lines(x, ajuste, col = "blue", lwd = 1.5)

  # Adiciona legenda limpa sem borda indesejada
  legend("topright",
         legend = c("Real", "Ajustado"),
         col = c("black", "blue"),
         lty = 1,
         lwd = 1.5,
         bty = "n")

  on.exit(par(oldpar))

  return(grafico)
}


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
