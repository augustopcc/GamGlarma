#' @title Diagnóstico de Resíduos
#' @description Gera gráficos e testes estatísticos (Box-Pierce e Shapiro-Wilk) para os resíduos do modelo, além de recomendar impulsos.
#' @param modelo Objeto retornado pela função Gam_Glarma.
#' @export
plotResiduos <- function(modelo = NULL) {
  if (is.null(modelo) || is.null(modelo$residuals)) {
    stop("Forneça um modelo válido que contenha 'residuals'.")
  }

  oldpar <- par(no.readonly = TRUE)

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  residuos <-  modelo$residuals

  if (!is.null(modelo$y) && inherits(modelo$y, "ts")) {
    eixo_x <- as.numeric(time(modelo$y))
    label_x <- "Tempo"
  } else {
    eixo_x <- seq_along(residuos)
    label_x <- "Índice"
  }

  plot(
    y = residuos,
    x = eixo_x,
    type = "p",# "p" para pontos (padrão)
    pch = 16,# Preenchimento dos pontos (bolinhas sólidas)
    col = "black",# Cor dos pontos
    main = "Resíduos",
    xlab = label_x,
    ylab = ""
  )
  abline(
    h = 0,# h indica linha horizontal, posicionada no y = 0
    col = "red",
    lwd = 1,# Espessura da linha
    lty = 2# Tipo da linha (2 = tracejada))
  )
  hist(residuos)
  abline(
    v = 0,# h indica linha horizontal, posicionada no y = 0
    col = "red",
    lwd = 1,
    lty = 2,
    main = "Histograma"
  )
  acf(residuos, main = "ACF", ylab = "")
  pacf(residuos, main = "PACF", ylab = "")


  on.exit(par(oldpar))

  ### Printa os testes###

  teste_bp <- Box.test(residuos, lag = 12, type = "Box-Pierce")
  cat("Teste de Box-Pierce (Independência / Lag = 12):\n")
  cat(sprintf(" Estatística X-squared = %.4f | p-valor = %.4f\n",
              teste_bp$statistic, teste_bp$p.value))

  # Teste de Normalidade de Shapiro-Wilk
  # Usamos try() caso todos os resíduos sejam iguais ou N seja muito grande/pequeno
  sw_test <- try(shapiro.test(residuos), silent = TRUE)
  cat("\nTeste de Shapiro-Wilk (Normalidade):\n")
  if (!inherits(sw_test, "try-error")) {
    cat(sprintf(" Estatística W = %.4f | p-valor = %.4f\n",
                sw_test$statistic, sw_test$p.value))
  } else {
    cat(" (Não foi possível calcular o teste de normalidade pelo Shapiro Wilk)\n")
  }

  ### Verifica Outlier ###
  res_padronizados <- as.numeric(scale(residuos))

  idx_outliers <- which(abs(res_padronizados) > 3)

  qtd_novos_outliers <- length(idx_outliers)

  # Faz a leitura do modelo para ver quantos impulsos já foram criados
  qtd_impulsos_atuais <- length(grep("^Impulso_", names(modelo$delta)))

  total_recomendado <- qtd_novos_outliers - qtd_impulsos_atuais

  if (total_recomendado > 0) {
    cat(sprintf(" Foram encontrados %d ponto(s) com desvio maior que %d desvios padrões.\n",
                length(idx_outliers), 3))
    cat(" Índice(s) na série temporal: ", paste(idx_outliers, collapse = ", "), "\n")
    cat(sprintf("\n Recomenda-se adicionar %d 'impulso' para melhor ajuste (considerando se já existe Impulso).",
                length(total_recomendado)))
  }

}
