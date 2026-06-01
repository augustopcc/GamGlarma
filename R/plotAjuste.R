#' @title Gráfico de Ajuste do Modelo
#' @description Gera um gráfico comparando os valores reais e os valores ajustados (fitted values) ao longo do tempo.
#' @param modelo Objeto retornado pela função Gam_Glarma.
#' @export
plotAjuste <- function(modelo = NULL){
  if (is.null(modelo) || is.null(modelo$residuals)) {
    stop("Forneça um modelo válido que contenha 'fitted.values'.")
  }

  real   <- modelo$y
  ajuste <- modelo$fitted.values
  n      <- length(real)

  if (isTRUE(modelo$is_ts_y)) {
    # Reconstrói a indexação de tempo através dos metadados guardados
    y_ts_temp <- ts(real, start = modelo$ts_start, frequency = modelo$ts_freq)
    x <- as.numeric(time(y_ts_temp))
    label_x <- "Tempo"
  } else {
    x <- 1:n
    label_x <- "Índice"
  }
  oldpar <- par(no.readonly = TRUE)

  # 3. Define a disposição dos gráficos com margens otimizadas (evita erros de tamanho)

  par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))

  ylim_range <- range(c(real, ajuste), na.rm = TRUE)

  plot(x, real,
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

  invisible(modelo)
}
