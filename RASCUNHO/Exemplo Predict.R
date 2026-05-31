

# Carrega os dados reais
BH5 <- readr::read_delim("Fwd_ Dados/BH5.csv", delim = ";",
                  escape_double = FALSE, locale = readr::locale(decimal_mark = ",",
                                                         grouping_mark = "."),
                  trim_ws = TRUE)

# Definição do horizonte de previsão
h_prev <- 12
n_total <- nrow(BH5)
n_treino <- n_total - h_prev

# Divisão temporal (Treino vs Teste)
dados_treino <- BH5[1:n_treino, ]
dados_teste  <- BH5[(n_treino + 1):n_total, ]

# O Y real do futuro para podermos calcular o erro no final
y_real_futuro <- dados_teste$PM10

# =========================================================================
# 2. AJUSTE DO MODELO APENAS NOS DADOS DE TREINO

modelo_treino <- Gam_Glarma(
  formula = PM10 ~ RH + NO,
  data = dados_treino,         # Usando apenas o passado
  ts_start = c(2007, 1),
  ts_frequency = 12,
  impulso = 1,
  sen_cos = c(6),
  type = "Gamma",
  link = "log",
  method = "FS",
  residuals = "Score",
  phiLags = NULL,
  thetaLags = 1,               # MA(1)
  n_spline = 2,
  auto_spline = T,
  spline_cols = "NO",          # Nome direto, sem o "X"
  maxit = 50,
  grad = 1e-6,
  patience = 3,
  trace = FALSE
)

# 3. CENÁRIO 1: PREVISÃO CEGA (O MODELO ADIVINHA O X VIA ARIMA)
cat("\n[2] Executando Cenário 1: Previsão Cega (auto.arima para X)...\n")

prev_obj_cenario1 <- predict(modelo_treino, n.ahead = h_prev)

# Extrai a série temporal prevista
Y_hat_cenario1 <- prev_obj_cenario1

# 4. CENÁRIO 2: PREVISÃO INFORMADA (CONHECEMOS O CLIMA/X DO FUTURO)

# Passamos a matriz/dataframe do futuro. O n.ahead é descoberto automaticamente.
prev_obj_cenario2 <- predict(modelo_treino, novos_dados = dados_teste)

# Extrai a série temporal prevista
Y_hat_cenario2 <- prev_obj_cenario2

# =========================================================================
# 5. AVALIAÇÃO DE DESEMPENHO (MÉTRICAS DE ERRO)
# =========================================================================
# EQM: Erro Quadrático Médio
eqm_c1 <- mean((y_real_futuro - as.numeric(Y_hat_cenario1))^2)
eqm_c2 <- mean((y_real_futuro - as.numeric(Y_hat_cenario2))^2)

cat("RESULTADOS DA VALIDAÇÃO OUT-OF-SAMPLE (12 Meses à frente)\n")
cat(sprintf("EQM Cenário 1 (Previsão Cega de X): %.4f\n", eqm_c1))
cat(sprintf("EQM Cenário 2 (Conhecendo o X real): %.4f\n", eqm_c2))


# 6. VISUALIZAÇÃO GRÁFICA DA COMPARAÇÃO
# Para alinhar as datas no gráfico, transformamos o Y real num objeto ts
start_teste <- start(Y_hat_cenario1)
y_real_ts <- ts(y_real_futuro, start = start_teste, frequency = 12)

# Junta tudo do Treino + Teste para um plot contínuo do último ano de treino e o futuro
plot(y_real_ts, type = "l", lwd = 2, col = "black",
     main = "Validação da Previsão GAM-GLARMA: PM10",
     ylab = "Concentração PM10", xlab = "Tempo",
     ylim = range(c(y_real_futuro, Y_hat_cenario1, Y_hat_cenario2)))

# Linha do Cenário 1
lines(Y_hat_cenario1, col = "red", lwd = 2, lty = 2)

# Linha do Cenário 2
lines(Y_hat_cenario2, col = "blue", lwd = 2, lty = 2)

legend("topleft",
       legend = c("Real", "Cenário 1 (X Previsto por ARIMA)", "Cenário 2 (X Real Fornecido)"),
       col = c("black", "red", "blue"),
       lty = c(1, 2, 2), lwd = 2, bty = "n")
