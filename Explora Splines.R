
# expolorando splines 
BH5 <- readr::read_delim("Fwd_ Dados/BH5.csv", delim = ";", 
                  escape_double = FALSE, locale = readr::locale(decimal_mark = ",", 
                                                         grouping_mark = "."), trim_ws = TRUE)



women$height |> plot()

basis1 <- splines::ns(BH5$RH, df = 1, intercept = F)
basis2 <- splines::bs(BH5$RH, df = 2,degree = 2, intercept = F)

# Ajustando os modelos
mod_ns <- glm(BH5$PM10 ~ basis1,family = Gamma(link = "inverse"))
mod_bs <- glm(BH5$PM10 ~ basis2,family = Gamma(link = "inverse"))

summary(mod_bs)
summary(mod_ns)
# Comparando as previsões finais
fitted(mod_ns)
fitted(mod_bs)


newX <- seq(58, 80, length.out = 51)
plot(newX)
# evaluate the basis at the new data
predict1 <- predict(basis1, newX)
predict2 <- predict(basis2, newX)

plot(predict1)
plot(predict2)


# ==============================================================
plot(BH5$RH, BH5$PM10, 
     pch = 16, col = "darkgray", cex = 1.2,
     main = "Previsões: ns() vs bs()",
     xlab = "Altura (X)", ylab = "Peso (Y)")

# Adiciona a linha do modelo com Spline Natural (Azul)
lines(BH5$RH, fitted(mod_ns), col = "blue", lwd = 3)

# Adiciona a linha do modelo com B-Spline (Vermelha e Tracejada)
lines(BH5$RH, fitted(mod_bs), col = "red", lwd = 3, lty = 2)

# Adiciona uma marcação onde a ns() colocou a sua "dobradiça" oculta
abline(v = attr(basis1, "knots"), col = "blue", lty = 3, lwd = 1.5)

legend("topleft", 
       legend = c("ns (Curva Articulada - 1 nó)", "bs (Parábola Rígida - 0 nós)"),
       col = c("blue", "red"), lwd = 3, lty = c(1, 2), bty = "n")

# =========================================================================
# GRÁFICO 2: A Diferença Visual das Bases (As colunas da Matriz X)
# =========================================================================
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)) # Divide o painel inferior em dois

# Desenha as 2 colunas da matriz ns()
matplot(BH5$RH, basis1, 
        type = "l", lty = 1, lwd = 3, col = c("darkblue", "dodgerblue"),
        main = "Geometria da ns()", xlab = "Rh", ylab = "Valor na Matriz X")
# Linha do nó interno
abline(v = attr(basis1, "knots"), col = "gray", lty = 3, lwd = 2)

# Desenha as 2 colunas da matriz bs()
matplot(BH5$RH, basis2, 
        type = "l", lty = 1, lwd = 3, col = c("darkred", "tomato"),
        main = "Geometria da bs()", xlab = "RH", ylab = "")

# Restaura o layout original
par(mfrow = c(1, 1))



n <- 150
x <- seq(1, 10, length.out = n)
# Criando um sinal com várias curvas para testar a flexibilidade
y_sinal <- sin(x) + 0.1 * x^2 + rnorm(n, sd = 0.5)

# Criar a matriz de entrada (X_in)
X_exemplo <- matrix(x, ncol = 1)
colnames(X_exemplo) <- "Tempo"
pred_linear <- Build_Preditor(Y_in = y_sinal, X_in = X_exemplo, 
                              n_spline = 3, spline_cols = 1)
matplot(x, pred_linear[, 1:3], type = "l", lty = 1, lwd = 2,
        main = "Bases: Grau 1 (Máx 5 col)", ylab = "Valor da Base")
# Gráfico 3: Ajuste Linear (Como o modelo 'enxerga' os dados)
mod_lin <- lm(y_sinal ~ pred_linear)
plot(x, y_sinal, col = "gray", pch = 19, main = "Ajuste Final (Grau 1)")
lines(x, fitted(mod_lin), col = "blue", lwd = 2)




### P4revendo proximo splines
### 

library(splines)

n <- 150
x_new <- seq(10, 11, length.out = 5)

spline <- splines::ns(seq(1, 10, length.out = n), df = 3)
spline
# O jeito correto de aplicar (use newx)
matriz_nova <- predict(spline, newx = x_new)

print(matriz_nova)
j=1
spline_j = splines::ns(seq(1, 10, length.out = n), df = 3)
formulas <-NULL
formulas <- list(formulas,spline_j)
