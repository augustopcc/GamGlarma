# =========================================================================
# PASSO 1: CONFIGURAÇÃO INICIAL E COVARIÁVEL
# =========================================================================
n <- 400
set.seed(42) # Fixa a semente para o resultado ser reprodutível

# Gera a covariável (X1) que vai alimentar o modelo
X1 <- rnorm(n, mean = 0, sd = 1.5)

# =========================================================================
# PASSO 2: A BASE DETERMINÍSTICA (O "X * Beta")
# (Escolha apenas UMA das opções abaixo e comente a outra)
# =========================================================================

# OPÇÃO A: MLG Linear Clássico
# preditor_linear <- 1.0 + 0.5 * X1

# OPÇÃO B: GAM (A Spline com alta variabilidade e intensidade)
preditor_linear <- 1.0 + 2.5 * sin(1.2 * X1) + 1.5 * cos(2.5 * X1)

# =========================================================================
# PASSO 3: INTERVENÇÕES ESTRUTURAIS
# (Rode as linhas abaixo apenas se quiser adicionar esses choques)
# =========================================================================

# Adicionando um PASSO (Quebra de nível a partir do tempo 200)
preditor_linear[200:n] <- preditor_linear[200:n] + 1.2

# Adicionando IMPULSOS (Choques absurdos isolados nos tempos 100 e 300)
preditor_linear[c(100, 300)] <- preditor_linear[c(100, 300)] + 2.5

# =========================================================================
# PASSO 4: PARÂMETROS ARMA E INICIALIZAÇÃO DO ESTADO LATENTE
# (Se quiser testar um MLG puro, mude phi1 e theta1 para 0)
# =========================================================================
phi1   <- 0.45
theta1 <- 0.25

# Cria os vetores vazios para guardar a história do tempo
Y  <- numeric(n)
Z  <- numeric(n)
e  <- numeric(n)
mu <- numeric(n)

# O primeiro dia (Tempo 1) não tem passado para olhar, então é calculado isolado
eta_1 <- preditor_linear[1]
mu[1] <- exp(eta_1)
Y[1]  <- rpois(1, lambda = mu[1])
e[1]  <- (Y[1] - mu[1]) / mu[1] # Resíduo Score (Poisson)

# =========================================================================
# PASSO 5: O LOOP DO TEMPO (A Mágica do GLARMA)
# =========================================================================
for (t in 2:n) {
  # 1. Atualiza a memória (Z) com base na inércia de ontem
  Z[t] <- phi1 * (Z[t-1] + e[t-1]) + theta1 * e[t-1]

  # 2. Soma o clima de hoje (Preditor Linear) com a inércia de ontem (Z)
  eta_t <- preditor_linear[t] + Z[t]

  # 3. Transforma na escala real (Link Log)
  mu[t] <- exp(eta_t)

  # 4. A Natureza age: Sorteia o valor real com base na distribuição Poisson
  Y[t] <- rpois(1, lambda = mu[t])

  # 5. Calcula o erro de hoje (que vai alimentar o loop de amanhã)
  e[t] <- (Y[t] - mu[t]) / mu[t]
}

# =========================================================================
# PASSO 6: EMPACOTAMENTO E VISUALIZAÇÃO
# =========================================================================
dados_simulados <- data.frame(Tempo = 1:n, Y = Y, X1 = X1)

# Plota a série real final que o modelo vai enxergar
plot(dados_simulados$Y, type = "l", col = "black", lwd = 1,
     main = "Série Simulada (DGP)", ylab = "Y", xlab = "Tempo")

# Se quiser ver apenas o formato determinístico oculto por baixo do ruído:
# lines(exp(preditor_linear), col = "red", lwd = 2)
