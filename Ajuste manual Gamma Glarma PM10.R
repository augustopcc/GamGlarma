#tentnado glarma com Gamma no na base BH5


library(readr)
BH5 <- read_delim("Fwd_ Dados/BH5.csv", delim = ";", 
                  escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                         grouping_mark = "."), trim_ws = TRUE)


library(corrplot)

par(mfrow=c(1,1))
M <- cor(BH5[2:15]) # Calcula a matriz de correlação de Pearson

# Gerando o gráfico de matriz
corrplot(M, method = "circle", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)



#CO#
par(mfrow=c(2,2))
serieBH5 <- ts(BH5[,3],start=1970,frequency=12) #confirmar ano que ocmeçou
n <- length(serieBH5)   # Tamanho da s�rie
ts.plot(serieBH5, xlab = NULL,ylab="Numero de Casos de Poliomielite") 
acf(serieBH5)
pacf(serieBH5)
hist(serieBH5)

#PM10#
par(mfrow=c(2,2))
serieBH5 <- ts(BH5[,4],start=1970,frequency=12) #confirmar ano que ocmeçou
n <- length(serieBH5)   # Tamanho da s�rie
ts.plot(serieBH5, xlab = NULL,ylab="Numero de Casos de Poliomielite") 
acf(serieBH5)
pacf(serieBH5)
hist(serieBH5)

#NO#
par(mfrow=c(2,2))
serieBH5 <- ts(BH5[,5],start=1970,frequency=12) #confirmar ano que ocmeçou
n <- length(serieBH5)   # Tamanho da s�rie
ts.plot(serieBH5, xlab = NULL,ylab="Numero de Casos de Poliomielite") 
acf(serieBH5)
pacf(serieBH5)
hist(serieBH5)



######### - Seguindo com o PM10 - #######

serieBH5 <- ts(BH5[,4],frequency=12) #confirmar ano que ocmeçou
n <- length(serieBH5)   # Tamanho da s�rie

# Valores para serem informados
H <- 12   		# Numero de previsoes
nT <- n-H		# Numero de observacoes para o conjunto treinamento
alfa <- 0.05  	# Valor de alfa para o intervalo de previsao
Real <- serieBH5[(n-H+1):n]	# Valores reais para previsao

X_MLG<-cbind(as.matrix(BH5[,c(6,7,8,10,15)]))

MLG <- glm(serieBH5[1:nT] ~ X_MLG[1:nT,], family = Gamma(link = "inverse"))
summary(MLG)

##Adicionar as funções impulsos e passo para pegar as outliers e remover a "tendencia" que nao existe.
##Adicionar função que cria a tendencia linear e quadratica.

plot(x = BH5$NOX,BH5$PM10)


########## AJUSTE MANUAL ############


serieBH5 <- ts(BH5[,4],frequency=12) #confirmar ano que ocmeçou
n <- length(serieBH5)   # Tamanho da s�rie

# Valores para serem informados
H <- 12   		# Numero de previsoes
nT <- n-H		# Numero de observacoes para o conjunto treinamento
alfa <- 0.05  	# Valor de alfa para o intervalo de previsao
Real <- serieBH5[(n-H+1):n]	# Valores reais para previsao

X_MLG<-cbind(as.matrix(BH5[,c(5,6,8,10,15)]))

###############################################################
# 2) Preparação da amostra de treino
###############################################################

y  <- serieBH5
yT <- y[1:nT]

###############################################################
# 3) Matriz X (com ou sem splines)
###############################################################

splines <- FALSE

if (splines == FALSE) {
  
  # Sem splines
  X <- cbind(Intercept = 1, X_MLG[1:nT, ])
  r <- ncol(X)
  
} else {
  
  # Com splines
  library(splines)
  
  X_cont <- as.matrix(BH5[, c(5, 6, 8, 10, 15)])
  
  df_spline <- 4
  
  B_list <- lapply(seq_len(ncol(X_cont)), function(j) {
    bs(X_cont[, j], df = df_spline)
  })
  
  B <- do.call(cbind, B_list)
  
  X <- cbind(
    Intercept = 1,
    B[1:nT, ]
  )
  
  r <- ncol(X)
}

###############################################################
# 4) Offset
###############################################################

offset_vec <- rep(0, nT)

###############################################################
# 5) GLM inicial - Gamma com link inverso
###############################################################

# glm.fit é mais seguro porque X já contém intercepto
MLG <- glm.fit(
  x = X,
  y = yT,
  family = Gamma(link = "inverse"),
  offset = offset_vec
)

Betas_iniciais <- as.numeric(MLG$coefficients)
Betas_iniciais[is.na(Betas_iniciais)] <- 0

# estima o shape alpha a partir da dispersão do GLM
# var(Y) = phi * mu^2  => alpha = 1/phi
MLG_aux <- glm(
  yT ~ X[, -1],
  family = Gamma(link = "inverse"),
  offset = offset_vec
)

alpha_shape <- 1 / summary(MLG_aux)$dispersion

cat("Shape alpha estimado via GLM:", alpha_shape, "\n")

###############################################################
# 6) Especificação GLARMA
###############################################################

phiLags   <- c(1L)          # AR(1)
thetaLags <- integer(0)     # MA(0)

p <- length(phiLags)
q <- length(thetaLags)

# parâmetros iniciais
beta  <- Betas_iniciais
phi   <- if (p > 0) rep(0.1, p) else numeric(0)
theta <- if (q > 0) rep(0.1, q) else numeric(0)

delta <- c(beta, phi, theta)

names(delta) <- c(
  colnames(X),
  if (p > 0) paste0("phi", phiLags) else NULL,
  if (q > 0) paste0("theta", thetaLags) else NULL
)

###############################################################
# 7) Memória dos lags
###############################################################

mpq <- if (p + q > 0) {
  max(c(
    if (p > 0) phiLags else 0L,
    if (q > 0) thetaLags else 0L
  ))
} else {
  0L
}

###############################################################
# 8) Objetos do algoritmo
###############################################################

Z  <- numeric(nT + mpq)
e  <- numeric(nT + mpq)
W  <- numeric(nT + mpq)
mu <- numeric(nT + mpq)

s <- length(delta)

Z.d <- matrix(0, nrow = s, ncol = nT + mpq)
W.d <- matrix(0, nrow = s, ncol = nT + mpq)
e.d <- matrix(0, nrow = s, ncol = nT + mpq)

tol <- 1e-6
iter <- 0
###############################################################
# 9) UMA ITERAÇÃO de Fisher Scoring
###############################################################
iter <- iter + 1
# re-particiona delta atual
r <- ncol(X)
beta  <- delta[1:r]
idx   <- r
phi   <- if (p > 0) delta[(idx + 1):(idx + p)] else numeric(0)
idx   <- idx + p
theta <- if (q > 0) delta[(idx + 1):(idx + q)] else numeric(0)

# zera acumuladores
Z[]  <- 0
e[]  <- 0
W[]  <- 0
mu[] <- 0

Z.d[,] <- 0
W.d[,] <- 0
e.d[,] <- 0

ll    <- 0
ll.d  <- matrix(0, nrow = s, ncol = 1)
ll.dd <- matrix(0, nrow = s, ncol = s)

###############################################################
# 10) Loop temporal
###############################################################

for (t in 1:nT) {
  tt <- t + mpq
  
  ## derivada wrt beta
  W.d[, tt] <- 0
  W.d[1:r, tt] <- X[t, ]
  
  ## ---------------- AR(p) ----------------
  if (p > 0) {
    Z.d[(r + 1):(r + p), tt] <- Z[tt - phiLags] + e[tt - phiLags]
    
    for (i in seq_len(p)) {
      lag_i <- phiLags[i]
      
      Z[tt] <- Z[tt] + phi[i] * (Z[tt - lag_i] + e[tt - lag_i])
      
      Z.d[, tt] <- Z.d[, tt] +
        phi[i] * (Z.d[, tt - lag_i] + e.d[, tt - lag_i])
    }
  }
  
  ## ---------------- MA(q) ----------------
  if (q > 0) {
    Z.d[(r + p + 1):(r + p + q), tt] <- e[tt - thetaLags]
    
    for (j in seq_len(q)) {
      lag_j <- thetaLags[j]
      
      Z[tt] <- Z[tt] + theta[j] * e[tt - lag_j]
      
      Z.d[, tt] <- Z.d[, tt] +
        theta[j] * e.d[, tt - lag_j]
    }
  }
  
  #############################################################
  # Parte dependente da distribuição - Gamma(link = inverse)
  #############################################################
  
  ## eta_t = 1/mu_t
  eta_t <- as.numeric(drop(X[t, ] %*% beta) + Z[tt] + offset_vec[t])
  
  ## proteção numérica: ligação inversa exige eta_t > 0
  if (!is.finite(eta_t) || eta_t <= 0) {
    eta_t <- 1e-8
  }
  
  W[tt] <- eta_t
  
  ## derivada final de eta
  W.d[, tt] <- W.d[, tt] + Z.d[, tt]
  
  ## média condicional
  mu[tt] <- 1 / eta_t
  
  ## resíduo score escalado pela Fisher
  ## e_t = (mu_t - y_t)/mu_t^2
  e[tt] <- (mu[tt] - yT[t]) / (mu[tt]^2)
  
  ## derivada do resíduo:
  ## e_t = eta_t - y_t * eta_t^2
  ## de/deta = 1 - 2 y_t eta_t
  e.d[, tt] <- (1 - 2 * yT[t] * eta_t) * W.d[, tt]
  
  ## log-verossimilhança Gamma
  ## shape = alpha_shape, scale = mu/alpha_shape
  ll <- ll + dgamma(
    yT[t],
    shape = alpha_shape,
    scale = mu[tt] / alpha_shape,
    log = TRUE
  )
  
  ## escore wrt delta:
  ## dℓ/dη = alpha_shape * (mu_t - y_t)
  ll.d <- ll.d + alpha_shape * (mu[tt] - yT[t]) * matrix(W.d[, tt], ncol = 1)
  
  ## informação de Fisher:
  ## I_t = alpha_shape * mu_t^2
  ll.dd <- ll.dd - alpha_shape * (mu[tt]^2) *
    (W.d[, tt, drop = FALSE] %*% t(W.d[, tt, drop = FALSE]))
}


##############################################################
# 11) Passo de Fisher Scoring
###############################################################

A <- -ll.dd
step <- try(solve(A, ll.d), silent = TRUE)

if (inherits(step, "try-error")) {
  stop("Info(Fisher) nao inversivel; ajuste palpites, p/q ou regularize.")
}

delta_old <- delta
delta <- as.numeric(delta + step)
names(delta) <- names(delta_old)

###############################################################
# 12) Relatório
###############################################################

crit <- max(abs(ll.d))

K <- length(delta)
AIC_atual <- -2 * ll + 2 * K

cat(sprintf(
  "\n[ITER %d] logLik = %.6f | AIC = %.6f | max|score| = %.3e\n",
  iter, ll, AIC_atual, crit
))

cat("Parâmetros atualizados:\n")
print(delta)

if (p > 0) {
  cat("\nCoeficientes AR (phi):\n")
  print(delta[(r + 1):(r + p)])
}

if (q > 0) {
  cat("\nCoeficientes MA (theta):\n")
  print(delta[(r + p + 1):(r + p + q)])
}
