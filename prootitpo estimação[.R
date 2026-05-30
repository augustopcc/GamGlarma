



#Esboço do entendimento do codigo:
#Repetir a estimação do Glarma 
##Usamos Verossimilhança, o calculo de residuos tem 2 tipos
##Person com lambda igual a .5 e score com lambda igual a 1
##
##Para inicializar o método recursivo de Newton-Raphson na maximização numérica da log 
##verossimilhança 𝑙 𝜑,𝑦 , Davis et al. (2003) sugerem que os valores obtidos das estimativas do  
##GLARMA sem os termos auto-regressivos média móveis sejam utilizados como valores iniciais. A 
##convergência, na maioria dos casos, ocorre após 10 iterações





#Base para teste
## Leitura dos dados
sPolio <- read.table('Base_Polio 1.txt',header = TRUE)

Polio <- ts(sPolio[,2],start=1970,frequency=12) 
n <- length(Polio)   # Tamanho da s�rie

# Valores para serem informados
H <- 12   		# Numero de previsoes
nT <- n-H		# Numero de observacoes para o conjunto treinamento
alfa <- 0.05  	# Valor de alfa para o intervalo de previsao
Real <- Polio[(n-H+1):n]	# Valores reais para previsao

ts.plot(Polio, xlab = NULL,ylab="Numero de Casos de Poliomielite") 

#################################################################
# MODELO MLG

# Calculo da posicao do primeiro outlier
max<-Polio[1]
for(i in 2:n){
  if(Polio[i] > max){
    max<-Polio[i]
    ind1=i}
}
Polio[ind1]

# Funcao Impulso1
I1 <- rep(0,n)
I1[ind1]=1

# Calculo da posicao do segundo outlier
max<-Polio[1]
for(i in (1):(ind1-1)){
  if(Polio[i] > max){
    max<-Polio[i]
    ind2=i}
}
Polio[ind2]

# Funcao Impulso2
I2 <- rep(0,n)
I2[ind2]=1

# Calculo da posicao do terceiro outlier
max<-Polio[ind1+1]
for(i in (ind1+1):(n)){
  if(Polio[i] > max){
    max<-Polio[i]
    ind3=i}
}
Polio[ind3]

# Funcao Impulso3
I3 <- rep(0,n)
I3[ind3]=1

aux <- as.matrix(sPolio[, 4:8])
X_MLG<-cbind(aux,I1,I2,I3)
## Criação da Base com os indicadores

#A criação dos indicadores de outlier não é algo que estara imbutido no pacote



#Como é feito o ajuste:
#O chute inicial vem de um GLM, este não sera implementado no pacote, sera usado a propria função antivad do R, já otimizada


## =============================================
## Prep. treino e objetos compatíveis com nT
## =============================================
y  <- sPolio$Polio
yT <- y[1:nT]



#### -------------------------------------------------------------- ####
## X com INTERCEPTO para casar com os Betas do GLM (Aqui serao inseridos os SPLINES)

splines <-  T

if (splines == F) {

  #Sem Splines
  X  <- cbind(Intercept = 1, X_MLG[1:nT, ])
  r  <- ncol(X)
    
} else {
  #Com Splines
  
  X_cont <- as.matrix(sPolio[, 4:8])   # covariáveis contínuas
  X_imp  <- cbind(I1, I2, I3)          # impulsos (sem spline)
  
  library(splines)
  
  df_spline <- 2  # graus de liberdade por covariável
  
  B_list <- lapply(seq_len(ncol(X_cont)), function(j) {
    bs(X_cont[, j], df = df_spline)
  })
  
  B <- do.call(cbind, B_list)
  
  X <- cbind(
    Intercept = 1,
    X_imp[1:nT, ],      # impulsos lineares
    B[1:nT, ]           # splines automáticos
  )  
  }




#### ---------------------------------------------------- ####

## offset numérico (evitar conflito com a função stats::offset)
offset_vec <- rep(0, nT)

## Ajuste GLM com o MESMO design de treino (opcionalmente com offset)
MLG <- glm(yT ~ X[1:nT, ], family = poisson(link = "log"), offset = offset_vec)
Betas_iniciais <- coef(MLG)

## GLARMA: AR(4), MA(0)
phiLags   <- c(1L,2L,3L,4L)
thetaLags <- integer(0)

p <- length(phiLags)
q <- length(thetaLags)

## Vetor de parâmetros: delta = (beta, phi, theta)
beta  <- as.numeric(Betas_iniciais)      # já inclui intercepto
phi   <- if (p > 0) rep(0.1, p) else numeric(0)
theta <- if (q > 0) rep(0.1, q) else numeric(0)

delta <- c(beta, phi, theta)
names(delta) <- c(colnames(X),
                  if (p>0) paste0("phi", phiLags) else NULL,
                  if (q>0) paste0("theta", thetaLags) else NULL)

## Largura de memória para lags
mpq <- if (p + q > 0) max(c(if (p>0) phiLags else 0L,
                            if (q>0) thetaLags else 0L)) else 0L

## Objetos do algoritmo com nT + mpq
Z  <- numeric(nT + mpq)
e  <- numeric(nT + mpq)
W  <- numeric(nT + mpq)
mu <- numeric(nT + mpq)

s   <- length(delta)
Z.d <- matrix(0, nrow = s, ncol = nT + mpq)
W.d <- matrix(0, nrow = s, ncol = nT + mpq)
e.d <- matrix(0, nrow = s, ncol = nT + mpq)

tol <- 1e-6



## ============================================================
## 4) UMA ITERAÇÃO de Fisher Scoring (rode este bloco 1x por passo)
##    - Reexecutar ESTE trecho = 1 nova iteração
## ============================================================

# (re)particiona delta atual
r <- ncol(X)
beta  <- delta[1:r]
idx   <- r
phi   <- if (p > 0) delta[(idx+1):(idx+p)] else numeric(0); idx <- idx + p
theta <- if (q > 0) delta[(idx+1):(idx+q)] else numeric(0)

## zera acumuladores desta iteração
Z[] <- 0; e[] <- 0; W[] <- 0; mu[] <- 0
Z.d[,] <- 0; W.d[,] <- 0; e.d[,] <- 0

ll    <- 0
ll.d  <- matrix(0, nrow = s, ncol = 1)
ll.dd <- matrix(0, nrow = s, ncol = s)
## -------------------------------------------------------------------------------- ##
## ---------------- laço temporal no TREINO: t = 1..nT ---------------
## -------------------------------------------------------------------------------- ##
for (t in 1:nT) {
  tt <- t + mpq
  
  ## d eta / d beta começa em X[t,]
  W.d[, tt] <- 0
  W.d[1:r, tt] <- X[t, ]
  
  ## AR(4): contribuições em Z e em derivadas
  if (p > 0) {
    Z.d[(r+1):(r+p), tt] <- Z[tt - phiLags] + e[tt - phiLags]
    for (i in seq_len(p)) {
      lag_i    <- phiLags[i]
      Z[tt]     <- Z[tt] + phi[i] * (Z[tt - lag_i] + e[tt - lag_i])
      Z.d[, tt] <- Z.d[, tt] + phi[i] * (Z.d[, tt - lag_i] + e.d[, tt - lag_i])
    }
  }
  
  ## MA(0):
  
  if (q > 0) {
    Z.d[(r+p+1):(r+p+q), tt] <- e[tt - thetaLags]
    for (j in seq_len(q)) {
      lag_j <- thetaLags[j]
      Z[tt] <- Z[tt] + theta[j] * e[tt - lag_j]
      Z.d[,tt] <- Z.d[,tt] + theta[j] * e.d[,tt - lag_j]
    }
  }
  ## -------------------------------------------------------------------------------- ##
  ## ----------- parte dependente de distribuição -----------------
  ## -------------------------------------------------------------------------------- ##
  
  ## eta, mu e resíduo score (Poisson)
  eta_t  <- as.numeric( drop(X[t, ] %*% beta) + Z[tt] + offset_vec[t] )
  W[tt]  <- eta_t
  
  ####Ligação exponencial
  mu[tt] <- exp(eta_t)
  ####Residuo score
  e[tt]  <- (yT[t] - mu[tt]) / mu[tt]
  
  ## derivadas finais
  W.d[, tt] <- W.d[, tt] + Z.d[, tt]
  e.d[, tt] <- -(1 + e[tt]) * W.d[, tt]
  
  ## loglik, escore e info esperada, mudar a verossimilhança de acordo com a distribuição
  ll    <- ll + (yT[t] * eta_t - mu[tt] - lfactorial(yT[t]))
  ll.d  <- ll.d + (yT[t] - mu[tt]) * matrix(W.d[, tt], ncol = 1)
  ll.dd <- ll.dd - mu[tt] * (W.d[, tt, drop = FALSE] %*% t(W.d[, tt, drop = FALSE]))
}

## passo de Fisher Scoring
A    <- -ll.dd
step <- try(solve(A, ll.d), silent = TRUE)
if (inherits(step, "try-error")) stop("Info(Fisher) não inversível; ajuste palpites, p/q ou regularize.")

delta_old <- delta
delta     <- as.numeric(delta + step)
names(delta) <- names(delta_old)

## Relatório
crit <- max(abs(ll.d))
cat(sprintf("\n[PASSO] logLik=%.6f | max|score|=%.3e\n", ll, crit))
cat("Parâmetros atualizados:\n"); print(delta)



########### esboço #################
#
Teste <- 1:36
Teste2 <- rep(1:6,length.out = length(Teste))
teste3 <- sin((2*pi*Teste2)/6)


### Inserire senos e cossenos

periodicidade <- 12

if (sen_cos == TRUE) {
  
  X_sincos <- cbind(
  sin((
    2 * pi * rep(1:periodicidade, length.out = length(nrow(X_in)))
  ) / periodicidade),
  cos((
    2 * pi * rep(1:periodicidade, length.out = length(nrow(X_in)))
  ) / periodicidade)
  )
  
}


### Insere o impulso

n_impulsos = 2
sPolio <- sPolio$Polio

y <- sPolio

y_padronizado <- abs(scale(y))
idx <- order(y_padronizado, decreasing = TRUE)[1:n_impulsos]
matriz_impulso <- NULL

for (i in idx) {
  impulso_i <- rep(0,n)
  impulso_i[i] <- 1
  matriz_impulso <- cbind(matriz_impulso,impulso_i)
}
#Renomeia as colunas de senos e cossenos e junta a matriz de preditor Linear
colnames(matriz_impulso)  <- paste0("Impulso_",1:n_impulsos)

X_in<- cbind(X_in,matriz_impulso)


### Insere o passo
### 
if (isTRUE(Passo)) {
  
  changepoint::cpt.mean(x,class = T,param.estimates = T)
  
  
}
require("changepoint")

x=c(rnorm(50,0,1),rnorm(50,2,1),rnorm(50,0,1),rnorm(50,2,1),rnorm(50,4,1))
idx <- changepoint::cpt.mean(x,class = F)
n=length(x)

matriz_impulso <- NULL

for (i in idx) {
  impulso_i <- rep(0,n)
  impulso_i[i] <- 1
  matriz_impulso <- cbind(matriz_impulso,impulso_i)
}

#Renomeia as colunas de senos e cossenos e junta a matriz de preditor Linear
colnames(matriz_impulso)  <- paste0("Impulso_",1:n_impulsos)

X_in<- cbind(X_in,matriz_impulso)




y,
                                    cpt_obj, 
                                    tol = 1e-6) {
  
  
  x=c(rnorm(50,0,1),rnorm(50,2,1),rnorm(50,0,1),rnorm(50,2,1),rnorm(50,4,1))
  idx <- changepoint::cpt.mean(x,class = F)

  n <- length(x)
  pontos <- c(0, idx)
  n_seg <- length(pontos) - 1
  
  # médias por segmento
  medias <- sapply(seq_along(pontos[-1]), function(i) {
    mean(x[(pontos[i]+1):pontos[i+1]])
  })
  
  # identificar níveis únicos (com tolerância)
  niveis <- numeric(0)
  dummies <- list()
  
  tol= .5
  
  grupo  <- integer(n_seg)
  
  #
  for (i in seq_len(n_seg)) {
    if (length(niveis) == 0) {
      niveis <- medias[i]
      grupo[i] <- 1
    } else {
      diffs <- abs(niveis - medias[i])
      if (any(diffs <= tol)) {
        grupo[i] <- which.min(diffs)
      } else {
        niveis <- c(niveis, medias[i])
        grupo[i] <- length(niveis)
      }
    }
  }
  
  # nivel base = primeiro nivel (nao gera dummy)
  niveis_dummy <- 2:length(niveis)
  
  if (length(niveis_dummy) == 0) {
    return(list(
      medias_segmentos = medias,
      niveis = niveis,
      dummies = NULL
    ))
  }
  
  # criar dummies
  dummies <- matrix(0, nrow = n, ncol = length(niveis_dummy))
  colnames(dummies) <- paste0("DUMMY_CENTROIDE_", seq_along(niveis_dummy))
  
  for (k in seq_along(niveis_dummy)) {
    nivel_k <- niveis_dummy[k]
    segs_k  <- which(grupo == nivel_k)
    
    for (s in segs_k) {
      ini <- pontos[s] + 1
      fim <- pontos[s + 1]
      dummies[ini:fim, k] <- 1
    }
  }
  
  

  
  list(
    medias_segmentos = medias,
    niveis_detectados = niveis,
    dummies = dummies_mat
  )
}

#tendencia
base_bh5 <- readr::read_delim("Fwd_ Dados/BH5.csv", delim = ";", 
           escape_double = FALSE, col_types = readr::cols(trend = readr::col_number()), 
           locale = readr::locale(decimal_mark = ",", grouping_mark = "."), 
           trim_ws = TRUE)

plot(base_bh5$trend)

n <- length(base_bh5$trend)
tendencia <- 1:n*(1/n)


vetor_tendencia <- as.matrix(1:n * (1/n))

#Renomeia as colunas de senos e cossenos e junta a matriz de preditor Linear
colnames(vetor_tendencia)  <- "Tendencia"

X_in<- cbind(X_in,vetor_tendencia)

########

Y_in <- c(rnorm(50,0,1),rnorm(50,2,1),rnorm(50,3,1))
plot(ts(Y_in))
idx <- changepoint::cpt.mean(Y_in,class = F,Q = 1,method = "BinSeg")
idx
