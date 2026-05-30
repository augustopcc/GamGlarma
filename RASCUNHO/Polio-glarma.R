# Leitura dos dados
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

MLG <- glm(Polio[1:nT] ~ X_MLG[1:nT,], family = poisson())
summary(MLG)


# Residuos de Pearson
r_MLG=residuals(MLG,type='pearson')

# Grafico de residuos
par(mfrow=c(2,2))
plot(r_MLG,main="Residuos de Pearson")
acf(r_MLG)
pacf(r_MLG)


# Teste Box-pierce até a ordem 12
Box.test(r_MLG, lag = 12, type = c("Box-Pierce", "Ljung-Box"))
Box.test(r_MLG, lag = 12, type = c("Box-Pierce", "Ljung-Box"))




########################################################################
# GLARMA

install.packages("glarma")
require(glarma)

aux2 <- as.matrix(sPolio[, 3:8])
X_GLARMA<-cbind(aux2,I1,I2,I3)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Modelo usando Residuos de Pearson

  # AR(1) - Ajuste do modelo AR(1) com residuos de Pearson (lambda=0,5) 
  MP.AR1 <- glarma(Polio[1:nT], X_GLARMA[1:nT,], phiLags = c(1), thetaLags = NULL, 
			  type = "Poi", residuals = "Pearson")
  summary(MP.AR1)


  # MA(1) - Ajuste do modelo MA(1) com residuos de Pearson (lambda=0,5) 
  MP.MA1 <- glarma(Polio[1:nT], X_GLARMA[1:nT,], phiLags = NULL, thetaLags = c(1), 
			  type = "Poi", residuals = "Pearson")
  summary(MP.MA1)

# Nem Phi1, nem Theta1 sao significativos


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Modelo usando Residuos Escore

  ## Ajuste do modelo AR(1) com residuos Escore (lambda=1)
  MS.AR1 <- glarma(Polio[1:nT], X_GLARMA[1:nT,], phiLags = c(1), thetaLags = NULL, 
			  type = "Poi", residuals = "Score")
  summary(MS.AR1)

  MS.MA1 <- glarma(Polio[1:nT], X_GLARMA[1:nT,], phiLags = NULL, thetaLags = c(1), 
			  type = "Poi", residuals = "Score")
  summary(MS.MA1)

# Melhor AR. Theta1 nao e significativo


  ## Ajuste do modelo AR(2) com residuos Escore (lambda=1)
  MS.AR2 <- glarma(Polio[1:nT], X_GLARMA[1:nT,], phiLags = c(1,2), thetaLags = NULL, 
			  type = "Poi", residuals = "Score")
  summary(MS.AR2)

  ## Ajuste do modelo AR(2) com residuos Escore (lambda=1) - SEM I3
  MS.AR3 <- glarma(Polio[1:nT], X_GLARMA[1:nT,-9], phiLags = c(1,2), thetaLags = NULL, 
			  type = "Poi", residuals = "Score")
  summary(MS.AR3)

  ## Ajuste do modelo AR(3) com residuos Escore (lambda=1) - SEM I3
  MS.AR3 <- glarma(Polio[1:nT], X_GLARMA[1:nT,-9], phiLags = c(1,2,3), thetaLags = NULL, 
			  type = "Poi", residuals = "Score")
  summary(MS.AR3)

  ## Ajuste do modelo AR(4) com residuos Escore (lambda=1) - SEM I3
  MS.AR4 <- glarma(Polio[1:nT], X_GLARMA[1:nT,-9], phiLags = c(1,2,3,4), thetaLags = NULL, 
			  type = "Poi", residuals = "Score")
  summary(MS.AR4)

  

  ## Ajuste do modelo AR(5) com residuos Escore (lambda=1) - SEM I3
  MS.AR5 <- glarma(Polio[1:nT], X_GLARMA[1:nT,-9], phiLags = c(1,2,3,4,5), thetaLags = NULL, 
			  type = "Poi", residuals = "Score")
  summary(MS.AR5)  # AR(5) nao significativo

# Melhor modelo:
GLARMA <- MS.AR4



#####################################
# Comparando AIC de MLG, AR(2)-Pearson e AR(4)-Escore
cbind(AIC(MLG),GLARMA$aic)


Ajuste_MLG<-ts(fitted(MLG),start=1970,frequency=12)
Ajuste_GLARMA<-ts(fitted(GLARMA),start=1970,frequency=12)
par(mfrow=c(2,1))
plot(Polio,type='l',xlab='tempo',ylab='Polio', main="Ajuste MLG")
lines(Ajuste_MLG, col='blue')
plot(Polio,type='l',xlab='tempo',ylab='Polio', main="Ajuste GLARMA")
lines(Ajuste_GLARMA, col='red')
#leg.txt <- c("Pearson", "Escore")
#legend("topright", leg.txt, pch = "PS", col = c("blue", "red"))

# Analise de residuos
r_GLARMA=residuals(GLARMA,type='Pearson')
par(mfrow=c(2,2))
plot(r_GLARMA,main="Residuos Pearson")
acf(r_GLARMA)
pacf(r_GLARMA)

# Teste Box-pierce até a ordem 1
Box.test(r_GLARMA, lag = 1, type = c("Box-Pierce", "Ljung-Box"))

# Teste Box-pierce até a ordem 12
Box.test(r_GLARMA, lag = 12, type = c("Box-Pierce", "Ljung-Box"))



# BINOMIAL NEGATIVA 
M<- glarma(Polio[1:nT], X_GLARMA[1:nT,-9], phiLags = c(1,2,3,4), thetaLags = NULL, 
			  type = "NegBin", residuals = "Score")
summary(M)



