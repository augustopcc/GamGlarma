# Leitura dos dados
library(readr)
BH5 <- read_delim("Fwd_ Dados/BH5.csv", delim = ";", 
                  escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                         grouping_mark = "."), trim_ws = TRUE)

serieBH5 <- ts(BH5[,2],start=1970,frequency=12) #confirmar ano que ocmeçou
n <- length(serieBH5)   # Tamanho da s�rie

# Valores para serem informados
H <- 12   		# Numero de previsoes
nT <- n-H		# Numero de observacoes para o conjunto treinamento
alfa <- 0.05  	# Valor de alfa para o intervalo de previsao
Real <- serieBH5[(n-H+1):n]	# Valores reais para previsao

ts.plot(serieBH5, xlab = NULL,ylab="Numero de Casos de Poliomielite") 


#################################################################
# MODELO MLG

 
X_MLG<-cbind(as.matrix(BH5[, 3:15]))

MLG <- glm(serieBH5[1:nT] ~ X_MLG[1:nT,], family = poisson())
summary(MLG)

library(corrplot)

M <- cor(X_MLG) # Calcula a matriz de correlação de Pearson

# Gerando o gráfico de matriz
corrplot(M, method = "circle", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)



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


X_GLARMA<-cbind(as.matrix(BH5[, 3:15]))

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Modelo usando Residuos de Pearson

# AR(1) - Ajuste do modelo AR(1) com residuos de Pearson (lambda=0,5) 
MP.AR1 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,], phiLags = c(1), thetaLags = NULL, 
                 type = "Poi", residuals = "Pearson")
summary(MP.AR1)


# MA(1) - Ajuste do modelo MA(1) com residuos de Pearson (lambda=0,5) 
MP.MA1 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,], phiLags = NULL, thetaLags = c(1), 
                 type = "Poi", residuals = "Pearson")
summary(MP.MA1)

#  Phi1 e Theta1 sao significativos

#menor aic do AR com 778

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Modelo usando Residuos Escore

## Ajuste do modelo AR(1) com residuos Escore (lambda=1)
MS.AR1 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,], phiLags = c(1), thetaLags = NULL, 
                 type = "Poi", residuals = "Score")
summary(MS.AR1)

MS.MA1 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,], phiLags = NULL, thetaLags = c(1), 
                 type = "Poi", residuals = "Score")
summary(MS.MA1)

# Melhor AR. Theta1 explodiu com Score




## Ajuste do modelo AR(2) com residuos Pearson (lambda=1)
MS.AR1 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,-c(1,2,3,12)], phiLags = c(1), thetaLags = NULL, 
                 type = "Poi", residuals = "Pearson")
summary(MS.AR1)

## Ajuste do modelo AR(2) com residuos Escore (lambda=1) - SEM I3
MS.AR2.2 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,-c(1,2,3,12)], phiLags = c(1,2), thetaLags = NULL, 
                 type = "Poi", residuals = "Pearson")
summary(MS.AR2.2)

## Ajuste do modelo AR(3) com residuos Escore (lambda=1) - SEM I3
MS.AR3 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,-c(1,2,3,12)], phiLags = c(1,2,3), thetaLags = NULL, 
                 type = "Poi", residuals = "Pearson")
summary(MS.AR3)



#ficou sem significancia no 3
##Comparando com MA
## Ajuste do modelo AR(2) com residuos Pearson (lambda=1)
MP.MA1 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,-c(1,2,3,12)], phiLags = NULL, thetaLags = c(1), 
                 type = "Poi", residuals = "Pearson")
summary(MP.MA1)

## Ajuste do modelo AR(2) com residuos Escore (lambda=1) - SEM I3
MP.MA2.2 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,-c(1,2,3,12)], phiLags = NULL, thetaLags = c(1,2), 
                   type = "Poi", residuals = "Pearson")
summary(MP.MA2.2)

#Theta nao singnificativo
## Ajuste do modelo AR(3) com residuos Escore (lambda=1) - SEM I3
MP.MA3 <- glarma(serieBH5[1:nT], X_GLARMA[1:nT,-c(1,2,3,12)], phiLags = NULL, thetaLags = c(1,2,3), 
                 type = "Poi", residuals = "Pearson")
summary(MP.MA3)




# Melhor modelo:
GLARMA <- MS.AR1



#####################################
# Comparando AIC de MLG, AR(2)-Pearson e AR(4)-Escore
cbind(AIC(MLG),GLARMA$aic)


Ajuste_MLG<-ts(fitted(MLG),start=1970,frequency=12)
Ajuste_GLARMA<-ts(fitted(GLARMA),start=1970,frequency=12)
par(mfrow=c(2,1))
plot(serieBH5,type='l',xlab='tempo',ylab='Polio', main="Ajuste MLG")
lines(Ajuste_MLG, col='blue')
plot(serieBH5,type='l',xlab='tempo',ylab='Polio', main="Ajuste GLARMA")
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



