#### Exemplo de uso do GamGlarma com poisson ####

sPolio <- read.table('Base_Polio 1.txt',header = TRUE)
Polio <- ts(sPolio[,2],start=1970,frequency=12) 
n <- length(Polio)   # Tamanho da s�rie

# Valores para serem informados
H <- 12   		# Numero de previsoes
nT <- n-H		# Numero de observacoes para o conjunto treinamento
alfa <- 0.05  	# Valor de alfa para o intervalo de previsao
Real <- Polio[(n-H+1):n]	# Valores reais para previsao


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
X_GLARMA<-cbind(aux,I1,I2,I3)

gamglarma <- fit_gam <- Gam_Glarma(
  y = sPolio$Polio[1:nT],
  X = X_GLARMA[1:nT, ],
  offset = rep(0, nT),
  type = "Poisson",
  method = "FS",
  residuals = "Score",
  phiLags = c(1L, 2L, 3L, 4L),
  thetaLags = integer(0),
  spline = FALSE,
  #spline_cols = 2:6,   # apenas as contínuas recebem spline
  #df_spline = 4,
  maxit = 100,
  grad = 1e-6,
  patience = 3,
  trace = F
)
gamglarma
summary(gamglarma)
##Compara com pacote glarma:
MS.AR4 <- glarma::glarma(Polio[1:nT], X_GLARMA[1:nT,-9], phiLags = c(1,2,3,4), thetaLags = NULL, 
                         type = "Poi", residuals = "Score")
summary(MS.AR4)
MS.AR4
#valores bateram

############ -------------------- ############
#### Exemplo Gam Glarma BH5 ####
############ -------------------- ############
ambiente <- setdiff (ls(),
                     c("Build_Preditor","Gam_Glarma","Iteracao_IRLS","parte_distribuicao","print.gamglarma","print.summary.gamglarma","summary.gamglarma"))
rm(list = ambiente)


BH5 <- readr::read_delim("Fwd_ Dados/BH5.csv", delim = ";", 
                  escape_double = FALSE, locale = readr::locale(decimal_mark = ",", 
                                                         grouping_mark = "."), trim_ws = TRUE)


library(corrplot)

par(mfrow=c(1,1))
M <- cor(BH5[2:15]) # Calcula a matriz de correlação de Pearson
n <- length(BH5$PM10)   # Tamanho da s�rie

# Valores para serem informados
# H <- 12   		# Numero de previsoes
	
# Gerando o gráfico de matriz
corrplot(M, method = "circle", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

y = BH5$PM10
X = cbind(as.matrix(BH5[,c("RH","NOX")]))

plot(ts(y))
hist(y)
pacf(y)

## Verificando o MLG ##
MLG <- glm(PM10 ~RH  + NO,data = BH5, family = Gamma(link = "log"))
summary(MLG)
MLG
ts(MLG$residuals) |> acf()
ts(MLG$residuals) |> pacf()


gamglarma <- Gam_Glarma(
  ts(PM10,start = 2010, frequency = 12) ~ RH  + NO,
  data = BH5,
  impulso = 1,
  sen_cos = c(6),
  offset = rep(0, n),
  type = "Gamma",link = "log",
  method = "FS",
  residuals = "Score",
  phiLags = NULL,
  thetaLags = 1,
  n_spline = 2, auto_spline = T,spline_cols = "RH",
  maxit = 50,
  grad = 1e-6,
  patience = 3,
  trace = F
)
gamglarma
summary(gamglarma)

plotAjuste(gamglarma)
 plotResiduos(gamglarma)

