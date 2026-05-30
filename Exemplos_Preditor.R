#

x=c(rnorm(50,0,1),rnorm(50,2,1),rnorm(50,0,1),rnorm(50,2,1),rnorm(50,4,1))


## Exemplos de matriz preditora ####

##### caso com impulso #####

x_impulso = cbind(rpois(n = 150,lambda = 5),rgamma(150,shape = 2,rate = 1))
y_impulso= rnorm(150,0,1) + c(rep(0,80),5,rep(0,69))
plot(ts(y_impulso))

preditor <- Build_Preditor(X_in = x_impulso,Y_in = y_impulso, impulso = 1)
View(preditor)

x_impulso = cbind(rpois(n = 150,lambda = 5),rgamma(150,shape = 2,rate = 1))
y_impulso= rnorm(150,0,1) + c(rep(0,40),6,rep(0,39),-7,rep(0,69))
plot(ts(y_impulso))

preditor <- Build_Preditor(X_in = x_impulso,Y_in = y_impulso, impulso = 2)
View(preditor)

##### caso com tendencia #####

x_tendencia = cbind(rpois(n = 150,lambda = 5),rgamma(150,shape = 2,rate = 1))
y_tendencia= rnorm(150,0,1)+ seq(from = 0, to = 4, length.out = 150)
plot(ts(y_tendencia))

preditor <- Build_Preditor(X_in = x_tendencia,Y_in = y_tendencia,tendencia = T)
View(preditor)


##### caso seno cossenos #####

x_sencos = cbind(rpois(n = 150,lambda = 5),rgamma(150,shape = 2,rate = 1))
y_sencos= rnorm(50,0,1) + rep(c(0,0,0,0,2,3,1,0,0,0,1,1), length.out = 150)
plot(ts(y_sencos))

preditor <- Build_Preditor(X_in = x_sencos,Y_in = y_sencos,sen_cos = c(6,12))
View(preditor)
plot(ts(preditor[,3]))
plot(ts(preditor[,4]))

##### caso com passo #####
#Niveis diferentes
x_passo = cbind(rpois(n = 150,lambda = 5),rgamma(150,shape = 2,rate = 1))
y_passo=c(rnorm(50,0,1),rnorm(50,1,1),rnorm(50,3,1))
plot(ts(y_passo))
preditor <- Build_Preditor(X_in = x_passo,Y_in = y_passo,passo = T,tol_passo = .5)
View(preditor)

#Niveis iguais
x_passo = cbind(rpois(n = 250,lambda = 5),rgamma(250,shape = 2,rate = 1))
y_passo=c(rnorm(50,0,1),rnorm(50,2,1),rnorm(50,0,1),rnorm(50,2,1),rnorm(50,4,1))
plot(ts(y_passo))
preditor <- Build_Preditor(X_in = x_passo,Y_in = y_passo,passo = T,tol_passo = .5)
View(preditor)

##### Caso com Spline #####
n <- 150
x <- 1:n
x_spline = cbind(coluna1= rpois(n = 150,lambda = 5),
                 coluina2=ifelse(x <= 50, 
                        0.2 * x,                    # Sobe até o 50
                        ifelse(x <= 100, 
                               10 - 0.15 * (x - 50), # Desce levemente entre 51 e 100
                               0 + 0.3 * (x - 100) # Sobe forte após o 100
                        )) + rnorm(n, mean = 0, sd = 3))
tendencia_real <- x_spline[,2] 
# Forçando valores positivos para a média da Gamma
mu <- exp(tendencia_real / 10) 
y_spline <- rgamma(n, shape = 2, scale = mu/2)
plot(x_spline[,2])
plot(y_spline,x_spline[,2])

preditor <- Build_Preditor(X_in = x_spline,Y_in = y_spline,spline = T,spline_cols = c(2),nos_spline = 2)
View(preditor)
plot(y_spline,preditor[,2])



n <- 150
x <- seq(1, 10, length.out = n)
# Criando um sinal com várias curvas para testar a flexibilidade
y_sinal <- sin(x) + 0.1 * x^2 + rnorm(n, sd = 0.5)

# Criar a matriz de entrada (X_in)
X_exemplo <- matrix(x, ncol = 1)
colnames(X_exemplo) <- "Tempo"

# 2. Rodar a função com diferentes configurações
# Teste A: Grau 1 (Linear por partes) com muitos nós (deve travar em 4 nós para dar 5 colunas)
pred_linear <- Build_Preditor(Y_in = y_sinal, X_in = X_exemplo, 
                              spline = TRUE, spline_cols = 1, 
                              nos_spline = 4, grau_spline = 1)

# Teste B: Grau 3 (Cúbico) com muitos nós (deve travar em 2 nós para dar 5 colunas)
pred_cubico <- Build_Preditor(Y_in = y_sinal, X_in = X_exemplo, 
                              spline = TRUE, spline_cols = 1, 
                              nos_spline = 2, grau_spline = 3)

# 3. Visualização Comparativa
par(mfrow = c(2, 2))

# Gráfico 1: Bases do Spline Linear (Grau 1)
# Colunas 2 a 6 são os splines (já que a 1 é a original)
matplot(x, pred_linear[, 2:6], type = "l", lty = 1, lwd = 2,
        main = "Bases: Grau 1 (Máx 5 col)", ylab = "Valor da Base")
abline(v = quantile(x, probs = seq(0, 1, length.out = 6)[-c(1, 6)]), col = "red", lty = 3)

# Gráfico 2: Bases do Spline Cúbico (Grau 3)
matplot(x, pred_cubico[, 2:4], type = "l", lty = 1, lwd = 2,
        main = "Bases: Grau 3 (Máx 5 col)", ylab = "Valor da Base")

# Gráfico 3: Ajuste Linear (Como o modelo 'enxerga' os dados)
mod_lin <- lm(y_sinal ~ pred_linear)
plot(x, y_sinal, col = "gray", pch = 19, main = "Ajuste Final (Grau 1)")
lines(x, fitted(mod_lin), col = "blue", lwd = 2)

# Gráfico 4: Ajuste Cúbico
mod_cub <- lm(y_sinal ~ pred_cubico)
plot(x, y_sinal, col = "gray", pch = 19, main = "Ajuste Final (Grau 3)")
lines(x, fitted(mod_cub), col = "red", lwd = 2)

par(mfrow = c(1, 1))

##### Impulso e Passo ##### 

#Niveis diferentes
x_passo = cbind(rpois(n = 150,lambda = 5),rgamma(150,shape = 2,rate = 1))
y_passo=c(rnorm(50,0,1),rnorm(50,1,1),rnorm(50,3,1)) + c(rep(0,80),5,rep(0,69))
plot(ts(y_passo))
preditor <- Build_Preditor(X_in = x_passo,Y_in = y_passo,passo = T,tol_passo = .5,impulso = 1)
View(preditor)

###### ##########
