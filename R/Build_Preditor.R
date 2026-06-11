#' @title Construtor da Matriz de Desenho (Preditor Linear)
#'
#' @description
#' Esta função processa as variáveis explicativas e constrói a matriz de desenho final (\code{X})
#' incorporando efeitos lineares, suavizações não-lineares (Splines), componentes sazonais
#' trigonométricos (Seno e Cosseno), tendência determinística e variáveis de intervenção
#' detetadas automaticamente (Passos e Impulsos).
#'
#' Embora seja o motor interno da função \code{Gam_Glarma}, pode ser utilizada de forma
#' independente para engenharia de variáveis (Feature Engineering) em análises de séries temporais.
#'
#' @param X_in Matriz ou data.frame contendo as covariáveis originais lineares.
#' @param y Vetor numérico da variável resposta (utilizado pelo algoritmo para detetar os pontos exatos de intervenção de Passos e Impulsos).
#' @param n_spline Inteiro. Graus de liberdade (número de nós) para as splines. Se for 0, nenhuma transformação é aplicada.
#' @param spline_cols Vetor de caracteres ou índices indicando quais colunas de \code{X_in} receberão a transformação spline.
#' @param grau_spline Inteiro. Grau do polinômio da spline. Valores >= 3 ativam \code{Natural Cubic Splines} (mais suaves). Valores < 3 ativam \code{B-Splines} (maior variabilidade local).
#' @param sen_cos Vetor numérico indicando os ciclos sazonais a serem criados (ex: \code{c(6, 12)} para semestral e anual).
#' @param tendencia Lógico. Se \code{TRUE}, adiciona uma coluna de tendência linear baseada no tempo.
#' @param impulso Inteiro. Número de choques/picos atípicos (outliers) a serem detetados e incorporados como variáveis dummy no modelo.
#' @param passo Lógico. Se \code{TRUE}, ativa a deteção de quebras de nível (degraus/mudanças de patamar).
#' @param n_passos Inteiro. Número máximo de quebras de nível estruturais a serem procuradas na série.
#' @param tol_passo Numérico. Limite de tolerância para a deteção dos passos.
#'
#' @return Uma lista estruturada contendo:
#' \itemize{
#'   \item \strong{\code{X_out}}: A matriz de desenho completa e expandida, pronta para estimação estatística.
#'   \item \strong{\code{Formula_spline}}: Uma lista com as funções de base históricas das splines geradas (essencial para projetar novos dados na função \code{predict}).
#'   \item \strong{Indicadores}: Lógicos internos (\code{ind_spline}, \code{ind_impulso}, etc.) mapeando quais transformações foram ativadas.
#' }
#'
#' @examples
#' \dontrun{
#' # Exemplo de uso independente para Engenharia de Variáveis
#' X_base <- matrix(rnorm(100), ncol = 1)
#' colnames(X_base) <- "Temperatura"
#' Y_obs <- rpois(100, lambda = 5)
#'
#' matriz_expandida <- Build_Preditor(
#'   X_in = X_base,
#'   y = Y_obs,
#'   n_spline = 3,
#'   spline_cols = "Temperatura",
#'   sen_cos = c(12),
#'   impulso = 1
#' )
#'
#' head(matriz_expandida$X_out)
#' }
#'
#' @export
Build_Preditor <- function(Y_in,
                           X_in,
                           n_spline = FALSE,
                           spline_cols = NULL,
                           sen_cos= c(), # receber c(6,7, 12) -periodicidade
                           tendencia = F,
                           impulso = 0, # receber 1,2 ou 5 - Numero de 1uantos impulsos a pessoa quer
                           passo = F,    # Define se vai criar o passo
                           n_passos = 0,   #permite escolher quantos passos sera identificados
                           tol_passo = 1e-1) {

  X_in <- as.matrix(X_in)
  n <- nrow(X_in)
  Y_in <- as.vector(Y_in)
  X_out <- X_in

  ####---------------- Adiciona os splines ----------------####

  if (n_spline > 0) {

    if (is.character(spline_cols)) {
      spline_cols_idx <- match(spline_cols, colnames(X_in))
      if (any(is.na(spline_cols_idx))) {
        warning("Algumas colunas informadas em spline_cols não foram encontradas em X.")
        spline_cols_idx <- spline_cols_idx[!is.na(spline_cols_idx)]
      }
      spline_cols <- spline_cols_idx
    }

    # Detecção automática (se o usuário não passar as colunas)
    if (is.null(spline_cols) || length(spline_cols) == 0) {
      spline_cols <- integer(0)
      for (j in seq_len(ncol(X_in))) {
        xj <- X_in[, j]
        if (all(is.finite(xj))) {
          uniq <- length(unique(xj))
          is_intercept <-  "(Intercept)" %in% colnames(X_in)
          if (!is_intercept && uniq > (grau_spline + 1)) {
            spline_cols <- c(spline_cols, j)
          }
        }
      }
    }

    if (length(spline_cols) > 0) {
      if (!requireNamespace("splines", quietly = TRUE)) {
        stop("Para spline > 0, o pacote 'splines' precisa estar disponível.")
      }

      B_list <- vector("list", length(spline_cols))
      Formula_spline <- list()
      for (k in seq_along(spline_cols)) {
        j <- spline_cols[k]
        nome_original <- colnames(X_in)[j]
        dados_x <- X_in[, j]

        # Cria a base da spline
        B_j <- splines::ns(dados_x,
                           df = n_spline)
      
        B_j <- scale(B_j, center= TRUE, scale = FALSE)

        colnames(B_j) <- paste0(nome_original, "_s", seq_len(ncol(B_j)))

        # 2. CORREÇÃO DA COLINEARIDADE:
        # A spline SUBSTITUI a variável original.
        B_list[[k]] <- B_j

        #salva a spline para usar a forma de calculo em valores preditos:

        Formula_spline[[nome_original]] <- B_j

      }

      X_spline <- do.call(cbind, B_list)

      # Separa as colunas que NÃO receberam spline (para não perdê-las)
      linear_cols <- setdiff(seq_len(ncol(X_in)), spline_cols)
      X_linear <- if(length(linear_cols) > 0) X_in[, linear_cols, drop = FALSE] else NULL

      # Junta a matriz final (Linear + Splines)
      X_out <- cbind(X_linear, X_spline)
    }
  } else { Formula_spline  <-  NULL}
  ####---------------- Adiciona Tendencia ----------------####

  if (isTRUE(tendencia)) {

    #cria um vetor com pontos que formam uma reta de 0 a 1. O tamanho do vetor é n
    vetor_tendencia <- as.matrix(1:n * (1/n))

    #Renomeia como tendencia e junta a matriz de preditor Linear
    colnames(vetor_tendencia)  <- "Tendencia"

    X_out<- cbind(X_out,vetor_tendencia)

  }

  ####---------------- Adiciona Seno Cosseno ----------------####

  if(!is.null(sen_cos)){

    #cria matriz vazia que recebera os senos e cossenos
    X_sencos <- NULL

    #para cada periodicidade definida no Vetor Sen_Cos é feito o calculo da
    # coluna de seno e cosseno com periodo I
    for (i in sen_cos) {

      sencos_i  <- cbind(
        sin((
          2 * pi * rep(1:i, length.out = n)
        ) / i),
        cos((
          2 * pi * rep(1:i, length.out = n)
        ) / i)
      )

      X_sencos <- cbind(X_sencos,sencos_i)

    }

    #Renomeia as colunas de senos e cossenos e junta a matriz de preditor Linear
    colnames(X_sencos)  <- paste0(c("sen","cos"),rep(x = sen_cos,each = 2))

    X_out <- cbind(X_out,X_sencos)

  }

  ####---------------- Adiciona Passo ----------------####


  if (n_passos > 0) {
    #encontra os pontos em que as medias sao distintas

    #encontra o numero de passos pedidos
    idx <- changepoint::cpt.mean(Y_in,
                                 class = F,
                                 Q = n_passos,
                                 method = "BinSeg")

    n <- length(Y_in)

    pontos <- c(0, idx, n)

    n_seg <- length(pontos) - 2

    # determina a media por segmento
    medias <- sapply(seq_along(pontos[-1]), function(i) {
      mean(Y_in[(pontos[i]+1):pontos[i+1]])
    })

    # identificar níveis únicos (com tolerância)
    niveis <- numeric(0)
    dummies <- list()

    tol= tol_passo

    grupo  <- integer(n_seg)

    #Determina quais sequencias são iguais ou nao, para nao gerar dummys repetidas para o mesmo nivel
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

    # criar dummies
    matriz_passo <- matrix(0, nrow = n, ncol = length(niveis_dummy))
    colnames(matriz_passo) <- paste0("Passo_", seq_along(niveis_dummy))

    for (k in seq_along(niveis_dummy)) {
      nivel_k <- niveis_dummy[k]
      segs_k  <- which(grupo == nivel_k)

      for (s in segs_k) {
        ini <- pontos[s] + 1
        fim <- pontos[s + 1]
        matriz_passo[ini:fim, k] <- 1
      }
    }
    #junta os dummys de passo
    X_out <- cbind(X_out,matriz_passo)

  }

  ####---------------- Adiciona Impulso ----------------####

  if(impulso >=1){

    if (isFALSE(passo)) {


      y_padronizado <- abs(scale(Y_in))

      idx <- order(y_padronizado, decreasing = TRUE)[1:impulso]

      matriz_impulso <- NULL

      for (i in idx) {
        impulso_i <- rep(0,n)
        impulso_i[i] <- 1
        matriz_impulso <- cbind(matriz_impulso,impulso_i)
      }

      #Renomeia as colunas de senos e cossenos e junta a matriz de preditor Linear
      colnames(matriz_impulso)  <- paste0("Impulso_",1:impulso)

      X_out<- cbind(X_out,matriz_impulso)
    }

    if (isTRUE(passo)) {#considera cada trecho

      y_padronizado_local <- numeric(n)

      for (s in seq_len(n_seg)) {
        ini <- pontos[s] + 1
        fim <- pontos[s + 1]

        y_trecho <- Y_in[ini:fim]

        # Trava de segurança: só padroniza se o trecho tiver variação e mais de 1 ponto
        if (length(y_trecho) > 1 && sd(y_trecho) > 0) {
          # Calcula o desvio absoluto em relação à média apenas daquele trecho
          y_padronizado_local[ini:fim] <- abs(scale(y_trecho))
        } else {
          y_padronizado_local[ini:fim] <- 0
        }
      }

      # Encontra os maiores outliers baseados na dinâmica de seus próprios trechos
      idx_outliers <- order(y_padronizado_local, decreasing = TRUE)[1:impulso]

      # Constrói a matriz de dummies para os impulsos
      matriz_impulso <- matrix(0, nrow = n, ncol = impulso)
      for (i in seq_along(idx_outliers)) {
        matriz_impulso[idx_outliers[i], i] <- 1
      }

      colnames(matriz_impulso) <- paste0("Impulso_", 1:impulso)
      X_out <- cbind(X_out, matriz_impulso)
    }

  }


  #### ---------------- Retorna a matriz X ----------------
  return(list(X_out = X_out,
              Formula_spline = Formula_spline))
}
