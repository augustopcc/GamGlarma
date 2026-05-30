#' Ajuste de Modelos GAM-GLARMA
#'
#' @description GamGlarma é usado para ajustar modelos Glarmas com ou sem partes aditivas, fornecendo os preditores lineares, os coeficientes arma e a descrição dos resíduos.
#'
#' @param formula Um objeto da classe \code{formula} (ex: \code{y ~ x1 + x2}).
#' @param data Um \code{data.frame} opcional contendo as variáveis.
#' @param type Distribuição da família exponencial ("Poisson", "Normal", "Gamma").
#' @param link Função de ligação ("log", "identity", "inverse").
#' @param n_spline Número de graus de liberdade para as splines.
#' @param auto_spline Lógico. Se TRUE, busca o melhor grau de spline (1 a 5) baseado no AIC.
#' @param passo Lógico. Se TRUE, utiliza o pacote changepoint para detectar quebras de média.
#'
#' @return Uma lista da classe \code{gamglarma} contendo coeficientes, deviance, AIC e a série latente Z.
#' @export
Gam_Glarma <- function(
    formula,            # Primeiro argumento passa a ser a fórmula
    data = NULL,        # Segundo argumento passa a ser o data.frame
    offset = NULL,
    type = NULL,
    link = NULL,
    method = "FS",
    residuals = "Score",
    phiLags = NULL,
    thetaLags = NULL,
    phiInit = NULL,
    thetaInit = NULL,
    beta = NULL,
    alphaInit = NULL,
    alpha = 1,
    maxit = 30,
    grad = 2.22e-16,
    n_spline = 0,
    auto_spline = F,
    spline_cols = NULL,
    nos_spline = 3,
    grau_spline = 2,
    sen_cos = c(),
    tendencia = F,
    impulso = 0,
    passo = F,
    n_passos = 0,
    tol_passo = 1e-1,
    patience = 3,
    tol_aic_imp = 1e-8,
    trace = FALSE
) {

  # Captura a chamada original
  cl <- match.call()

  #Auto spline, detectar o melhor numero de spline segundo AIC e se algum é significativo. Roda a função gamglarma 5 vezes.
  if (isTRUE(auto_spline)) {
    cat("\n[Auto-Spline] Iniciando busca pelo melhor grau de spline (1 a 5)...\n")

    melhor_aic <- Inf
    melhor_modelo <- NULL
    melhor_nspl <- 0

    for (k in 1:5) {
      # 1. Clona a chamada original e substitui o n_spline
      call_k <- cl
      call_k$n_spline <- k
      call_k$auto_spline <- FALSE # Trava de segurança vital contra loop infinito
      call_k$trace <- FALSE       # Silencia o console para não poluir a tela

      # 2. Roda o modelo usando try() para evitar que falhas de matriz parem o código
      mod_k <- try(eval(call_k, parent.frame()), silent = TRUE)

      # 3. Avalia se o modelo rodou perfeitamente
      if (!inherits(mod_k, "try-error") && mod_k$converged) {

        # 4. Avalia a Significância (Teste de Wald: Z = Beta / Erro Padrão)
        idx_splines <- grep("_s[0-9]+$", names(mod_k$delta))
        significativo <- TRUE # Assume TRUE caso não ache splines (ex: n_spline=0)

        if (length(idx_splines) > 0) {
          z_values <- mod_k$delta[idx_splines] / mod_k$se[idx_splines]
          p_values <- 2 * (1 - pnorm(abs(z_values)))

          # Regra: Pelo menos UMA das bases geradas pela spline deve ser significativa (<0.05)
          significativo <- any(p_values < 0.05)
        }

        # 5. O Veredito: O AIC é o menor E tem relevância estatística?
        if (mod_k$aic < melhor_aic && significativo) {
          melhor_aic <- mod_k$aic
          melhor_modelo <- mod_k
          melhor_nspl <- k
        }

        cat(sprintf("  -> Testando n_spline = %d | AIC: %.2f | Significativo: %s\n",
                    k, mod_k$aic, ifelse(significativo, "Sim", "Não")))
      } else {
        cat(sprintf("  -> Testando n_spline = %d | Falhou em convergir.\n", k))
      }
    }

    if (is.null(melhor_modelo)) {
      stop("\n[Auto-Spline] Falha: Nenhum modelo com splines foi estatisticamente significativo ou convergiu. Recomenda-se ajustar o modelo linear padrão.")
    }

    cat(sprintf("\n[Auto-Spline] Melhor Modelo: n_spline = %d (AIC = %.2f)\n\n", melhor_nspl, melhor_aic))

    # Retorna o modelo vencedor instantaneamente, ignorando o resto do código!
    return(melhor_modelo)
  }

  #############################################################
  # 1) Destrincha a Fórmula e constrói as matrizes base
  if (missing(formula)) stop("O argumento 'formula' é obrigatório (ex: y ~ x1 + x2).")

  # Cria o ambiente de dados avaliando a fórmula
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # Extrai y e X puros da fórmula
  y <- model.response(mf, "numeric")
  X_in <- model.matrix(attr(mf, "terms"), mf)
  y_raw <- mf[[1]]

  #Extrai informação se é serie temporal
  is_ts_y <- inherits(y_raw, "ts")
  if (is_ts_y) {
    ts_start <- start(y_raw)
    ts_end   <- end(y_raw)
    ts_freq  <- frequency(y_raw)
  } else {
    ts_start <- NULL
    ts_end   <- NULL
    ts_freq  <- NULL
  }

  # Trata o offset
  if (is.null(offset)) {
    offset <- model.offset(mf)
    if (is.null(offset)) offset <- rep(0, length(y))
  }

  #############################################################
  # 2) Checagens básicas e padronizações
  #############################################################
  if (!(type %in% c("Poisson","Normal", "Gamma"))) stop("type deve ser 'Poisson', 'Normal' ou 'Gamma'.")
  if (!(method %in% c("FS", "NR"))) stop("method deve ser 'FS' ou 'NR'.")
  if (method == "NR") stop("Nesta versao, somente method = 'FS' esta implementado.")
  if (tolower(residuals) == "fisher") residuals <- "Score"
  if (!(residuals %in% c("Score", "Pearson"))) stop("residuals deve ser 'Score' ou 'Pearson'.")
  if (type == "Gamma" && any(y <= 0)) stop("Para type = 'Gamma', y deve ser estritamente positivo.")

  #############################################################
  # 3) Constrói a Matriz de Preditor passando os dados destrinchados
  #############################################################
  Preditor <- Build_Preditor(
    Y_in = y,              # Alimentado pela fórmula
    X_in = X_in,           # Alimentado pela fórmula
    n_spline = n_spline,
    #auto_spline = auto_spline,
    spline_cols = spline_cols,
    sen_cos = sen_cos,
    tendencia = tendencia,
    impulso = impulso,
    passo = passo,
    n_passos = n_passos,
    tol_passo = tol_passo
  )

  X_fit <- Preditor$X_out

  r <- ncol(X_fit)

  #############################################################
  # 4) Palpites iniciais (Definição de Link e Betas)
  #############################################################
  if (is.null(link)) {
    link <- if (type == "Poisson") "log" else if (type == "Gamma") "log" else "identity"
  }

  if (is.null(beta)) {
    if (type == "Poisson") {
      fit0 <- glm.fit(x = X_fit, y = y, family = poisson(link = link), offset = offset)
      beta <- as.numeric(fit0$coefficients)
    } else if (type == "Gamma") {
      fit0 <- glm.fit(x = X_fit, y = y, family = Gamma(link = link), offset = offset)
      beta <- as.numeric(fit0$coefficients)
    }else if (type == "Normal") {
      fit0 <- glm.fit(x = X_fit, y = y, family = gaussian(link = link), offset = offset)
      beta <- as.numeric(fit0$coefficients)
      beta[is.na(beta)] <- 0
    } else {
      beta <- as.numeric(beta)
      if (length(beta) != r) stop("Comprimento de 'beta' inicial incompatível com o modelo.")
    }
  }

  # phi inicial
  p <- length(phiLags)
  if (is.null(phiInit)) {
    phi <- if (p > 0) rep(0.1, p) else numeric(0)
  } else {
    phi <- as.numeric(phiInit)
    if (length(phi) != p) stop("Comprimento de phiInit incompativel com phiLags.")
  }

  # theta inicial
  q <- length(thetaLags)
  if (is.null(thetaInit)) {
    theta <- if (q > 0) rep(0.1, q) else numeric(0)
  } else {
    theta <- as.numeric(thetaInit)
    if (length(theta) != q) stop("Comprimento de thetaInit incompativel com thetaLags.")
  }

  # delta inicial
  delta <- c(beta, phi, theta)
  names(delta) <- c(
    colnames(X_fit),
    if (p > 0) paste0("phi", phiLags) else NULL,
    if (q > 0) paste0("theta", thetaLags) else NULL
  )

  #############################################################
  # 6) Memoria dos lags
  #############################################################
  mpq <- if (p + q > 0) {
    max(c(if (p > 0) phiLags else 0L,
          if (q > 0) thetaLags else 0L))
  } else {
    0L
  }

  #############################################################
  # 7) Controles de parada
  #############################################################
  tol_beta    <- grad
  tol_score   <- grad

  iter        <- 0
  convergiu   <- FALSE
  no_improve  <- 0

  best_delta  <- delta
  best_aic    <- Inf
  best_ll     <- -Inf
  best_iter   <- 0
  best_pass   <- NULL

  #############################################################
  #  Loop principal
  #############################################################
  while (!convergiu && iter < maxit && no_improve < patience) {

    iter <- iter + 1

    pass <- Iteracao_IRLS(
      delta = delta,
      Xf = X_fit,
      y = y,
      offset = offset,
      type = type,
      link = link,
      residuals = residuals,
      alpha = alpha,
      phiLags = phiLags,
      thetaLags = thetaLags,
      p = p,
      q = q,
      mpq = mpq
    )

    ll    <- pass$ll
    ll_d  <- pass$ll_d
    ll_dd <- pass$ll_dd
    r     <- pass$r

    A <- -ll_dd
    step <- try(solve(A, ll_d), silent = TRUE)

    if (inherits(step, "try-error")) {
      warning("Falha na inversao da informacao de Fisher. Ajuste interrompido.")
      break
    }

    delta_old <- delta
    delta_new <- as.numeric(delta + step)
    names(delta_new) <- names(delta_old)

    K <- length(delta_new)
    aic_new <- -2 * ll + 2 * K

    diff_beta  <- max(abs(delta_new - delta_old))
    diff_score <- max(abs(ll_d))

    if (trace) {
      cat(sprintf(
        "\n[ITER %d] logLik = %.8f | AIC = %.8f | max|score| = %.3e | max|delta_new-delta_old| = %.3e | no_improve = %d\n",
        iter, ll, aic_new, diff_score, diff_beta, no_improve
      ))
      cat("Parâmetros propostos nesta iteração:\n")
      print(delta_new)

      if (p > 0 || q > 0) {
        cat("\n--- Parte temporal proposta nesta iteracao ---\n")
        if (p > 0) {
          cat("Coeficientes AR (phi):\n")
          print(delta_new[(r + 1):(r + p)])
        }
        if (q > 0) {
          cat("Coeficientes MA (theta):\n")
          print(delta_new[(r + p + 1):(r + p + q)])
        }
        cat("----------------------------------------------\n")
      }
    }

    # guarda melhor AIC
    if (aic_new < best_aic - tol_aic_imp) {
      best_aic   <- aic_new
      best_ll    <- ll
      best_delta <- delta_new
      best_iter  <- iter
      best_pass  <- pass
      no_improve <- 0
    } else {
      no_improve <- no_improve + 1
    }

    delta <- delta_new

    # convergência numérica
    if (diff_beta < tol_beta || diff_score < tol_score) {
      convergiu <- TRUE
      if (trace) cat("\n*** Critério de convergência numérica atingido ***\n")
    }

    if (no_improve >= patience && trace) {
      cat("\n*** Patience atingido: sem melhora do AIC por ", patience, " iterações consecutivas ***\n", sep = "")
    }
  }

  #############################################################
  # 9) Recupera o melhor delta e recalcula o ajuste final
  #############################################################
  delta <- best_delta
  names(delta) <- c(
    colnames(X_fit),
    if (p > 0) paste0("phi", phiLags) else NULL,
    if (q > 0) paste0("theta", thetaLags) else NULL
  )

  # --- PASSO 9.1: Estimar o Parâmetro de Dispersão Real ---
  # Roda um passe rápido apenas para recuperar o vetor de médias (mu)
  temp_pass <- Iteracao_IRLS(
    delta = delta, Xf = X_fit, y = y, offset = offset,
    type = type, link = link, residuals = residuals,
    alpha = alpha, phiLags = phiLags, thetaLags = thetaLags,
    p = p, q = q, mpq = mpq
  )
  fitted_mu_temp <- temp_pass$mu[(mpq + 1):(mpq + length(y))]

  # Calcula o alpha (shape) correto baseado no resíduo de Pearson
  if (type == "Gamma") {
    df_residual <- length(y) - length(delta)
    phi_dispersion <- sum( ((y - fitted_mu_temp) / fitted_mu_temp)^2 ) / df_residual
    alpha_final <- 1 / phi_dispersion
  } else if (type == "Normal") {
    # Para a Normal, a dispersão (variância) é a média dos erros quadráticos
    # alpha_final aqui funcionará como o Inverso da Variância (1/sigma^2)
    variancia_est <- sum((y - fitted_mu_temp)^2) / (length(y) - length(delta))
    alpha_final <- 1 / variancia_est
  } else {
    alpha_final <- 1 # Poisson
  }

  # --- PASSO 9.2: Recalcular a Matriz de Fisher e Verossimilhança com alpha correto ---
  final_pass <- Iteracao_IRLS(
    delta = delta, Xf = X_fit, y = y, offset = offset,
    type = type, link = link, residuals = residuals,
    alpha = alpha_final, phiLags = phiLags, thetaLags = thetaLags,
    p = p, q = q, mpq = mpq
  )

  ll    <- final_pass$ll
  ll_d  <- final_pass$ll_d
  ll_dd <- final_pass$ll_dd

  K <- length(delta)

  # O AIC do R adiciona 2 unidades pelo parâmetro de dispersão estimado na Gamma
  AIC_final <- -2 * ll + 2 * (K + (if(type == "Gamma") 1 else 0))

  fitted_mu <- final_pass$mu[(mpq + 1):(mpq + length(y))]
  resids    <- final_pass$e[(mpq + 1):(mpq + length(y))]
  state_Z   <- final_pass$Z[(mpq + 1):(mpq + length(y))]
  eta_final <- final_pass$W[(mpq + 1):(mpq + length(y))]

  # --- PASSO 9.3: Cálculo Exato da Deviance Residual ---
  if (type == "Poisson") {
    res_dev <- 2 * sum(y * log(ifelse(y == 0, 1, y / fitted_mu)) - (y - fitted_mu))
  } else if (type == "Gamma") {
    res_dev <- 2 * sum(-log(y / fitted_mu) + (y - fitted_mu) / fitted_mu)
  } else if (type == "Normal") {
    # Deviance da normal é a Soma dos Quadrados dos Resíduos (RSS)
    res_dev <- sum((y - fitted_mu)^2)
  }

  vcov_FS <- try(solve(-ll_dd), silent = TRUE)

  if (!inherits(vcov_FS, "try-error")) {
    se_FS <- sqrt(diag(vcov_FS))
  } else {
    vcov_FS <- NULL
    se_FS <- rep(NA_real_, length(delta))
  }

  is_glm_only <- (length(phiLags) == 0 && length(thetaLags) == 0)

  has_intercept <- "(Intercept)" %in% colnames(X_in)

  if (has_intercept) {
    fit_null <- try(glm(y ~ 1, family = if(type=="Poisson") poisson(link=link) else if(type=="Gamma") Gamma(link=link) else gaussian(link=link), offset = offset), silent = TRUE)
  } else {
    fit_null <- try(glm(y ~ 0, family = if(type=="Poisson") poisson(link=link) else if(type=="Gamma") Gamma(link=link) else gaussian(link=link), offset = offset), silent = TRUE)
  }

  # 3. Salva os valores para a lista de saída
  if (!inherits(fit_null, "try-error")) {
    null_loglik_val <- as.numeric(logLik(fit_null))
    null_deviance_val <- fit_null$deviance
  } else {
    null_loglik_val <- NA
    null_deviance_val <- NA
  }
  #############################################################
  # 10) Saída no estilo objeto de modelo
  #############################################################
  out <- list(
    call = cl,
    delta = delta,
    logLik = ll,
    logLikDeriv = ll_d,
    logLikDeriv2 = ll_dd,
    eta = eta_final,
    W = final_pass$W,
    Z = final_pass$Z,
    mu = fitted_mu,
    fitted.values = fitted_mu,
    residuals = resids,
    cov = vcov_FS,
    se = se_FS,
    phiLags = phiLags,
    thetaLags = thetaLags,
    r = r,
    pq = p + q,
    y = y,
    X = X_fit,
    offset = offset,
    type = type,
    link= link,
    method = method,
    residType = residuals,
    iter = iter,
    bestIter = best_iter,
    bestAIC = best_aic,
    bestLogLik = best_ll,
    converged = convergiu || (no_improve >= patience),
    patience = patience,
    minScore = max(abs(ll_d)),
    aic = AIC_final,

    spline = spline,
    spline_cols = spline_cols,

    X_original = X_in,
    ind_spline = n_spline > 0,
    ind_tendencia= tendencia,
    ind_sen_cos = length(sen_cos) > 0,
    ind_impulso = impulso > 0 ,
    ind_passo = passo,
    sen_cos = sen_cos,                               # O vetor numérico de ondas (ex: c(6, 12))
    Formula_spline = Preditor$Formula_spline,        # As réguas matemáticas (ns)
    variaveis_lineares = setdiff(colnames(X_in), names(Preditor$Formula_spline)),

    alpha = alpha_final,            # Retorna o alpha atualizado
    nullLogLik = null_loglik_val,
    nullDev = null_deviance_val,
    resDev = res_dev,               # Passa a Deviance real para o Summary
    is_glm = is_glm_only,

    is_ts_y = is_ts_y,
    ts_start = ts_start,
    ts_end = ts_end,
    ts_freq = ts_freq
  )

  class(out) <- "gamglarma"
  return(out)
}





