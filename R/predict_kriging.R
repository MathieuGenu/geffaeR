#' @export


# fonction pour prédire
predict_kriging <- function(segdata_obs, esw, vario_model, distmat_x, distmat_xy,
                            fast_inversion = TRUE, saturate = TRUE) {
  col_needed <- c("y","n","Effort")
  if(!all(col_needed %in% colnames(segdata_obs))) {
    missing_col <- col_needed[!(col_needed %in% colnames(segdata_obs))]
    stop(paste("Les colonnes suivantes ne sont pas dans les colonnes de segdata_obs :",
               paste(missing_col, collapse = "\n"), sep="\n"))
  }
  y <- segdata_obs$y
  l <- segdata_obs$Effort

  # distmat_x = matrice des distances entre points du jeu de données
  # distmat_xy = matrice des distances entre points du jeu de données et points à prédire

  cov_matern <- function(h, sill, range) {
    return(sill *(1 + h * sqrt(3) / range) * exp(- h * sqrt(3) / range))
  }

  cov_exp <- function(h, sill, range) {
    return(sill * exp(- h / range))
  }

  # calcul du taux d'obs moyen m
  m <- TauxObs(y, esw, l)

  # on prend le modele de covariance
  form <- vario_model$cov.model

  # on prend l'ajustement 2 avec parametres:
  sill <- as.numeric(vario_model$cov.pars[1])
  range <- as.numeric(vario_model$cov.pars[2])

  # nombre de points(segments)
  nbpts <- length(y)

  # krigeage

  if(form == "matern") {
    # matrice des covariances entre toutes les données
    a <- cov_matern(distmat_x, sill, range)

    # matrice des covariances entre données et prédictions
    cmat <- cov_matern(distmat_xy, sill = sill, range = range)

  }

  if(form == "exponential") {
    # matrice des covariances entre toutes les données
    a <- cov_exp(distmat_x, sill, range)

    # matrice des covariances entre données et prédictions
    cmat <- cov_exp(distmat_xy, sill = sill, range = range)

  }

  if(fast_inversion) {
    diag(a) <- diag(a) + m /(l * 2 * esw)

    A1 <- chol2inv(chol(a))
    B1 <- matrix(rep(1, nrow(a)), ncol = 1)
    B2 <- A1 %*% B1
    A4 <- -1 / t(B2) %*% B1
    A2 <- B2 %*% A4
    A1 <- A1 + A2 %*% t(B2)

    inv_a <- rbind(cbind(A1, -A2), cbind(-t(A2), A4))
  } else {
    # construire la matrice de covariance + une colonne et une ligne de 1
    a <- rbind(cbind(a, rep(1, nbpts)), c(rep(1, nbpts), 0))

    #ajouter la constante sur la diagonale
    diag(a) <- diag(a) + c(m /(l * 2 * esw), 0) #on n'oublie pas le 0

    # ecrire inverse de a
    inv_a <- solve(a)
  }

  # equation de krigeage
  lambda <- cbind(cmat, rep(1, nrow(distmat_xy))) %*% inv_a

  # moyenne
  mean_pred <- as.numeric(lambda %*% matrix(c(y /(l * 2 * esw), 0), ncol = 1))
  mean_pred <- ifelse(mean_pred < 0, 0, mean_pred)
  if(saturate) {
    mean_pred <- ifelse(mean_pred > quantile(mean_pred, probs = 0.999),
                        quantile(mean_pred, probs = 0.999),
                        mean_pred
    )
  }
  # variance
  var_pred <- as.numeric(sill - diag(cmat %*% t(lambda[, 1:nbpts])) + lambda[, nbpts + 1])
  var_pred <- ifelse(var_pred < 0, 0, var_pred)

  predict_krig_df <- data.frame(mean_pred = round(mean_pred, 3),
                                var_pred = round(var_pred + vario_model$nugget, 3))

  return(predict_krig_df = predict_krig_df)
}
