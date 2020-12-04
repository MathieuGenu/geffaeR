#' Get spatial density / Abundance prediction.
#'
#' Give spatial density from adjusted covariance model.
#'
#'
#' @param segdata_obs data.frame. Output of \code{\link[geffaeR]{ajout_obs}}.
#' @param esw Effective (half) Strip-Width determined by \code{\link[geffaeR]{plot_detection}}.
#' @param vario_model output of \code{\link[geffaeR]{fit_variomodel}}.
#' @param predict_coord_xy data.frame with X and Y columns for longitude and latitude,
#' giving area where the kriging prediction is made.
#' @param fast_inversion
#' @param saturate If TRUE, allows to exclude extreme prediction by keeping values under 0.999 quantile
#' @param intraTransect Boolean. Consider only distance between point on the same segment for the kriging distance ?
#'
#' @return data.frame of mean prediction and its se (standard error) with coordinates (X and Y).
#'
#' @examples
#' @importFrom arm bayesglm
#' @importFrom fields rdist
#' @export
predict_kriging <- function(segdata_obs,
                            esw,
                            vario_model,
                            predict_coord_xy = NULL, # dataframe with columns X and Y (correctly projected)
                            fast_inversion = TRUE,
                            saturate = TRUE,
                            intraTransect = TRUE
                            ) {


  col_needed <- c("y", "n", "Effort")
  if(!all(col_needed %in% colnames(segdata_obs))) {
    missing_col <- col_needed[!(col_needed %in% colnames(segdata_obs))]
    stop(paste("Les colonnes suivantes ne sont pas dans les colonnes de segdata_obs :",
               paste(missing_col, collapse = "\n"), sep="\n"))
  }
  y <- segdata_obs$y
  l <- segdata_obs$Effort

  # distmat_x = matrice des distances entre points du jeu de données
  if(all(c("X", "Y") %in% names(segdata_obs))) {
    center <- segdata_obs[, c("X", "Y")]
  } else {
    stop("Must have X and Y coordinates (projected for accurate distance computation)")
  }
  D <- fields::rdist(center) / 1e3 # in kilometers
  # tenir compte du design de l'Ã©chantillonnage

  if(!intraTransect) {
    id_transect <- rep(1, length(y))
  } else {
    id_transect <- segdata_obs[, "Transect.Label"]
  }
  D_design <- fields::rdist(as.numeric(id_transect)) * 1000
  distmat_x <- D + D_design

  # distmat_xy = matrice des distances entre points du jeu de données et points à prédire
  if(is.null((predict_coord_xy))) {
    xlim <- 1e3 * c(floor(min(center[, 1]) / 1e3) - 1, ceiling(max(center[, 1]) / 1e3) + 1)
    ylim <- 1e3 * c(floor(min(center[, 2]) / 1e3) - 1, ceiling(max(center[, 2]) / 1e3) + 1)
    predict_coord_xy <- expand.grid(X = seq(xlim[1], xlim[2], length.out = 1e2),
                                    Y = seq(ylim[1], ylim[2], length.out = 1e2)
    )
  }

  cov_matern <- function(d, sill, range) {
    return(sill *(1 + d * sqrt(3) / range) * exp(- d * sqrt(3) / range))
  }

  cov_exp <- function(d, sill, range) {
    return(sill * exp(- d / range))
  }

  TauxObs <- function(y, esw, l){
    my_control <- list(epsilon = 1e-6, maxit = 100000, trace = FALSE)
    # utiliser quasi poisson pour eventuelle surdispersion
    model_m <- arm::bayesglm(y ~ 1 + offset(2 * esw * l),
                             prior.mean.for.intercept = -1.0, # prior is less than 4 obs/ind per km
                             prior.scale.for.intercept = 1.0,
                             prior.df.for.intercept = 7,
                             family = "quasipoisson", control = my_control
    )
    m <- as.numeric(exp(model_m$coefficients))
    return(m)
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
  if(nrow(predict_coord_xy) > 1e3) {
    n_break <- floor(nrow(predict_coord_xy) / 1e3)
    if(nrow(predict_coord_xy) %% 1e3 == 0) {
      distmat_xy <- vector(mode = 'list', length = n_break)
    } else {
      distmat_xy <- vector(mode = 'list', length = n_break + 1)
      distmat_xy[[n_break + 1]] <- fields::rdist(predict_coord_xy[c((1e3 * n_break + 1):nrow(predict_coord_xy)), ], center) / 1e3 # in kilometers
    }
    for(i in 1:n_break) {
      distmat_xy[[i]] <- fields::rdist(predict_coord_xy[1e3 * (i - 1) + c(1:1e3), ], center) / 1e3 # in kilometers
    }
  } else {
    distmat_xy <- fields::rdist(predict_coord_xy, center) / 1e3 # in kilometers
  }

  if(form == "matern") {
    # matrice des covariances entre toutes les données
    a <- cov_matern(distmat_x, sill, range)

    # matrice des covariances entre données et prédictions
    cmat <- lapply(distmat_xy, cov_matern, sill = sill, range = range)
    # cmat <- cov_matern(distmat_xy, sill = sill, range = range)

  }

  if(form == "exponential") {
    # matrice des covariances entre toutes les données
    a <- cov_exp(distmat_x, sill, range)

    # matrice des covariances entre données et prédictions
    cmat <- lapply(distmat_xy, cov_exp, sill = sill, range = range)
    # cmat <- cov_exp(distmat_xy, sill = sill, range = range)

  }

  if(fast_inversion) {
    diag(a) <- diag(a) + m /(l * 2 * esw)

    # tester si toutes les valeurs propres de a sont > 0
    if(any(eigen(a)$values < 0)) {
      stop("Non-invertible matrix, eigen values negative")
    }
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
  laplacien <- function(input_mat) {
    cbind(input_mat, rep(1, nrow(input_mat))) %*% inv_a
  }
  krig_pred <- function(input_mat) {
    as.numeric(input_mat %*% matrix(c(y /(l * 2 * esw), 0), ncol = 1))
  }
  krig_var <- function(covariance_mat, laplacian_mat) {
    if(nrow(covariance_mat) != nrow(laplacian_mat)) { stop("Dimension mismatch") }
    return(as.numeric(sill - diag(covariance_mat %*% t(laplacian_mat[, -ncol(laplacian_mat)])) + laplacian_mat[, ncol(laplacian_mat)]))
  }
  lambda <- lapply(cmat, laplacien) # matrix list

  mean_pred <- drop(do.call('c', lapply(lambda, krig_pred)))
  var_pred <- drop(do.call('c', lapply(1:length(lambda), function(id) {krig_var(covariance_mat = cmat[[id]], laplacian_mat = lambda[[id]])})))
  # lambda <- cbind(cmat, rep(1, nrow(distmat_xy))) %*% inv_a

  # moyenne
  # mean_pred <- as.numeric(lambda %*% matrix(c(y /(l * 2 * esw), 0), ncol = 1))
  mean_pred <- ifelse(mean_pred < 0, 0, mean_pred)
  if(saturate) {
    mean_pred <- ifelse(mean_pred > quantile(mean_pred, probs = 0.999, na.rm = TRUE),
                        quantile(mean_pred, probs = 0.999, na.rm = TRUE),
                        mean_pred
    )
  }
  # variance
  # var_pred <- as.numeric(sill - diag(cmat %*% t(lambda[, 1:nbpts])) + lambda[, nbpts + 1])
  var_pred <- ifelse(var_pred < 0, 0, var_pred)

  predict_krig_df <- data.frame(mean_pred = round(mean_pred, 3),
                                se_pred = round(sqrt(var_pred + vario_model$nugget), 3)
  )

  predict_krig_df <- cbind(predict_coord_xy, predict_krig_df)
  return(predict_krig_df = predict_krig_df)
}
