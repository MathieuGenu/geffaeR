#' @export

get_k_best_models <- function(tab_model, k = 5, use_AIC = TRUE) {
  if(use_AIC) {
    tab_model <- tab_model[order(tab_model$AIC, decreasing = FALSE), ]
  } else {
    writeLines("\tUsing Explained Deviance for model selection")
    tab_model <- tab_model[order(tab_model$ExpDev, decreasing = TRUE), ]
  }
  return(tab_model[1:k, "index"])
}


#' @export

fit_all_dsm <- function(distFit, segdata, obsdata, outcome, predictors,
                        likelihood = "negbin", esw = NULL,
                        max_cor = 0.5, nb_max_pred = 3, complexity = 4,
                        smooth_xy = TRUE, use_gam = FALSE,
                        k = 5, weighted = FALSE
                        ) {
  ## likelihood must be one of "negbin", "poisson" or "tweedie"
  ## default is "negbin" for negative binomial likelihood
  ## smooth_xy controls the intercept to include or not a bivariate smooth on x and y :
  # for prediction inside the prospected polygon,
  # should be set to TRUE to obtain stable estimates
  # for prediction outside, MUST be set to FALSE to keep extrapolation under control
  ## k is the number of model to return (based on AIC or Deviance)
  ## by default, use cubic B-splines with shrinkage
  rescale <- function(x) { (x - mean(x)) / sd(x) }

  ## raw data
  X <- segdata

  ## standardize
  segdata[, c(predictors, "longitude", "latitude", "X", "Y")] <- apply(segdata[, c(predictors, "longitude", "latitude", "X", "Y")], 2, rescale)

  ## prepare smooth terms
  smoothers <- paste("s(", predictors, ", k = ", complexity, ", bs = 'cs')", sep = "")
  intercept <- ifelse(smooth_xy, "~ s(X, Y)", "~ 1")

  ## all combinations among nb_max_pred
  all_x <- lapply(1:nb_max_pred, combn, x = length(predictors))

  ## check whether cross-correlation needs to be evaluated
  if(nb_max_pred == 1) {
    rm_combn <- c(rep(0, length(predictors) + 1))
  } else {
    ## identify which combination is to be removed
    rm_combn <- lapply(all_x[-1], function(mat) {
      sapply(1:ncol(mat), function(i) {
        rho <- cor(X[, predictors[mat[, i]]]) ; diag(rho) <- 0
        return(max(abs(as.numeric(rho))))
      })
    })
    rm_combn <- c(c(rep(0, length(predictors) + 1)), unlist(rm_combn))
  }

  ## list of models
  mlist <- function(n, y, predictors) {
    paste(y,
          apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + "),
          sep = paste(intercept, "+", sep = " ")
          )
  }
  all_mods <- c(paste(outcome, intercept, sep = " "),
                unlist(lapply(1:nb_max_pred, mlist, y = outcome, predictors = smoothers))
                )

  ## remove combinations of variables that are too correlated
  all_mods <- all_mods[which(rm_combn < max_cor)]
  ## compute weights
  if(weighted) {
    w <- lapply(all_x, function(tab) {
      sapply(1:ncol(tab), function(j) {
        make_cfact_2(calibration_data = segdata,
                     test_data = segdata,
                     var_name = covariable[tab[, j]],
                     percent = F,
                     near_by = T
                     )
      })
    })
    w <- cbind(rep(1, nrow(segdata)), do.call('cbind', w))
    w <- w[, which(rm_combn < max_cor)]

  } else {
    w <- matrix(1, nrow = nrow(segdata), ncol = length(which(rm_combn < max_cor)))
  }

  # suppress warnings
  # options(warn = -1)

  # relabel Sample.Label to have only unique match per segments
  # X$Sample.Label <- paste(X$Sample.Label, X$Seg, sep = "_")
  # segdata$Sample.Label <- paste(segdata$Sample.Label, segdata$Seg, sep = "_")
  # obsdata$Sample.Label <- paste(obsdata$Sample.Label, obsdata$Seg, sep = "_")

  ## fit the models
  if(!is.null(distFit)) {
    my_dsm_fct <- function(x, tab = TRUE, segdata) {
      model <- dsm(as.formula(all_mods[x]), 
                   ddf.obj = distFit, 
                   segment.data = segdata, 
                   observation.data = obsdata, 
                   strip.width = NULL, 
                   family = switch(likelihood, 
                                   negbin = nb(),
                                   poisson = poisson(),
                                   tweedie = tw()
                                   ), 
                   method = "REML", 
                   weights = w[, x]
                   )
      ### store some results in a data frame
      if(tab) {
        return(data.frame(model = all_mods[x],
                          index = x,
                          Convergence = ifelse(model$converged, 1, 0),
                          AIC = model$aic,
                          # GCV = model$gcv.ubre,
                          ResDev = model$deviance,
                          NulDev = model$null.deviance,
                          ExpDev = 100 * round(1 - model$deviance/model$null.deviance, 3)
                          )
        )
      } else { return(model) }
    }
  } else {
    if(is.null(esw)) {
      stop("Must Provide a value for esw")
    } else {
      if(length(esw) == 1 | length(esw) == nrow(segdata)) {
        if(use_gam == TRUE) {
          my_dsm_fct <- function(x, tab = TRUE, segdata) {
            model <- gam(as.formula(all_mods[x]), 
                         data = segdata, 
                         offset = 2 * esw * segdata$Effort, 
                         family = switch(likelihood, 
                                         negbin = nb(),
                                         poisson = poisson(),
                                         tweedie = tw()
                                         ), 
                         method = "REML", 
                         weights = w[, x]
                         )
            ### store some results in a data frame
            if(tab) {
              return(data.frame(model = all_mods[x],
                                index = x,
                                Convergence = ifelse(model$converged, 1, 0),
                                AIC = model$aic,
                                # GCV = model$gcv.ubre,
                                ResDev = model$deviance,
                                NulDev = model$null.deviance,
                                ExpDev = 100 * round(1 - model$deviance/model$null.deviance, 3)
                                )
              )
            } else { return(model) }
          }
        } else {
          my_dsm_fct <- function(x, tab = TRUE, segdata) {
            model <- dsm(as.formula(all_mods[x]), 
                         ddf.obj = NULL, 
                         strip.width = esw, 
                         segment.data = segdata, 
                         observation.data = obsdata, 
                         family = switch(likelihood, 
                                         negbin = nb(),
                                         poisson = poisson(),
                                         tweedie = tw()
                                         ), 
                         method = "REML", 
                         weights = w[, x]
                         )
            ### store some results in a data frame
            if(tab) {
              return(data.frame(model = all_mods[x],
                                index = x,
                                Convergence = ifelse(model$converged, 1, 0),
                                AIC = model$aic,
                                # GCV = model$gcv.ubre,
                                ResDev = model$deviance,
                                NulDev = model$null.deviance,
                                ExpDev = 100 * round(1 - model$deviance/model$null.deviance, 3)
                                )
              )
            } else { return(model) }
          }
        }
      } else {
        stop("Please check esw: provide either a single value or a vector with same length as segdata")
      }
    }
  }
  all_fits <- lapply(1:length(all_mods), my_dsm_fct, segdata = segdata)
  ## Collapse to a data frame
  all_fits_binded <- do.call('rbind', all_fits)

  ## select the n-best models
  best <- lapply(get_k_best_models(tab_model = all_fits_binded, k = k), my_dsm_fct, tab = FALSE, segdata = X)
  best_std <- lapply(get_k_best_models(tab_model = all_fits_binded, k = k), my_dsm_fct, tab = FALSE, segdata = segdata)

  ## wrap-up with the outputs
  return(list(all_fits_binded = all_fits_binded,
              best_models = best,
              best_models_std = best_std # pour le pred splines
              )
         )
}
