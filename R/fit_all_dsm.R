#' @export

get_k_best_models <- function(tab_model, k = 5, use_loo) {

  if(use_loo) {

    writeLines("\t* Using leave-one out for model selection")
    tab_model <- tab_model[order(tab_model$looic, decreasing = FALSE), ]

  } else {

    writeLines("\t* Using AIC for model selection")
    tab_model <- tab_model[order(tab_model$AIC, decreasing = FALSE), ]

  }

  tab_model <- tab_model[1:k, "index"]
  return(tab_model)
}

#' @importFrom arm rescale
#' @importFrom utils combn
#' @importFrom mgcv nb tw
#' @importFrom dsm dsm
#' @importFrom mvtnorm rmvnorm
#' @importFrom loo loo.matrix stacking_weights
#' @importFrom tweedie dtweedie
#' @export

fit_all_dsm <- function(distFit = NULL,
                        segdata_obs,
                        obsdata,
                        response = "ind",
                        predictors,
                        likelihood = "negbin",
                        esw = NULL,
                        max_cor = 0.5,
                        nb_max_pred = 3,
                        complexity = 4,
                        smooth_xy = TRUE,
                        k = 5,
                        splines_by = NULL,
                        weighted = FALSE,
                        random = NULL, # a character vector
                        soap = list(xt = NULL, knots = NULL), # a list with xt = list(bnd = ) and knots,
                        use_loo = FALSE # model selection with leave-one-out
) {
  ## likelihood must be one of "negbin", "poisson" or "tweedie"
  ## default is "negbin" for negative binomial likelihood
  ## smooth_xy controls the intercept to include or not a bivariate smooth on x and y :
  # for prediction inside the prospected polygon,
  # should be set to TRUE to obtain stable estimates
  # for prediction outside, MUST be set to FALSE to keep extrapolation under control
  ## soap is a list to pass in order to use a soap-film smooth: must be prepared outside this function
  ## k is the number of models to return for inference (based on AIC or Deviance)
  ## by default, use cubic B-splines with shrinkage
  ## arg 'splines_by' is to have an interaction with splines

  rescale <- function(x) { (x - mean(x)) / sd(x) }
  ## raw data
  X <- segdata_obs

  ## prepare smooth terms
  if(is.null(splines_by)) {
    smoothers <- paste("s(", predictors, ", k = ", complexity, ", bs = 'cs')", sep = "")
  } else {
    if(!is.character(splines_by)) {
      stop("Arg. 'splines_by' must be character")
    }
    if(!any(names(segdata_obs) == splines_by)) {
      stop(paste("No column named '", splines_by, "' in dataframe 'segdata_obs'", sep = ""))
    }
    smoothers <- paste("s(", predictors, ", k = ", complexity, ", bs = 'cs', by = ", splines_by, ")", sep = "")
  }

  ## need soap?
  if(is.null(soap$xt) && is.null(soap$knots)) {
    intercept <- ifelse(smooth_xy, "~ te(longitude, latitude, bs = 'cs')", "~ 1")
    ## standardize
    segdata_obs[, c(predictors, "longitude", "latitude", "X", "Y")] <- apply(segdata_obs[, c(predictors, "longitude", "latitude", "X", "Y")], 2, rescale)
  } else {
    intercept <- ifelse(smooth_xy, "~ s(longitude, latitude, bs = 'so', xt = list(bnd = bnd))", "~ 1")
    ## standardize: do not standardize long/lat with soap
    segdata_obs[, predictors] <- apply(segdata_obs[, predictors], 2, rescale)
  }

  ## include random effects: must be factors for mgcv
  if(!is.null(random)) {
    if(!all(random %in% names(segdata_obs))) {
      stop("Check random effect: no matching column in table 'segdata_obs'")
    } else {
      if(!is.factor(segdata_obs[, random[1]])) {
        segdata_obs[, random[1]] <- factor(as.character(segdata_obs[, random[1]]), levels = c(unique(segdata_obs[, random[1]]), "new_level"))
        X[, random[1]] <- factor(as.character(X[, random[1]]), levels = c(unique(X[, random[1]]), "new_level"))
      }
      intercept <- paste(intercept, " + s(", random[1], ", bs = 're')", sep = "")
      if(length(random) > 1) {
        for(k in 1:lenght(random)) {
          if(!is.factor(segdata_obs[, random[1]])) {
            segdata_obs[, random[k]] <- factor(as.character(segdata_obs[, random[k]]), levels = c(unique(segdata_obs[, random[k]]), "new_level"))
            X[, random[k]] <- factor(as.character(X[, random[k]]), levels = c(unique(X[, random[k]]), "new_level"))
          }
          intercept <- paste(intercept, " + s(", random[k], ", bs = 're')", sep = "")
        }
      }
    }
  }

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
  all_mods <- c(paste("count", intercept, sep = " "),
                unlist(lapply(1:nb_max_pred, mlist, y = "count", predictors = smoothers))
  )

  ## remove combinations of variables that are too correlated
  all_mods <- all_mods[which(rm_combn < max_cor)]
  ## compute weights
  if(weighted) {
    w <- lapply(all_x, function(tab) {
      sapply(1:ncol(tab), function(j) {
        make_cfact_2(calibration_data = segdata_obs,
                     test_data = segdata_obs,
                     var_name = covariable[tab[, j]],
                     percent = FALSE,
                     near_by = TRUE
        )
      })
    })
    w <- cbind(rep(1, nrow(segdata_obs)), do.call('cbind', w))
    w <- w[, which(rm_combn < max_cor)]

  } else {
    w <- matrix(1, nrow = nrow(segdata_obs), ncol = length(which(rm_combn < max_cor)))
  }

  # suppress warnings
  # options(warn = -1)

  # relabel Sample.Label to have only unique match per segments
  X$Sample.Label <- paste(X$Sample.Label, X$Seg, sep = "_")
  segdata_obs$Sample.Label <- paste(segdata_obs$Sample.Label, segdata_obs$Seg, sep = "_")
  obsdata$Sample.Label <- paste(obsdata$Sample.Label, obsdata$Seg, sep = "_")

  ## response variable is either n (nb of observations) or y (nb of individuals)
  if(response != "ind") {
    writeLines("\t* Response variable is the number of observations")
    obsdata$size <- 1
  } else {
    writeLines("\t* Response variable is the number of individuals")
  }

  ## detection
  if(is.null(distFit) && is.null(esw)) {
    stop("Must provide either a detection function as 'distFit', or esw")
  } else {
    if(!is.null(distFit)) {
      writeLines("\t* Detection function provided")
      esw <- NULL
    } else {
      if(length(esw) == nrow(segdata_obs)) {
        writeLines("\t* esw provided for each segment")
      } else {
        esw <- esw[1]
        writeLines(paste("\t* esw set to", esw, sep = " "))
      }
    }
  }
  ## fit the models
  my_dsm_fct <- function(x, tab = TRUE, segdata_obs, loo = FALSE, bnd = soap$xt, knots = soap$knots) {
    model <- dsm(as.formula(all_mods[x]),
                 ddf.obj = distFit,
                 strip.width = esw,
                 segment.data = segdata_obs,
                 observation.data = obsdata,
                 family = switch(likelihood,
                                 negbin = nb(),
                                 poisson = poisson(),
                                 tweedie = tw()
                 ),
                 method = "REML",
                 weights = w[, x],
                 knots = knots
    )
    ### leave-one-out cross-validation
    if(loo) {

      # if loo, do not return tab
      tab <- FALSE
      # approximate posterior distribution with MV normal
      beta <- mvtnorm::rmvnorm(1e3, mean = model$coefficients, sigma = model$Vp)
      Z <- predict(model,
                   newdata = model$model,
                   off.set = model$offset,
                   type = "lpmatrix"
      )
      mu <- exp(as.matrix(beta %*% t(Z)))
      # response variable from fitted dsm object
      y <- model$model$count
      # log pointwise posterior density
      lppd = switch(likelihood,
                    negbin = { apply(mu, 1, function(iter) { w[, x] * dnbinom(y, size = model$family$getTheta(trans = TRUE), mu = exp(model$offset) * iter, log = TRUE) }) },
                    poisson = { apply(mu, 1, function(iter) { w[, x] * dpois(y, lambda = exp(model$offset) * iter, log = TRUE) }) },
                    tweedie = { apply(mu, 1, function(iter) { w[, x] * log(tweedie::dtweedie(y, xi = model$family$getTheta(trans = TRUE), mu = exp(model$offset) * iter, phi = model$sig2)) }) }
      )
      out <- loo::loo.matrix(t(lppd), save_psis = TRUE)

    } else {

      out <- model

    }
    ### store some results in a data frame
    if(tab) {

      return(
        data.frame(model = all_mods[x],
                   index = x,
                   Convergence = ifelse(model$converged, 1, 0),
                   AIC = model$aic,
                   # GCV = model$gcv.ubre,
                   ResDev = model$deviance,
                   NulDev = model$null.deviance,
                   ExpDev = 100 * round(1 - model$deviance/model$null.deviance, 3)
        )
      )

    } else {

      return(out)

    }
  }

  writeLines("\t* Fitting all possible models, please wait")
  all_fits <- lapply(1:length(all_mods), my_dsm_fct, segdata_obs = segdata_obs)

  ## Collapse to a data frame
  all_fits <- do.call('rbind', all_fits)


  # 5 best looic else 5 best AIC
  if(use_loo) {

    ## leave-one-out cross-validation using Pareto Smoothing Importance Sampling
    writeLines("\t* Estimating loocv on all models: please wait")
    all_psis <- lapply(1:length(all_mods), my_dsm_fct, segdata_obs = segdata_obs, loo = TRUE)

    # get loo_ic and se_looic
    loo_ic_coefs <- do.call(rbind,lapply(all_psis, function(x){x$estimates["looic",]}))
    colnames(loo_ic_coefs) <- c("looic","se_looic")
    all_fits <- cbind(all_fits, loo_ic_coefs)

    # initiate stacking_weights column
    all_fits$stacking_weights <- NA

    # all_fits filtered 5 best
    all_fits_best <- all_fits %>%
      arrange(looic) %>%
      top_n(n=-k, wt=looic)

    # get best loo order
    index_order_best <- get_k_best_models(tab_model = all_fits_best, k = k, use_loo = T)

    # estimating stacking weights
    get_elpd_loo <- do.call('cbind', lapply(index_order_best, function(l) {all_psis[[l]]$pointwise[, "elpd_loo"]}))
    loow <- as.numeric(loo::stacking_weights(get_elpd_loo))
    all_fits$stacking_weights[index_order_best] <- loow

    all_fits_best_sw <- all_fits %>%
      arrange(desc(stacking_weights)) %>%
      top_n(n=k, wt=stacking_weights)

    index_order_sw <- all_fits_best_sw$index

    ## select the n-best models
    best <- lapply(index_order_sw, my_dsm_fct, tab = FALSE, segdata_obs = X)
    best_std <- lapply(index_order_sw, my_dsm_fct, tab = FALSE, segdata_obs = segdata_obs)



  } else {

    # all_fits filtered 5 best
    all_fits_best <- all_fits %>%
      arrange(AIC) %>%
      top_n(n=-k, wt=AIC)

    # initiate stacking_weights column
    all_fits$stacking_weights <- NA

    index_order_best <- get_k_best_models(tab_model = all_fits_best, k = k, use_loo = F)

    # estimating stacking weights
    writeLines("\t* Estimating loocv on k best models: please wait")
    all_psis <- lapply(index_order_best, my_dsm_fct, segdata_obs = segdata_obs, loo = TRUE)

    get_elpd_loo <- do.call('cbind', lapply(1:k, function(l) {all_psis[[l]]$pointwise[, "elpd_loo"]}))
    loow <- as.numeric(loo::stacking_weights(get_elpd_loo))
    all_fits$stacking_weights[index_order_best] <- loow

    all_fits_best_sw <- all_fits %>%
      arrange(desc(stacking_weights)) %>%
      top_n(n=k, wt=stacking_weights)

    index_order_sw <- all_fits_best_sw$index

    ## select the n-best models
    best <- lapply(index_order_sw, my_dsm_fct, tab = FALSE, segdata_obs = X)
    best_std <- lapply(index_order_sw, my_dsm_fct, tab = FALSE, segdata_obs = segdata_obs)

  }


  ## wrap-up with the outputs
  return(
    list(all_fits_binded = all_fits,
         best_models = best,
         best_models4plotting = best_std # pour le pred splines
    )
  )
}
