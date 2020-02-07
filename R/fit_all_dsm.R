


#' @export

fit_all_dsm <- function(distFit, segdata, obsdata, outcome, predictors,
                        family = "negative binomial", esw = NULL,
                        max_cor = 0.7, nb_max_pred = 3, complexity = 3,
                        longlat = TRUE, use_gam = FALSE) {
  ## family must be one of "negative binomial", "poisson" or "tweedie"
  ## default is "negative binomial"
  ## longlat controls what sort of intercept is being fitted :
  # for prediction inside the prospected polygon,
  # should be set to TRUE to obtain stable estimates
  # for prediction outside, MUST be set to FALSE to keep extrapolation under control
  rescale <- function(x) { (x - mean(x))/sd(x) }

  ## design matrix
  X <- segdata[, c(predictors, "longitude", "latitude")]

  ## standardize
  segdata[, c(predictors, "longitude", "latitude")] <- apply(X, 2, rescale)

  ## prepare smooth terms
  smoothers <- paste("s(", predictors, ", k = ", complexity, ")", sep = "")
  intercept <- ifelse(longlat, "~ s(longitude, latitude)", "~ 1")

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

  ## Create list of models
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

  # suppress warnings
  options(warn=-1)

  ## fit the models
  if(!is.null(distFit)) {
    my_dsm_fct <- function(x) {
      if(family == "negative binomial") {
        model <- dsm(as.formula(x), distFit, segdata, obsdata, family = nb(), method = "REML")
      }
      if(family == "poisson") {
        model <- dsm(as.formula(x), distFit, segdata, obsdata, family = poisson(), method = "REML")
      }
      if(family == "tweedie") {
        model <- dsm(as.formula(x), distFit, segdata, obsdata, family = tw(), method = "REML")
      }
      ### store some results in a data frame
      data.frame(model = x,
                 Convergence = ifelse(model$converged, 1, 0),
                 AIC = model$aic,
                 # GCV = model$gcv.ubre,
                 ResDev = model$deviance,
                 NulDev = model$null.deviance,
                 ExpDev = 100 * round(1 - model$deviance/model$null.deviance, 3)
      )
    }
  } else {
    if(is.null(esw)) {
      stop("Must Provide a value for esw")
    } else {
      if(length(esw) == 1 | length(esw) == nrow(segdata)) {
        if(use_gam == TRUE) {
          if(any(names(segdata) == outcome) == FALSE) {
            countdata <- as.data.frame(obsdata[, c("Sample.Label", "size")] %>% group_by(Sample.Label) %>% summarize(n = n(), count = sum(size)))
            segdata$n <- merge(segdata, countdata, by  = "Sample.Label", all.x = TRUE)[, "n"]
            segdata$n <- ifelse(is.na(segdata$n), 0, segdata$n)
            segdata$count <- merge(segdata, countdata, by  = "Sample.Label", all.x = TRUE)[, "count"]
            segdata$count <- ifelse(is.na(segdata$count), 0, segdata$count)
          }
          my_dsm_fct <- function(x) {
            if(family == "negative binomial") {
              model <- gam(as.formula(x), data = segdata, offset = 2*esw*segdata$Effort, family = nb(), method = "REML")
            }
            if(family == "poisson") {
              model <- gam(as.formula(x), data = segdata, offset = 2*esw*segdata$Effort, family = poisson(), method = "REML")
            }
            if(family == "tweedie") {
              model <- gam(as.formula(x), data = segdata, offset = 2*esw*segdata$Effort, family = tw(), method = "REML")
            }
            ### store some results in a data frame
            data.frame(model = x,
                       Convergence = ifelse(model$converged, 1, 0),
                       AIC = model$aic,
                       # GCV = model$gcv.ubre,
                       ResDev = model$deviance,
                       NulDev = model$null.deviance,
                       ExpDev = 100*round(1 - model$deviance/model$null.deviance, 3)
            )
          }
        } else {
          my_dsm_fct <- function(x) {
            if(family == "negative binomial") {
              model <- dsm(as.formula(x), ddf.obj = distFit, strip.width = esw, segment.data = segdata, observation.data = obsdata, family = nb(), method = "REML")
            }
            if(family == "poisson") {
              model <- dsm(as.formula(x), ddf.obj = distFit, strip.width = esw, segment.data = segdata, observation.data = obsdata, family = poisson(), method = "REML")
            }
            if(family == "tweedie") {
              model <- dsm(as.formula(x), ddf.obj = distFit, strip.width = esw, segment.data = segdata, observation.data = obsdata, family = tw(), method = "REML")
            }
            ### store some results in a data frame
            data.frame(model = x,
                       Convergence = ifelse(model$converged, 1, 0),
                       AIC = model$aic,
                       # GCV = model$gcv.ubre,
                       ResDev = model$deviance,
                       NulDev = model$null.deviance,
                       ExpDev = 100*round(1 - model$deviance/model$null.deviance, 3)
            )
          }
        }
      } else {
        stop("Please check esw: provide either a single value or a vector with same length as segdata")
      }
    }
  }
  all_fits <- lapply(all_mods, my_dsm_fct)
  ## Collapse to a data frame
  all_fits_binded <- do.call(rbind, all_fits)
  return(all_fits_binded = all_fits_binded)
}
