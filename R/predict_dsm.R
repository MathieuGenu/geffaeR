#' Predict density from dsm models estimated with fit_all_dsm.
#'
#' @param listdsm Whole object \code{fit_all_dsm}.
#' @param predata One of object built from \code{prep_predata} function corresponding to the data from which to predict.
#'
#' @return
#'
#' @examples
#' @import dplyr
#'
#' @export
predict_all <- function(listdsm, predata) {

  # listdsm is the output from fit_all_dsm
  # predata is the dataframe from which to predict
  ## get models
  dsmodels <- listdsm$best_models
  k <- length(dsmodels)
  ## same order as in fit_all_dsm to ensure indexes are correct
  w <- listdsm$all_fits_binded %>%
    arrange(desc(stacking_weights)) %>%
    top_n(n = k, wt = stacking_weights)
  ## predict
  writeLines(paste("\t* Predicting from ", k, " best models: please wait", sep = ""))
  out <- lapply(dsmodels, predict, newdata = predata, off.set = 1, se = TRUE)

  mean_pred <- do.call('cbind', lapply(out, function(l) { as.numeric(l$fit) }))
  se_pred <- do.call('cbind', lapply(out, function(l) { as.numeric(l$se.fit) }))

  ## predict
  writeLines("\t* Stacking predictions and done ;)")
  stack_pred <- data.frame(mean = drop(mean_pred %*% matrix(w$stacking_weights, ncol = 1)),
                           se = sqrt(drop(se_pred^2 %*% matrix(w$stacking_weights, ncol = 1)))) %>%
    mutate(cv = round(se / mean, 3)) %>%
    as.data.frame()

  # Output
  return(
    list(modelpredictions = list(mean = mean_pred, se = se_pred, cv = round(se_pred / mean_pred, 3)),
         stackedprediction = stack_pred
    )
  )
}
