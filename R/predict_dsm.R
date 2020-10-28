#' Predict density from dsm models estimated with fit_all_dsm.
#'
#' @param listdsm Whole object \code{fit_all_dsm}.
#' @param predata One of object built from \code{prep_predata} function corresponding to the data from which to predict.
#'
#' @return
#'
#' @examples
#' @import dplyr tidyr
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


  # reorganize best model values to a data.frame in stack_pred format to merge to it
  matrix_pred <- function(list_out, var) {

    temp_pred <- do.call('cbind', lapply(list_out, function(l) { as.numeric(l[[var]]) }))
    return(temp_pred)

  }

  mean_pred <- matrix_pred(out, var = "fit")
  se_pred <- matrix_pred(out, var = "se.fit")

  list_to_dataFrame <- function(temp_pred, var, k, new_var_name) {

    temp_df <- temp_pred %>%
      `colnames<-`(paste("best", seq(1,k,1), sep = "_")) %>%
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "model", values_to = new_var_name) %>%
      arrange(model)

    return(temp_df)

  }

  mean_df <- list_to_dataFrame(temp_pred = mean_pred,
                               var = "fit",
                               k = length(dsmodels),
                               new_var_name = "mean")

  se_df <- list_to_dataFrame(temp_pred = se_pred,
                             var = "se.fit",
                             k = length(dsmodels),
                             new_var_name = "se")

  best_of_k_models_df <- cbind(mean_df, se_df[,"se"]) %>%
    mutate(cv = round(se / mean, 3)) %>%
    select(mean,se,cv,model)


  ## predict
  writeLines("\t* Stacking predictions and done ;)")
  stack_pred <- data.frame(mean = drop(mean_pred %*% matrix(w$stacking_weights, ncol = 1)),
                           se = sqrt(drop(se_pred^2 %*% matrix(w$stacking_weights, ncol = 1)))) %>%
    mutate(cv = round(se / mean, 3), model = "stacking") %>%
    as.data.frame()

  # merge all models
  df_predict_all <- rbind(stack_pred, best_of_k_models_df)


  # Output
  return(df_predict_all)
}
