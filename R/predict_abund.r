#' Predict density from dsm models estimated with fit_all_dsm.
#'
#' @param listdsm Whole object \code{fit_all_dsm}.
#' @param predata One of object built from \code{prep_predata} function corresponding to the data from which to predict.
#' @param n_sim number of draws to approximate posterior distribution with MV normal distribution.
#' @param abund_by grouping factor for abundance estimation.
#' @param alpha coverage level for confidence interval.
#'
#' @return
#'
#' @examples
#' @import dplyr tidyr
#' @importFrom coda HPDinterval as.mcmc
#' @importFrom mvtnorm rmvnorm
#'
#' @export
predict_abund <- function(listdsm, predata, n_sim = 5e2, abund_by = NULL, alpha = 0.8) {

  # listdsm is the output from fit_all_dsm
  # predata is the dataframe from which to predict

  ## check
  if(!any(names(predata) == "Area")) {
  	stop("\t 'predata' must include a column 'Area' with area information")
  }
  if(!is.null(abund_by)) {
  	if(!any(names(predata) == abund_by)) {
  	  stop("\t 'predata' must include a column with same name as 'abund_by'")
    } else {
      if(!is.factor(predata[, abund_by])) {
      	writeLines("\t\t forcing conversion to factor: please check levels")
      	predata[, abund_by] <- factor(predata[, abund_by], levels = unique(predata[, abund_by]))
      }
      nom <- levels(predata[, abund_by])
      X <- sapply(levels(predata[, abund_by]), function(x) { ifelse(predata[, abund_by] == x, 1, 0) })
    }
  } else {
  	nom <- "abundance"
  	X <- rep(1, nrow(predata))
  }
  
  ## useful functions
  lower <- function(x) {
    if(all(is.na(x))) {
      return(NA)
    } else {
      x <- x[!is.na(x)]
      return(as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[1]) )
    }
  }
  upper <- function(x) {
    if(all(is.na(x))) {
      return(NA)
    } else {
      x <- x[!is.na(x)]
      return(as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[2]) )
    }
  }

  ## get models
  dsmodels <- listdsm$best_models
  k <- length(dsmodels)
  
  ## same order as in fit_all_dsm to ensure indexes are correct
  w <- listdsm$all_fits_binded %>%
    arrange(desc(stacking_weights)) %>%
    top_n(n = k, wt = stacking_weights)

  ## predict
  writeLines(paste("\t* Predicting from ", k, " best models: please wait", sep = ""))
  ### approximate posterior distribution with MV normal
  # output is a dataframe with iterations as rows and three columns 'grouping', 'estimate' and 'iter'
  out <- lapply(dsmodels, function(dsm_model) {
    dd <- exp(mvtnorm::rmvnorm(n_sim, mean = dsm_model$coefficients, sigma = dsm_model$Vp) %*% 
                t(predict(dsm_model, newdata = predata, off.set = 1, type = "lpmatrix")))
    # handles NA
    rm_na <- sapply(1:ncol(dd), function(j) { !all(is.na(dd[, j])) })
    dd <- dd[, rm_na] %*% 
      diag(predata$Area[rm_na]) %*% # scale to area
      X[rm_na, ] %>% # grouping
      as.matrix() %>%
      `colnames<-`(nom) %>%
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "grouping", values_to = "estimate") %>%
      arrange(grouping) %>%
      mutate(iter = rep(1:n_sim, times = length(nom))) %>% # check this is correct
      as.data.frame()
    return(dd)
  } 
  )
  # reorganize to a long format data.frame
  out <- do.call('rbind', out) %>%
    mutate(model = rep(1:k, each = n_sim * length(nom))) %>%
    left_join(w %>% mutate(model = 1:k) %>% select(model, stacking_weights), # need to specify correct column in w here, must be ordered correctly
              by = "model"
              ) %>% 
    mutate(stacked_estimate = stacking_weights * estimate) %>%
    as.data.frame()

  ## predict
  writeLines("\t* Stacking predictions and done ;)")
  best_of_k_models_df <- out %>%
    group_by(grouping, model, iter) %>%
    summarize(estimate = sum(estimate)) %>%
    group_by(grouping, model) %>%
    summarize(mean = round(mean(estimate), 0),
              se = round(sd(estimate), 0),
              lower = round(lower(estimate), 0),
              upper = round(upper(estimate), 0)
              ) %>%
    as.data.frame()

  stack_pred <- out %>%
    group_by(grouping, iter) %>%
    summarize(estimate = sum(stacked_estimate)) %>%
    mutate(model = "stacked") %>%
    group_by(grouping, model) %>%
    summarize(mean = round(mean(estimate), 0),
              se = round(sd(estimate), 0),
              lower = round(lower(estimate), 0),
              upper = round(upper(estimate), 0)
              ) %>%
    as.data.frame()

  # merge all models
  df_predict_all <- rbind(stack_pred, best_of_k_models_df)

  # Output
  return(df_predict_all)
}
