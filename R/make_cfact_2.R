#' @export

make_cfact_2 <- function(calibration_data,
                         test_data,
                         var_name = NULL,
                         howmany = 1,
                         eps = 6,
                         rm.dup.test = FALSE,
                         percent = TRUE,
                         near_by = FALSE
                         ) {

  ### custom code
  if(is.null(var_name)) { var_name = names(calibration_data) }

  ## standardize new data to predict from
  ### useful functions
  rescale <- function(x) { return((x - mean(x)) / sd(x)) }
  rescale2 <- function(ynew, y) { return((ynew - mean(y, na.rm = TRUE)) /(sd(y, na.rm = TRUE))) }
  # this simplifies computation A LOT!
  make_X <- function(calibration_data, test_data, var_name){
    X <- sapply(var_name,
                function(k) { rescale2(ynew = test_data[, k],
                                       y = calibration_data[, k]
                )}
    )
    X <- as.data.frame(X)
    names(X) <- var_name
    return(X)
  }
  ### standardize
  Xcal = make_X(calibration_data = calibration_data, test_data = calibration_data, var_name)
  Xtest = make_X(calibration_data = calibration_data, test_data = test_data, var_name)

  # Round the standardized values
  Xcal <- round(Xcal, eps) ; Xtest <- round(Xtest, eps)

  # Remove duplicates
  dup <- duplicated(Xcal[, var_name])
  Xcal <- Xcal[dup == FALSE, ] ; rm(dup)
  if(rm.dup.test) {
    dup <- duplicated(Xtest[, var_name])
    Xtest <- Xtest[dup == FALSE, ] ; rm(dup)
  }

  # rename rows
  if(is.null(dim(Xcal))) {
    Xcal <- as.data.frame(Xcal)
    names(Xcal) <- var_name
  }
  if(is.null(dim(Xtest))) {
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- var_name
  }
  row.names(Xcal) <- 1:nrow(Xcal)
  row.names(Xtest) <- 1:nrow(Xtest)

  # compute counterfactuals
  if(near_by) {
    cf <- whatif(formula = NULL, data = Xcal, cfact = Xtest, choice = "distance", nearby = howmany, mc.cores = 1)
    if(percent) {
      out <- round(mean(cf$geom.var), 3)
    } else {
      out <- as.numeric(cf$sum.stat)
      out <- out / mean(out)
    }
  } else {
    cf <- whatif(formula = NULL, data = Xcal, cfact = Xtest, choice = "hull", mc.cores = 1)
    cf <- ifelse(cf$in.hull, 0, 1)
    if(percent) {
      out <- round(mean(cf), 3)
    } else {
      out <- cf
    }
  }
  return(out)
}
