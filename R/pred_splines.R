#' @export

### spline curves
pred_splines <- function(segdata, dsm_model, remove_intercept = FALSE) {
  # segdata is the original dataset used to calibrate the model
  # dsm_model is the dsm model you want to use
  lower <- function(x) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = 0.95)[1]) }
  upper <- function(x) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = 0.95)[2]) }

  var_name <- as.character(dsm_model$pred.formula)[grep(" + ", as.character(dsm_model$pred.formula), fixed = "TRUE")]
  var_name <- strsplit(var_name, split = " + ", fixed = TRUE)[[1]]
  if(any(var_name %in% c("off.set", "offset"))) { var_name <- var_name[-which(var_name %in% c("off.set", "offset"))] }

  nx <- 1e3
  n_sim <- 1e4

  Xnew <- as.data.frame(sapply(var_name, function(id) { rep(0, nx) }))
  df_splines <- NULL

  ### approximate posterior distribution with MV normal
  beta <- mvtnorm::rmvnorm(n_sim, mean = dsm_model$coefficients, sigma = dsm_model$Vp)

  if(remove_intercept) { beta[, 1] <- 0; writeLines("\tRemoving intercept") }

  rm_spatial <- grep("X,Y", names(dsm_model$coefficients))
  if(length(rm_spatial) != 0) { beta[, rm_spatial] <- 0.0 }
  if(any(var_name %in% c("X", "Y"))) { var_name <- var_name[-which(var_name %in% c("X", "Y"))] }

  for(j in var_name) {
    Z <- Xnew
    Z[, j] <- (seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx) - mean(segdata[, j], na.rm = TRUE)) / sd(segdata[, j], na.rm = TRUE)
    Z <- predict(dsm_model, newdata = Z, off.set = 1, type = "lpmatrix")
    linpred <- beta %*% t(Z); rm(Z)
    df_splines <- rbind(df_splines,
                        data.frame(x = seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx),
                                   y = apply(linpred, 2, mean),
                                   lower = apply(linpred, 2, lower),
                                   upper = apply(linpred, 2, upper),
                                   param = rep(j, nx),
                                   scale = rep("log", nx)
                        ),
                        data.frame(x = seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx),
                                   y = apply(exp(linpred), 2, mean),
                                   lower = apply(exp(linpred), 2, lower),
                                   upper = apply(exp(linpred), 2, upper),
                                   param = rep(j, nx),
                                   scale = rep("natural", nx)
                        )
    )
  }
  g_splines <- ggplot(data = df_splines,
                      aes(x = x, y = y, ymin = lower, ymax = upper)
                      ) +
    geom_ribbon(alpha = 0.3, fill = "midnightblue") +
    geom_line(color = "midnightblue") +
    facet_grid(scale ~ param, scales = "free") +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = "Covariate") +
    theme_bw()

  return(list(df_splines = df_splines,
              g_splines = g_splines
  )
  )
}
