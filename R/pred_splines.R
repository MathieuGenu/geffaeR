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

  # rm_spatial <- grep("X,Y", names(dsm_model$coefficients))
  rm_spatial <- grep("longitude,latitude", names(dsm_model$coefficients))
  if(length(rm_spatial) != 0) {
    # gridxy <- expand.grid(X = (seq(min(segdata[, "X"], na.rm = TRUE), max(segdata[, "X"], na.rm = TRUE), length.out = 1e2) - mean(segdata[, "X"], na.rm = TRUE)) / sd(segdata[, "X"], na.rm = TRUE),
    #                       Y = (seq(min(segdata[, "Y"], na.rm = TRUE), max(segdata[, "Y"], na.rm = TRUE), length.out = 1e2) - mean(segdata[, "Y"], na.rm = TRUE)) / sd(segdata[, "Y"], na.rm = TRUE)
    #                       )
    gridxy <- expand.grid(longitude = (seq(min(segdata[, "longitude"], na.rm = TRUE), max(segdata[, "longitude"], na.rm = TRUE), length.out = 1e2) - mean(segdata[, "longitude"], na.rm = TRUE)) / sd(segdata[, "longitude"], na.rm = TRUE),
                          latitude = (seq(min(segdata[, "latitude"], na.rm = TRUE), max(segdata[, "latitude"], na.rm = TRUE), length.out = 1e2) - mean(segdata[, "latitude"], na.rm = TRUE)) / sd(segdata[, "latitude"], na.rm = TRUE)
                          )
    Z <- predict(dsm_model,
                 newdata = cbind(gridxy,
                                 # as.data.frame(sapply(var_name[-which(var_name %in% c("X", "Y"))], function(id) { rep(0, 1e4) }))
                                 as.data.frame(sapply(var_name[-which(var_name %in% c("longitude", "latitude"))], function(id) { rep(0, 1e4) }))
                                 ),
                 off.set = 1,
                 type = "lpmatrix"
                 )
    slinpred <- beta %*% t(Z); rm(Z)
    # gridxy <- expand.grid(X = seq(min(segdata[, "X"], na.rm = TRUE), max(segdata[, "X"], na.rm = TRUE), length.out = 1e2),
    #                       Y = seq(min(segdata[, "Y"], na.rm = TRUE), max(segdata[, "Y"], na.rm = TRUE), length.out = 1e2)
    #                       )
    gridxy <- expand.grid(longitude = seq(min(segdata[, "longitude"], na.rm = TRUE), max(segdata[, "longitude"], na.rm = TRUE), length.out = 1e2),
                          latitude = seq(min(segdata[, "latitude"], na.rm = TRUE), max(segdata[, "latitude"], na.rm = TRUE), length.out = 1e2)
                          )
    gridxy <- rbind(gridxy, gridxy)
    gridxy$spatial <- c(apply(slinpred, 2, median), apply(exp(slinpred), 2, mean))
    gridxy$lower <- c(apply(slinpred, 2, lower), apply(exp(slinpred), 2, lower))
    gridxy$upper <- c(apply(slinpred, 2, upper), apply(exp(slinpred), 2, upper))
    gridxy$scale <- rep(c("log", "natural"), each = 1e4)
    ### set to 0 for remaining effects
    beta[, rm_spatial] <- 0.0
  } else {
    gridxy <- NULL
  }
  if(any(var_name %in% c("X", "Y", "longitude", "latitude"))) { var_name <- var_name[-which(var_name %in% c("X", "Y", "longitude", "latitude"))] }

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
              g_splines = g_splines,
              spatial = gridxy
              )
         )
}


#' @export

rootogram_nb <- function(model_fit, n_rep = 1e3, min_obs = 0) {
  ### posterior predictive checks
  beta <- mvtnorm::rmvnorm(n_rep, mean = model_fit$coefficients, sigma = model_fit$Vp)
  Z <- predict(model_fit, type = "lpmatrix")
  linpred <- beta %*% t(Z)
  transfo_overdispersion <- exp(model_fit$family$getTheta())
  y_rep <- t(apply(exp(linpred), 1,
                   function(x) { MASS::rnegbin(n = length(x), mu = x, theta = transfo_overdispersion) }
                   )
             )

  ### rootogram
  countdata <- model_fit$data$count
  f_histogram <- function(x, max_obs) { table(factor(x, levels = min_obs:max_obs)) }
  max_obs <- max(countdata, as.numeric(y_rep))

  ### earth-mover distance
  cantonnier <- function(x, y, breaks, remove_zeroes = FALSE) {
    if(remove_zeroes) {
      x <- x[which(x != 0)]
      y <- y[which(x != 0)]
    }
    x <- hist(x, breaks, plot = FALSE)$density
    y <- hist(y, breaks, plot = FALSE)$density
    return(emd = sum(abs(cumsum(x) - cumsum(y))))
  }
  gof <- apply(y_rep, 1, cantonnier, y = countdata, breaks = min_obs:max_obs, remove_zeroes = TRUE)

  ### plot
  theme_set(theme_bw(base_size = 16))
  g <- ggplot(data.frame(mids = (min_obs:max_obs + (min_obs + 1):(max_obs + 1))/2,
                         y_obs = as.numeric(f_histogram(x = countdata, max_obs = max_obs)),
                         y_rep = apply(apply(y_rep, 1, f_histogram, max_obs = max_obs), 1, mean)
                         ),
              aes(x = mids, y = y_rep)
              ) +
    geom_rect(aes(xmin = mids - 0.5, xmax = mids + 0.5,
                  ymax = y_obs, ymin = 0),
              fill = "lightgrey", color = grey(0.6), alpha = 0.5
              ) +
    geom_line(aes(x = mids, y = y_rep), color = "black", size = 1.0) +
    scale_y_sqrt(name = "Count") +
    scale_x_sqrt(name = quote(y[obs]), breaks = c(c(1, 5, 10, 50), seq(100, 100 * ceiling(max_obs / 100), 100))) +
    guides(size = "none") +
    theme(plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          strip.background = element_rect(fill = grey(0.95))
          )

  ### wrap-up
  return(list(rootogram = g,
              earthMover = gof
              )
         )
}
