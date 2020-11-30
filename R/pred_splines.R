#' Prediction diagnostic.
#'
#' Visualize prediction fiiting and get
#'
#' @param segdata segdata data.frame built with \code{prepare_data_effort}
#' @param dsm_model The dsm model you want to use, built with \code{fit_all_dsm}.
#' @param remove_intercept Remove intercept of model (\code{beta[, 1] <- 0}).
#' @param random Variable to include random effect.
#' @param splines_by Interaction with splines given by one variable of segdata.
#' @param alpha
#'
#' @return \enumerate{
#'           \item df_splines : \code{data.frame}.
#'           \item g_splines : \code{ggplot}.
#'           \item spatial : \code{data.frame}.
#'         }
#'
#' @import ggplot2
#' @import glue
#' @importFrom coda HPDinterval as.mcmc
#' @importFrom mvtnorm rmvnorm
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid grobTree
#' @importFrom cowplot draw_grob ggdraw
#' @export
#' @examples
pred_splines <- function(segdata, dsm_model, remove_intercept = FALSE, random = NULL, splines_by = NULL, alpha = 0.8) {
  # segdata is the original dataset used to calibrate the model
  # dsm_model is the dsm model you want to use
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

  var_name <- as.character(dsm_model$pred.formula)[grep(" + ", as.character(dsm_model$pred.formula), fixed = "TRUE")]
  var_name <- strsplit(var_name, split = " + ", fixed = TRUE)[[1]]
  if(any(var_name %in% c("off.set", "offset", random, splines_by))) { var_name <- var_name[-which(var_name %in% c("off.set", "offset", random, splines_by))] }

  nx <- 1e3
  n_sim <- 1e4

  Xnew <- as.data.frame(sapply(var_name, function(id) { rep(0, nx) }))
  df_splines <- NULL

  ### handles splines_by arg.
  if(!is.null(splines_by)) {
    if(!any(names(dsm_model$data) == splines_by)) {
      stop(paste("No column named '", splines_by, "' in data used for model calibration", sep = ""))
    }
    fact_by <- factor(rep(levels(dsm_model$data[, splines_by])[1], nx),
                      levels = levels(dsm_model$data[, splines_by])
                      )
    Xnew <- cbind(fact_by, Xnew)
    names(Xnew)[1] <- splines_by
  }

  ### handles random effects
  if(!is.null(random)) {
    if(!all(random %in% names(dsm_model$data))) {
      stop("Check random effect: no corresponding match in data used for model calibration")
    } else {
      rand <- as.data.frame(array(NA, dim = c(nx, length(random))))
      names(rand) <- random
      for(j in random) {
        rand[, j] <- factor(rep("new_level", nx), levels = levels(dsm_model$data[, j]))
      }
      Xnew <- cbind(Xnew, rand)
    }
  }
  ### approximate posterior distribution with MV normal
  beta <- mvtnorm::rmvnorm(n_sim, mean = dsm_model$coefficients, sigma = dsm_model$Vp)

  if(remove_intercept) { beta[, 1] <- 0; writeLines("\tRemoving intercept") }

  # rm_spatial <- grep("X,Y", names(dsm_model$coefficients), fixed = TRUE)
  rm_spatial <- grep("longitude,latitude", names(dsm_model$coefficients), fixed = TRUE)
  # check for soap
  soap <- length(grep(pattern = "bnd", x = as.character(formula(dsm_model))[3], fixed = TRUE)) != 0
  if(length(rm_spatial) != 0) {
    if(soap) {
      gridxy <- expand.grid(longitude = seq(min(segdata[, "longitude"], na.rm = TRUE), max(segdata[, "longitude"], na.rm = TRUE), length.out = 1e2),
                            latitude = seq(min(segdata[, "latitude"], na.rm = TRUE), max(segdata[, "latitude"], na.rm = TRUE), length.out = 1e2)
                            )
    } else {
      gridxy <- expand.grid(longitude = (seq(min(segdata[, "longitude"], na.rm = TRUE), max(segdata[, "longitude"], na.rm = TRUE), length.out = 1e2) - mean(segdata[, "longitude"], na.rm = TRUE)) / sd(segdata[, "longitude"], na.rm = TRUE),
                            latitude = (seq(min(segdata[, "latitude"], na.rm = TRUE), max(segdata[, "latitude"], na.rm = TRUE), length.out = 1e2) - mean(segdata[, "latitude"], na.rm = TRUE)) / sd(segdata[, "latitude"], na.rm = TRUE)
                            )
    }
    Z <- cbind(gridxy,
               as.data.frame(sapply(var_name[-which(var_name %in% c("longitude", "latitude"))], function(id) { rep(0, 1e4) }))
               )

    if(!is.null(splines_by)) {
      Z <- cbind(fact_by, Z)
      names(Z)[1] <- splines_by
    }

    if(!is.null(random)) {
      Z <- cbind(Z, rand)
    }

    Z <- predict(dsm_model, newdata = Z, off.set = 1, type = "lpmatrix")

    if(!is.null(splines_by)) {
      Z[, grep(splines_by, names(Z), fixed = TRUE)] <- 0
    }
    slinpred <- beta %*% t(Z); rm(Z)
    # gridxy <- expand.grid(X = seq(min(segdata[, "X"], na.rm = TRUE), max(segdata[, "X"], na.rm = TRUE), length.out = 1e2),
    #                       Y = seq(min(segdata[, "Y"], na.rm = TRUE), max(segdata[, "Y"], na.rm = TRUE), length.out = 1e2)
    #                       )
    gridxy <- expand.grid(longitude = seq(min(segdata[, "longitude"], na.rm = TRUE), max(segdata[, "longitude"], na.rm = TRUE), length.out = 1e2),
                          latitude = seq(min(segdata[, "latitude"], na.rm = TRUE), max(segdata[, "latitude"], na.rm = TRUE), length.out = 1e2)
                          )
    gridxy <- rbind(gridxy, gridxy)
    gridxy$spatial <- c(apply(slinpred, 2, median), apply(exp(slinpred), 2, median))
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
    if(!is.null(splines_by)) {
      for(k in 1:length(levels(dsm_model$data[, splines_by]))) {
        Z <- Xnew
        Z[, j] <- (seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx) - mean(segdata[, j], na.rm = TRUE)) / sd(segdata[, j], na.rm = TRUE)
        Z[, splines_by] <- factor(rep(levels(dsm_model$data[, splines_by])[k], nx),
                                  levels = levels(dsm_model$data[, splines_by])
                                  )

        Z <- predict(dsm_model, newdata = Z, off.set = 1, type = "lpmatrix")
        Z[, grep("longitude,latitude", names(as.data.frame(Z)), fixed = TRUE)] <- 0
        linpred <- beta %*% t(Z)
        df_splines <- rbind(df_splines,
                            data.frame(x = seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx),
                                       y = apply(linpred, 2, median),
                                       lower = apply(linpred, 2, lower),
                                       upper = apply(linpred, 2, upper),
                                       param = rep(j, nx),
                                       scale = rep("log", nx),
                                       level = rep(levels(dsm_model$data[, splines_by])[k], nx)
                                       ),
                            data.frame(x = seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx),
                                       y = apply(exp(linpred), 2, median),
                                       lower = apply(exp(linpred), 2, lower),
                                       upper = apply(exp(linpred), 2, upper),
                                       param = rep(j, nx),
                                       scale = rep("natural", nx),
                                       level = rep(levels(dsm_model$data[, splines_by])[k], nx)
                                       )
                            )
      }
    } else {
      Z <- Xnew
      Z[, j] <- (seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx) - mean(segdata[, j], na.rm = TRUE)) / sd(segdata[, j], na.rm = TRUE)
      Z <- predict(dsm_model, newdata = Z, off.set = 1, type = "lpmatrix")
      Z[, grep("longitude,latitude", names(as.data.frame(Z)), fixed = TRUE)] <- 0
      linpred <- beta %*% t(Z); rm(Z)
      df_splines <- rbind(df_splines,
                          data.frame(x = seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx),
                                     y = apply(linpred, 2, median),
                                     lower = apply(linpred, 2, lower),
                                     upper = apply(linpred, 2, upper),
                                     param = rep(j, nx),
                                     scale = rep("log", nx)
                                     ),
                          data.frame(x = seq(min(segdata[, j], na.rm = TRUE), max(segdata[, j], na.rm = TRUE), length.out = nx),
                                     y = apply(exp(linpred), 2, median),
                                     lower = apply(exp(linpred), 2, lower),
                                     upper = apply(exp(linpred), 2, upper),
                                     param = rep(j, nx),
                                     scale = rep("natural", nx)
                                     )
                          )
    }
  }

  # graph part : g_splines
  g_fun <- function(df, theme, transformation) {

    g <- ggplot(data = df_splines %>%
                  filter(scale == transformation & param == unique(df_splines$param)[p]),
                aes(x = x, y = y, ymin = lower, ymax = upper)
    ) +
      geom_ribbon(alpha = 0.3, fill = "midnightblue") +
      geom_line(color = "midnightblue") +
      facet_grid(level ~ param, scales = "free_x") +
      scale_y_continuous(name = "Density") +
      scale_x_continuous(name = "Covariate") +
      theme_bw() +
      theme

    return(g)
  }

  # log and natural
  if(!is.null(splines_by)) {

    list_g_param_log <- list()
    list_g_param_natural <- list()

    len_param <- length(unique(df_splines$param))

    for(p in 1:len_param) {

      # custom theme to create only one faceting in y axis
      if(p > 1 & p < len_param){
        theme_custom <- theme(
          strip.background.y  = element_blank(),
          strip.text.y = element_blank(),
          axis.title.y = element_blank()
        )
      }
      if(p == 1){
        theme_custom <- theme(
          strip.background.y  = element_blank(),
          strip.text.y = element_blank(),
        )
      }
      if(p == len_param){
        theme_custom <- theme(
          axis.title.y = element_blank()
        )
      }

      g_param_log <- NULL
      g_param_natural <- NULL

      g_param_log <- g_fun(df = df_splines,
                       theme = theme_custom,
                       transformation = "log")

      g_param_natural <- g_fun(df = df_splines,
                               theme = theme_custom,
                               transformation = "natural")

      list_g_param_log[[p]] <- g_param_log
      list_g_param_natural[[p]] <- g_param_natural
    }

    g1 <- do.call(arrangeGrob,
                  c(list_g_param_natural,
                    ncol = len_param,
                    nrow = 1,
                    top = "natural"))
    g1 <- ggdraw() + draw_grob(grobTree(g1))

    g2 <- do.call(arrangeGrob,
                  c(list_g_param_log,
                    ncol = len_param,
                    nrow = 1,
                    top = "log"))
    g2 <- ggdraw() + draw_grob(grobTree(g2))


    g_splines <- list(g1, g2)

  } else {

    g_splines <- ggplot(data = df_splines,
                        aes(x = x, y = y, ymin = lower, ymax = upper)
                        ) +
      geom_ribbon(alpha = 0.3, fill = "midnightblue") +
      geom_line(color = "midnightblue") +
      facet_grid(scale ~ param, scales = "free") +
      scale_y_continuous(name = "Density") +
      scale_x_continuous(name = "Covariate") +
      theme_bw()

  }

  return(list(df_splines = df_splines,
              g_splines = g_splines,
              spatial = gridxy
              )
         )
}

#' @import ggplot2
#' @importFrom MASS rnegbin
#' @importFrom mvtnorm rmvnorm
#' @export

rootogram_nb <- function(model_fit, n_rep = 1e3, min_obs = 0, max_obs = 1e3, by = NULL) {
  ### posterior predictive checks
  beta <- mvtnorm::rmvnorm(n_rep, mean = model_fit$coefficients, sigma = model_fit$Vp)
  Z <- predict(model_fit, type = "lpmatrix")
  linpred <- beta %*% t(Z)
  transfo_overdispersion <- exp(model_fit$family$getTheta())
  y_rep <- t(apply(exp(linpred), 1,
                   function(x) { MASS::rnegbin(n = length(x), mu = x, theta = transfo_overdispersion) }
                   )
             )
  ### check for NA
  if(any(is.na(y_rep))){
    rm_row <- apply(y_rep, 1, function(row) {any(is.na(row))})
    y_rep <- y_rep[-rm_row, ]
  }
  y_rep <- ifelse(y_rep > max_obs, max_obs, y_rep)

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

  ### make df
  rootdf <- function(index) {
    data.frame(mids = (min_obs:max_obs + (min_obs + 1):(max_obs + 1))/2,
               y_obs = as.numeric(f_histogram(x = countdata[index], max_obs = max_obs)),
               y_rep = apply(apply(y_rep[, index], 1, f_histogram, max_obs = max_obs), 1, mean)
               )
  }

  ### by
  index <- list(1:length(countdata))
  out <- do.call('rbind', lapply(index, rootdf))
  if(!is.null(by)) {
    if(!any(names(model_fit$data) == by)) {
      writeLines("\t* Ignoring 'by' argument as no matching column name in data")
    }
    else {
      stratum <- model_fit$data[, by]
      index <- lapply(unique(stratum), function(x) { which(stratum == x) })
      out <- do.call('rbind', lapply(index, rootdf))
      out$by <- rep(unique(stratum), each = length((min_obs:max_obs + (min_obs + 1):(max_obs + 1))/2))
    }
  }

  ### goodness of fit
  gof <- do.call('cbind', lapply(index, function(x) {apply(y_rep[, x], 1, cantonnier, y = countdata, breaks = min_obs:max_obs, remove_zeroes = TRUE)}))

  ###
  ### plot
  theme_set(theme_bw(base_size = 16))
  g <- out %>%
    ggplot(aes(x = mids, y = y_rep)) +
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
  if(any(names(out) == "by")) {
    g <- g +
      facet_wrap(~ by)
  }
  ### wrap-up
  return(list(rootogram = g,
              earthMover = gof
              )
         )
}
