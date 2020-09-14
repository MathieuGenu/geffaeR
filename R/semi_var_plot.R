#' @import ggplot2
#' @export

# faire un plot de la semi-variance (semivarPlot)
semi_var_plot <- function(variogram, vario_model, distance) {

    semivario_function <- function(h, sill, range, form) {
      if(form == "matern") { f <- sill *(1 -(1 + h * sqrt(3) / range) * exp(- h * sqrt(3) / range)) }
      if(form == "exponential") { f <- sill *(1 - exp(- h / range)) }
      return(f)
    }


  obs_df <- with(variogram, data.frame(x = u, y = v, n = n))
  pred_df <- data.frame(x = distance,
                        y = vario_model$nugget +
                          semivario_function(distance,
                                             vario_model$cov.pars[1],
                                             vario_model$cov.pars[2],
                                             form = vario_model$cov.model
                                             )
                        )


  theme_set(theme_bw(base_size = 14))

  g <- ggplot() +
    geom_point(data = obs_df,
               aes(x = x, y = y, size = log1p(n)), color = "red") +
    geom_line(data = pred_df,
              aes(x = x, y = y),
              color = "midnightblue", size = 1) +
    xlab("Distance") +
    ylab("Semi-variance") +
    guides(size = "none") +
    theme(plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 12)
    )

  return(g)

}
