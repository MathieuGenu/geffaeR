#' distance smapling analysis
#'
#' This function calculate the ESW (effective-strip witdth) of a taxonomic group,
#' It also use the \link[Distance]{ds} function to estimate all the abondance, density encounter-rate,...
#' parameter and summarize it in the results of the function.
#' Finally i gives the detection plot of the distance analysis.
#'
#' @param distdata distadata data.frame previously built with  \code{\link{prepare_data_obs}}.
#' @param bin Min and max value of the distance detection and the step between
#' theese two.
#' @param key Choice bewtween two key function :
#'            \enumerate{
#'              \item "halfnorm" : Half-Normal.
#'              \item "hazard" : hazard-rate.
#'            }
#'
#' @param upper Maximum value of the distance range of detection estimation.
#' @return This function return list with :
#'         \enumerate{
#'           \item esw : effective strip width. with confidence interval at 95 percent.
#'           \item esw_cv : coefficient of variation of esw in percent.
#'           \item graph : Barplot of detections depending of distance, with confindence
#'           estimated with mcmc chains.
#'         }
#' @examples
#'
#'
#' @export


plot_detection <- function(distdata, bin, key, upper = NULL) {
  # renommer observer car capote la fonction ds() (nom attribué à une fonction pour ds)
  if(nrow(distdata) < 2) {
    stop(paste("il doit y avoir au moins 2 lignes d'observation dans distdata"))
  }

  if(!is.null(distdata$observer)) {
    distdata$observer <- distdata$observer_ID
  }

  #distance sampling
  fit = ds(distdata, key = "hn", adjustment = NULL)

  if(is.null(upper)) {
    upper <- max(bin)
  }

  ### Compute the effective strip width of the transect
  hn <- function(x, sigma) {
    return(exp(-x^2/(2 * sigma^2)))
  }
  hz <- function(x, sigma, nu) {
    return(1 - exp(-(x/sigma)^(-nu)))
  }
  get_summary <- function(x){
    c(mean(x), HPDinterval(as.mcmc(x), prob = 0.95))
  }

  x <- seq(min(bin), max(bin), length.out = 1000)

  ### approximation to the posterior distribution
  my_coef <- as.numeric(fit$ddf$par)
  my_mat <- as.matrix(solve(fit$ddf$hessian))

  if(key == "halfnorm") {
    simpleMC <- exp(rnorm(1000, my_coef, sqrt(as.numeric(my_mat))))
    esw <- sapply(simpleMC, function(q) {
      integrate(hn, 0, upper, sigma=q)$value
    })
    y <- sapply(x, function(i) {
      sapply(simpleMC, hn, x = i)
    })
  }

  if(key == "hazard") {
    simpleMC <- exp(rmvnorm(1000, my_coef, my_mat))
    esw <- sapply(1:1000,function(k){
                    integrate(hz, 0, upper, sigma = simpleMC[k, 2], nu = simpleMC[k, 1])$value
    })
    y <- sapply(x, function(i) {
      sapply(1:1000, function(k) {
        hz(x = i, sigma = simpleMC[k, 2], nu = simpleMC[k, 1])
      })
    })
  }

  y <- as.data.frame(cbind(x, t(apply(y, 2, get_summary))))
  names(y) <- c("distance", "moyenne", "binf", "bsup")

  # histogramme des distances
  dd <- hist(distdata$distance, breaks = bin, plot = FALSE)
  step <- diff(bin)
  area_hist <- sum(step*dd$density)

  drect <- data.frame(x1 = numeric(length(step)), x2 = numeric(length(step)),
                      y1 = numeric(length(step)), y2 = numeric(length(step)))

  for(i in 2:length(dd$breaks)){
    drect[i-1, ] <- c(dd$breaks[i-1], dd$breaks[i], 0, dd$density[i-1]/area_hist*mean(esw))
  }

  g_plot_detection <- ggplot() +
    geom_rect(data = drect, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
              colour = "black", fill = "white") +
    geom_vline(xintercept = mean(esw), col ="red", linetype = "dotted", size = 1.5) +
    geom_line(data = y, mapping = aes(x = distance, y = moyenne), size = 1, col = "midnightblue") +
    geom_line(data = y, mapping = aes(x = distance, y = binf), linetype = "dashed", col = "midnightblue",
              size = 1, alpha = 0.6) +
    geom_line(data = y, mapping = aes(x = distance, y = bsup), linetype = "dashed", col = "midnightblue",
              size = 1, alpha = 0.6) +
    scale_x_continuous(breaks = bin, limits = range(bin)) + xlab("Perpendicular Distance") +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) + ylab("Pr(Detection)") +
    theme(plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 10))

  esw = round(get_summary(esw), 3)
  esw_cv = 100 * round(sd(esw)/mean(esw), 3)

  return(list(graph = g_plot_detection,
              esw = esw,
              esw_cv = esw_cv,
              distFit = fit))
}
