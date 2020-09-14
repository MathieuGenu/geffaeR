#' @import dplyr
#' @import ggplot2
#' @importFrom arm bayesglm
#' @importFrom fields rdist
#' @importFrom reshape melt
#' @export

variogram_intraTransect <- function(segdata_obs,
                                    obs,
                                    esw,
                                    intraTransect = T,
                                    breaks = NULL,
                                    plot = TRUE,
                                    pairs.min = 100,
                                    n_sim = 100,
                                    distmat = NULL,
                                    region = NULL,
                                    esp = "species"  # only for title
) {

  # y = donnÃ©es de comptage
  # esw = effective half-strip width en km
  # l = longueur en km des legs ou segments en km
  # center = matrice de dim(length(y), 2) avec les coordonnÃ©es XY des centres de chaque segment
  # id_transect = numero de code du transect (pour en tenir sinon laisser null)
  # breaks = break points for empirical variogram
  # distmat = distance matrix
  col_needed <- c("y","n","X","Y","Transect.Label","Region.Label","Effort")
  if(!all(col_needed %in% colnames(segdata_obs))) {
    missing_col <- col_needed[!(col_needed %in% colnames(segdata_obs))]
    stop(paste("Les colonnes suivantes ne sont pas dans les colonnes de segdata_obs :",
               paste(missing_col, collapse = "\n"), sep="\n"))
  }
  if(obs == "n_ind") {
    y <- segdata_obs$y
  }
  if(obs == "n_detection") {
    y <- segdata_obs$n
  }
  l <- segdata_obs$Effort

  center <- segdata_obs[, c("X", "Y")]
  id_transect <- segdata_obs$Transect.Label
  region.label <- segdata_obs$Region.Label

  if(!is.null(region)){
    data=cbind(y,center, id_transect,region.label) %>% filter(region.label==region)
    y=data[,1]
    center=data[,2:3]
    id_transect=data[,4]
    rm(data)

  }
  if(is.null(breaks)) {
    writeLines("Automatic, but potentially crappy, choice for distance bins of 2.5 km, capped at 100 km")
    breaks <- c(-0.01, 0.001, 1, seq(2.5, 100, 2.5))
  }

  ### fonction pour calculer le taux de rencontre moyen
  TauxObs <- function(y, esw, l){
    my_control <- list(epsilon = 1e-6, maxit = 100000, trace = FALSE)
    # utiliser quasi poisson pour eventuelle surdispersion
    model_m <- arm::bayesglm(y ~ 1 + offset(2 * esw * l),
                             prior.mean.for.intercept = -1.0, # prior is less than 4 obs/ind per km
                             prior.scale.for.intercept = 1.0,
                             prior.df.for.intercept = 7,
                             family = "quasipoisson", control = my_control
                             )
    m <- mean(exp(rnorm(1000, as.numeric(coef(model_m)), as.numeric(sqrt(vcov(model_m))))))
    return(round(m, 6))
  }

  ### calculer la matrice de distance
  if(is.null(distmat)) {
    D <- fields::rdist(center) / 1e3 # in kilometers
    # tenir compte du design de l'Ã©chantillonnage
    if(!intraTransect) { id_transect <- rep(1, length(y)) }
    D_design <- fields::rdist(as.numeric(id_transect)) * 1000
    distmat <- D + D_design
  }
  distclas <- cut(distmat, breaks = breaks)

  # calcul de la moyenne du taux d'observation m pour correction Poisson du Variogramme
  m <- TauxObs(y, esw, l)

  # calcul du variogramme
  calc_vario <- function(y, l, esw, m) {
    return(0.5 *((outer(l, l, "*") * 2 * esw / outer(l, l, "+")) *(outer(y /(l * 2 * esw), y /(l * 2 * esw), "-"))^2 - m))
  }
  z2 <- calc_vario(y = y, l = l, esw = esw, m = m)

  # constante de pondÃ©ration
  cstpond <- tapply(outer(l, l, "*") * 2 * esw / outer(l, l, "+"),
                    distclas,
                    sum, na.rm = TRUE
  )

  #on divise par la constante de pondÃ©ration
  vario1 <- tapply(z2, distclas, sum, na.rm = TRUE) / cstpond
  distvario1 <- tapply(distmat, distclas, mean, na.rm = TRUE)
  nbcouples1 <- tapply(z2, distclas, function(x) { length(which(!is.na(x))) })

  # enlever les estimations trop bruitÃ©es, ie pas assez de couples de points
  noise <- c(which(is.na(nbcouples1)), which(nbcouples1 < pairs.min))
  noise <- as.numeric(noise[order(noise)])
  vario1 <- ifelse(vario1 < 0, 0, vario1)

  if(length(noise) != 0) {
    vario1 <- vario1[-noise]
    distvario1 <- distvario1[-noise]
    nbcouples1 <- nbcouples1[-noise]
  }

  # mettre les rÃ©sultats dans un objet vario de classe "variogram", package geoR
  vario <- list(u = as.numeric(distvario1), # distance
                v = as.numeric(vario1),     # estimated variogram values
                n = as.numeric(nbcouples1), # nombre de couples de points
                sd = rep(0, length(vario1)),
                bins.lim = breaks[-c(1, noise)],
                ind.bin = ifelse(as.numeric(nbcouples1) > pairs.min, TRUE, FALSE),
                var.mark = var(y /(2 * esw * l)),
                beta.ols = m,
                output.type = "bin",
                max.dist = max(breaks),
                estimator.type = "classical", # it should rather be Monestiezal
                n.data = length(y),
                lambda = 1,
                trend = "cte",
                pairs.min = pairs.min,
                nugget.tolerance = 1e-12,
                direction = "omnidirectional",
                tolerance = "none",
                uvec = as.numeric(distvario1),
                call = NULL
  )
  class(vario) <- "variogram"

  # test graphique
  if(plot) {
    # calculer une enveloppe nulle
    # Quelques fonctions simples pour faire du bootstrap et dÃ©terminer s'il faut utiliser une loi nÃ©gative binomiale
    overdispersion <- function(x) { return(var (x) / mean (x))}
    permute_index <- function(x) { return(x[sample(c(1:length(x)), length(x), replace = TRUE)])}
    bootstrap_pval <- function(x, simulator, statistic, n_sim, alpha, truth) {
      tboots <- replicate(n_sim, statistic(simulator (x)))
      ci_lower <- 2*statistic(x) - quantile(tboots, 1 - alpha/2, na.rm=TRUE)
      ci_upper <- 2*statistic(x) - quantile(tboots, alpha/2, na.rm=TRUE)
      return(as.numeric(ifelse(ci_lower > truth, 0, ifelse(ci_upper < truth, 0, 1))))
    }

    if(overdispersion(y) < 1) { negbin <- FALSE }
    else { negbin <- ifelse(bootstrap_pval(y, permute_index, overdispersion, 10000, 0.10, 1) == 0, TRUE, FALSE) }
    null_env <- function(distri) {
      # simuler des donnÃ©es
      if(distri) {
        x <- rnbinom(length(y), size = mean(y)/(overdispersion(y) - 1), mu = mean(y))
        model_x <- glm(x ~ 1 + offset(2 * esw * l), family = "quasipoisson", control = list(epsilon = 1e-6, maxit = 10000, trace = FALSE))
      }
      else{
        x <- rpois(length(y), mean(y))
        model_x <- glm(x ~ 1 + offset(2 * esw * l), family = "poisson", control = list(epsilon = 1e-6, maxit = 10000, trace = FALSE))
      }
      z0 <- calc_vario(y = x, l = l, esw = esw,
                       m = mean(exp(rnorm(10000, as.numeric(coef(model_x)), as.numeric(sqrt(vcov(model_x))))))
      )
      vario0 <- tapply(z0, distclas, sum, na.rm = TRUE) / cstpond
      if(length(noise) != 0) { vario0 <- vario0[-noise] }
      return(vario0)
    }
    null <- t(replicate(n_sim, null_env(distri = negbin), simplify = "array"))

    # resumer les simulations et rajouter les donnÃ©es sur la derniÃ¨re ligne
    sumnull <- rbind(apply(null, 2, quantile, probs = c(0.1, 0.5, 0.9), na.rm = TRUE),
                     as.numeric(vario1)
    )
    sumnull= sumnull %>% t() %>%data.frame() %>%  mutate(distcl=as.character(as.numeric(distvario1)))
    colnames(sumnull)<- c("lower", "median", "upper", "emp","distcl")

    # reprendre les simulations
    null <- as.data.frame(null)
    names(null) <- as.character(as.numeric(distvario1))
    null$sim <- as.character(1:n_sim)
    null <- melt(null, id = "sim")
    null$variable <- as.numeric(distvario1)[as.numeric(null$variable)]

    # graphique avec ggplot2
    theme_set(theme_bw(base_size = 14))
    g <- ggplot() +
      geom_line(data = null,
                aes(x = variable, y = value, group = sim),
                alpha = 0.05
      ) +
      geom_ribbon(data = sumnull,
                  aes(x = as.numeric(distcl), ymin = lower, ymax = upper),
                  fill = "midnightblue",
                  alpha = 0.3, linetype = "dashed", size = 1
      ) +
      geom_line(data = sumnull,
                aes(x = as.numeric(distcl), y = emp),
                color = "firebrick1", size = 1
      ) +
      labs(x="Distance",y="Semi-variance", title=paste0(esp," in ", region)) +
      theme(plot.title = element_text(lineheight = 0.8, face = "bold"),
            axis.text = element_text(size = 12)
      ) +coord_cartesian(ylim=c(0, max(sumnull$emp)))

    return(list(variogram = vario,
                g = g,
                distmat = distmat))
  } else {
    return(list(variogram = vario,
                distmat = distmat))
  }
}
