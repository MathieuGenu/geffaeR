#' \encoding{Modelisation Baésienne à partir d'un modèle stan.}
#'
#' \encoding{Cette fonction permet de d'appliquer un modèle bayésien sur les données
#' d'observation (distadata et segdata_obs) à partir d'un modèle stan}
#'
#' @param cds_data \encoding{Modèle écrit en "stan". Le modèle doit être renseigné sous forme de character string.}
#' @param segdata_obs \encoding{data.frame segdata contenant les infos par segment avec le nombre de d'individus et le
#'        nombre de détections. Provient de la fonction \code{\link{ajout_obs}}}
#' @param distdata \encoding{data.frame distdata. Provient de la fonction \code{\link{prepare_data_obs}}}
#' @return
#'
#' @examples
#'
#'
#' @export

stan_cds <- function(cds_data, segdata_obs, distdata) {

  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  stan_cds_hn <- stan_model(model_code = cds_data,
                            model_name = "quasi-Conventional Distance Sampling (Negative Binomial likelihood)"
  )

  ### initialize algorithm
  staninit_cds_hn = function(n_strata, n_obs) {
    return(list(unscaled_intercept = rnorm(n_strata),
                aux_intercept = rgamma(n_strata, 1.5, rate = 1.5),
                unscaled_alpha = rnorm(1),
                aux_alpha = rgamma(1, 1.5, rate = 1.5),
                unscaled_sigma2_esw = rgamma(1, 2.0, rate = 1.0),
                tau_esw = rgamma(1, 2.0, 1.0),
                logit_omega = rnorm(1)
                )
          )
  }


  standata_cds <- function(segdata_obs_, distdata_, truncation_in_km = 1.0) {
    datalist <- list(n_obs = nrow(segdata_obs),
                     n_strata = length(unique(segdata_obs$session)),
                     OBS = segdata_obs$y,
                     STRATUM = as.numeric(as.factor(segdata_obs$session)),
                     EFFORT = segdata_obs$Effort,
                     W = segdata_obs$seaState + 1,
                     n_beaufort = max(c(distdata$seaState, segdata_obs$seaState)) + 1,
                     n_detection = nrow(subset(distdata, detected == 1)),
                     TRUNCATION = truncation_in_km,
                     BEAUFORT = subset(distdata, detected == 1)$seaState + 1,
                     DISTANCE = subset(distdata, detected == 1)$distance,
                     prior_scale_distance = log(2)
                     )
    return(datalist)
  }

  n_chains <- 4
  n_iter <- 750
  n_warm <- 250
  n_thin <- 1
  mod <- sampling(obj = stan_cds_hn,
                         data = standata_cds(segdata_obs_ = segdata_obs, distdata_ = distdata),
                         pars = c('intercept', 'alpha', 'inv_omega', 'esw', 'log_lik'),
                         chains = n_chains,
                         iter = n_iter,
                         warmup = n_warm,
                         thin = n_thin)

  chain_plot <- plot(mod,
                     pars = c('intercept', 'alpha', 'inv_omega', 'esw'),
                     plotfun = "trace",
                     inc_warmup = TRUE)


  # pairs_plot <- pairs(mod, pars = c('intercept', 'alpha', 'inv_omega', 'esw', 'density'))

  # print(mod, digits = 3)

  intercept <- summary(exp(rstan::extract(mod, 'intercept')$intercept))

  esw <- summary(rstan::extract(mod, 'esw')$esw)

  # summary(1/rstan::extract(mod, 'inv_omega')$inv_omega)
  return(list(
    model = mod,
    chain_plot = chain_plot,
    intercept = intercept,
    esw = esw))

}
