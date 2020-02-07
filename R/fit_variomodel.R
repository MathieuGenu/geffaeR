
#' @export

# fitter le model "matern" ou "exponential", qui sert pour les autres fonctions
fit_variomodel <- function(variogram, form = c("matern", "exponential")) {
  vario_model <- geoR::variofit(variogram,
                                cov.model = form,
                                fix.nugget = FALSE,
                                nugget = mean(variogram$v) / 4,
                                fix.kappa = TRUE,
                                kappa = ifelse(form == "matern", 1.5, 1),
                                weights = "npairs",
                                minimisation.function = "optim",
                                ini.cov.pars = expand.grid(
                                  seq(0, mean(variogram$v)*2, l = 10),
                                  seq(as.numeric(quantile(variogram$u, prob = 0.1)),
                                      as.numeric(quantile(variogram$u, prob = 0.8)),l = 10)
                                ),
                                messages = F)
  return(vario_model = vario_model)
}
