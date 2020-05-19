#' @export

get_density <- function(distdata = NULL,
                        truncation = NULL,
                        legdata,
                        esw = NULL,
                        esw_cv = NULL,
                        g_0 = NULL,
                        g_0_cv = NULL,
                        by_leg = TRUE,
                        jackknife = FALSE,
                        alpha = 0.05
) {
  ### must provide the legdata
  if(is.null(legdata)) {
    stop("Must provide a legdata dataframe")
  }

  ### some sanity checks
  if(is.null(distdata) && is.null(esw)) {
    stop("Must provide either a distdata dataframe or an esw")
  }

  if(is.null(esw) && !is.null(distdata)) {

    if(is.null(truncation)) {

      stop("Must provide a truncation value for distance data")

    } else {

      writeLines("\t Half-normal key for detection function")
      mod <- crch::crch(distance ~ -1 | 1, data = subset(distdata, distance <= truncation),
                        link.scale = "log", dist = "gaussian", left = 0, right = truncation,
                        truncated = TRUE, control = crch::crch.control(maxit = 1e4)
      )
      sigma <- exp(as.numeric(coef(mod)) + sqrt(as.numeric(vcov(mod))) * rnorm(1e6))
      esw <- (pnorm(truncation, 0, sigma) - 0.5) / dnorm(0, 0, sigma)
      esw_cv <- sd(esw) / mean(esw)
      esw <- mean(esw)

    }
  }

  ### estimate density
  estimate_density <- function(leg) {

    stat <- leg %>%
      summarize(N = sum(n_detected),
                L = sum(Effort),
                Es = mean(n_ind / n_detected, na.rm = TRUE)
      ) %>%
      mutate(esw = esw,
             g_0 = ifelse(is.null(g_0), 1, g_0),
             D_bar = N * Es / (2 * esw * L * g_0),
             n_cv = sqrt(L / (nrow(leg) - 1) * sum(leg$Effort * (leg$n_detected / leg$Effort - N / L)^2)) / N,
             Es_cv = sd(leg$n_ind / leg$n_detected, na.rm = TRUE) / Es,
             esw_cv = ifelse(is.null(esw_cv), 0, esw_cv),
             g_0_cv = ifelse(is.null(g_0_cv), 0, g_0_cv)
      ) %>%
      as.data.frame()
    return(stat)
  }
  ### point estimate
  if(by_leg) {

    sumstat <- estimate_density(leg = legdata)

  } else {

    legdata <- legdata %>%
      group_by(Transect.Label) %>%
      summarize(Effort = sum(Effort),
                n_detected = sum(n_detected),
                n_ind = sum(n_ind)
      ) %>%
      as.data.frame()
    sumstat <- estimate_density(leg = legdata)

  }

  if(jackknife) {

    if(by_leg) {

      jack_leg <- lapply(unique(legdata$Sample.Label), function(id) { subset(legdata, Sample.Label != id) })
      jack_leg <- do.call('rbind', lapply(jack_leg, estimate_density)) %>%
        mutate(Sample.Label = legdata$Sample.Label) %>%
        left_join(legdata %>% select(Sample.Label, Effort),
                  by = "Sample.Label"
        )

    } else {

      jack_leg <- lapply(unique(legdata$Transect.Label), function(id) { subset(legdata, Transect.Label != id) })
      jack_leg <- do.call('rbind', lapply(jack_leg, estimate_density)) %>%
        mutate(Transect.Label = legdata$Transect.Label) %>%
        left_join(legdata %>% select(Transect.Label, Effort),
                  by = "Transect.Label"
        )

    }
    jack_leg <- jack_leg %>%
      mutate(dj = sumstat$L * (sumstat$D_bar - D_bar) / Effort + D_bar,
             Dj = sum(jack_leg$Effort * jack_leg$dj) / sumstat$L,
             delta = Effort * (dj - Dj)^2
      ) %>%
      as.data.frame()

    sumstat <- sumstat %>%
      mutate(D_bar = sum(jack_leg$Effort * jack_leg$dj) / L,
             D_cv = sqrt((sum(jack_leg$delta) / (L * (nrow(jack_leg) - 1))) / D_bar^2 + Es_cv^2 + esw_cv^2 + g_0_cv^2),
             D_lower = D_bar * exp(qt(alpha / 2, df = nrow(legdata) - 1) * sqrt(log1p(D_cv^2))),
             D_upper = D_bar * exp(qt(1 - alpha / 2, df = nrow(legdata) - 1) * sqrt(log1p(D_cv^2)))
      ) %>%
      as.data.frame()

  } else {

    # (see 1st edition distance sampling page 89-90)
    sumstat <- sumstat %>%
      mutate(D_cv = sqrt(n_cv^2 + Es_cv^2 + esw_cv^2 + g_0_cv^2),
             D_df = D_cv^4 / (n_cv^4 / (nrow(legdata) - 1) + esw_cv^4 / N),
             D_lower = D_bar * exp(qt(alpha / 2, df = D_df) * sqrt(log1p(D_cv^2))),
             D_upper = D_bar * exp(qt(1 - alpha / 2, df = D_df) * sqrt(log1p(D_cv^2)))
      ) %>%
      as.data.frame()

  }

  return(sumstat)
}

