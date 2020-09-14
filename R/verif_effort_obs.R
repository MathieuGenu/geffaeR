#' Compatibility verification between effort and observation data.
#'
#' Matching verification of columns values of prepared effort and prepared observation.
#'
#' @param var Variables for which matching between effort and observation data have to be verified.
#' @param standard_obs Prepared observation data. e.g. Output of \code{\link[geffaeR]{prepare_data_obs}}.
#' @param standard_effort Prepared effort data. e.g. Output of \code{\link[geffaeR]{prepare_data_effort}}.
#' @import cli
#' @rawNamespace import(crayon, except = `%+%`)
#' @export

verif_effort_obs <- function(var, standard_obs, standard_effort) {

  # error messages
  error <- red $ bold
  warn <- magenta
  note <- cyan

  NA_match <- function(var) {

    nb_error <- 0

    if(any(var %in% "session") &
       !all("session" %in% colnames(standard_obs) , "session" %in% colnames(standard_effort) )) {

      cli_li(error('There are no column session in standard_obs nor standard_effort'))

    } else {

      # NA in effort
      if(any(is.na(standard_effort[,var]))) {
        cli_li(error('There are NA(s) in standard_effort${var}'))
        nb_error <- nb_error + 1
      }

      # NA in obs
      if(any(is.na(standard_obs[,var]))) {
        cli_li(error('There are NA(s) in standard_obs${var}'))
        nb_error <- nb_error + 1
      }

      # mathch obs & effort
      if(any(!(unique(standard_obs[,var]) %in% unique(standard_effort[,var])))) {
        cli_li(error("Some {var} of standard_obs don't match {var} of standard_effort"))
        nb_error <- nb_error + 1
      }

    }

    return(nb_error)

  }

  are_there_error <- lapply(var, NA_match)
  recap_error <- do.call('sum', are_there_error)
  if (recap_error < 1) {
    cli_alert_success("All variable tested are complete with no NA and observation matches effort")
  }

}

