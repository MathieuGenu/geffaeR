% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verif_effort_obs.R
\name{verif_effort_obs}
\alias{verif_effort_obs}
\title{Compatibility verification between effort and observation data.}
\usage{
verif_effort_obs(var, standard_obs, standard_effort)
}
\arguments{
\item{var}{Variables for which matching between effort and observation data have to be verified.}

\item{standard_obs}{Prepared observation data. e.g. Output of \code{\link[geffaeR]{prepare_data_obs}}.}

\item{standard_effort}{Prepared effort data. e.g. Output of \code{\link[geffaeR]{prepare_data_effort}}.}
}
\description{
Matching verification of columns values of prepared effort and prepared observation.
}
