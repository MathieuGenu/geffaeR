% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ajout_obs.R
\name{ajout_obs}
\alias{ajout_obs}
\title{Adding observation on effort.}
\usage{
ajout_obs(countdata_leg, legdata, countdata_seg, segdata)
}
\arguments{
\item{countdata_leg}{data.frame containing the detection number and the total number of
individuals for legs for which observation have beebn made.}

\item{legdata}{"Legdata" data.frame, built with \code{\link{prepare_data_effort}}.}

\item{countdata_seg}{data.frame containing the detection number and the total number of
individuals for segments for which observation have beebn made.}

\item{segdata}{"Segdata" data.frame, built with \code{\link{prepare_data_effort}}.}
}
\value{
This function return a list containing :
        \enumerate{
          \item legdata_obs : data.frame corresponding to legdata with observation information.
          \item segdata_obs : data.frame corresponding to segdata with observation information.
        }
}
\description{
Adding the number of individuals and the group size on prepared effort data
(e.g. on legdata and segdata).
}
\examples{


}
