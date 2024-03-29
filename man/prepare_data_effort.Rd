% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_prep_effort.R
\name{prepare_data_effort}
\alias{prepare_data_effort}
\title{Preparation of effort data.}
\usage{
prepare_data_effort(
  effort_base,
  covariable = NULL,
  block_area,
  shape,
  shape_layer,
  New_projection,
  optimal = FALSE,
  col2keep = NULL
)
}
\arguments{
\item{effort_base}{\code{data.frame} with effort data.}

\item{covariable}{\code{vector} of covariable names to keep in output of the function.}

\item{block_area}{Data.frame with 2 colnames :
\enumerate{
  \item Block.
  \item Area.
}}

\item{shape}{Shapefile of the study area. It can be either a \code{SpatialPolygonsDataFrame} class object,
in this case it is not necessary to give shape_layer argument. Or it can be the name of the shape
object with its extension ".shp" (ex : "data/studyAreaShapefile.shp").}

\item{shape_layer}{Layer of the shapefile if shape is a \code{character string}.}

\item{New_projection}{New projection of longitude and latitude of columns
(POINT_X et POINT_Y) in \code{Proj4String} format,
see : \code{\link[sp]{CRS}} for more infos.}

\item{optimal}{Argument which allows to keep data sampled in optimal conditions.
Defaults settings are all data are kept. In case of optimal = T, indexes
\code{"c("GG", "GM", "MG", "EG", "GE", "EE", "ME", "EM", "MM")"} are kept.}

\item{col2keep}{Columns to keep from effort_base in output of the function.}
}
\value{
This function return a list containing:
        \enumerate{
          \item legdata : \code{data.frame} with infos at leg scale.
          \item segdata : \code{data.frame} with infos at segment scale.
        }
}
\description{
Transform raw effort data into multiple sub data.frame for next analysis of
other functions of the package.
}
\examples{

}
