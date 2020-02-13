# WARNING #
# this package must not be run completly but in case by case
# all function for package editing are here but must be run in different optic
# ex : document is run when we want to refresh documentation of functions
# except for the fisrt part, it have to be run in prelude


# FIRST PART : ----
# script for package writing
library(devtools)
library(roxygen2)
library(usethis)
library(available)
library(rstantools)

use_devtools()


# check if the name we want to assign to the new package is available
# available("geffaeR")

# just when we have to create a package
# create_package("~/Documents/Projet/package/geffaeR")

setwd("C:/Users/mgenu.RATUFA/Documents/Projet/package/geffaeR")


# END OF FISRT PART ----

### Modifier le fichier DESCRIPTION du package ###
#------------------------------------------------#
use_description(fields = list(`Authors@R` = NULL))
use_gpl3_license("Mathieu Genu & Matthieu Authier")
use_description(fields = list(Author = "Mathieu Genu [aut, cre],
                              Matthieu Authier [aut, cre]",
                              Title = "Generic function for abundance estimation with R",
                              Language = "fr",
                              Description = "The package transforms SAMMOA output and applies different analysis and give graphical output."))
use_description(fields = list(
                              Language = "fr"
))
# import of functions in the description field of the package
packg_to_import <- c("rgdal", "maptools",  "sp", "raster","foreign", "lubridate", "reshape",
                     "dplyr", "ggplot2", "cowplot", "mvtnorm", "Distance", "captioner",
                     "Rdistance", "dsm", "knitr", "fields", "ggthemes", "broom", "coda",
                     "purrr","VIM","Amelia","missMDA","FactoMineR", "rstan", "Rcpp",
                     "geoR", "MASS","WhatIf","arm","viridisLite","stats","utils",
                     "mgcv","cli","crayon")
for(pckg in packg_to_import){
  use_package(pckg, type="Imports")
}



# check package (regarder les erreurs et oublis dans documentation)
check()

# mettre ? jour la documentation des fonction (les fichiers .Rd)
document()

# load tout le package (tout les fields)
devtools::load_all("./", quiet = TRUE)

# modify the version of the package
use_version(which = "dev") # which = "major"/"minor"/"patch"/"dev"

# export package in .tar.gz
build()

# Vignette creation
usethis::use_vignette("Package_example")

# mettre des data dans le package
load(paste("C:/Users/mgenu.RATUFA/Documents/Projet/donnees/data_package/Europe_fond2carte.RData"))
usethis::use_data(cds)
usethis::use_data(cds_krig)
usethis::use_data(Europe)
usethis::use_data(NEA)
usethis::use_data(NEA_isobath100)
usethis::use_data(NEA_isobath200)
usethis::use_data(lbrt93_proj)

# # pour rstan
# use_rstan()


source("C:/Users/mgenu.RATUFA/Documents/Projet/Build_example_for_package/res/01_build_study_area_with_effort_obs/")

# create a shapefile for vignette example
load("C:/Users/mgenu.RATUFA/Documents/Projet/Build_example_for_package/res/01_build_study_area_with_effort_obs/effort_observation.RData")
usethis::use_data(spdf)

# create observation table
usethis::use_data(observation)

# create effort table
usethis::use_data(effort)
