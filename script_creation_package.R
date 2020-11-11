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
library(pkgdown)

use_devtools()


# check if the name we want to assign to the new package is available
# available("geffaeR")

# just when we have to create a package
# create_package("~/Documents/Projet/package/geffaeR")

# setwd("C:/Users/mgenu.RATUFA/Documents/Projet/package/geffaeR")


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
                     "geoR", "MASS","arm","viridisLite","stats","utils",
                     "mgcv","cli","crayon","crch","sf","ggspatial","loo","tweedie", "glue",
                     "DT", "htmltools", # pour la vignette
                     "methods", "gridExtra" # suggested by package checking
                     )
for(pckg in packg_to_import){
  use_package(pckg, type="Imports")
}

use_package("WhatIf", min_version = T)
use_git_remote(name = "WhatIf",
               url = "https://github.com/IQSS/WhatIf.git",
               overwrite = T)


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
devtools::build_vignettes()
devtools::build()

# mettre des data dans le package
load(paste("C:/Users/mgenu.RATUFA/Documents/Projet/donnees/data_package/Europe_fond2carte.RData"))
load(paste("C:/Users/mgenu.RATUFA/Documents/Projet/donnees/data_package/NEA.RData"))

# version simplifié NEA
library(magrittr)
library(sf)
library(rmapshaper)
sf_study_area <- NEA %>%
  st_as_sf(crs=4326) %>%
  st_crop(xmin = -10,
          xmax = 15,
          ymin = 35,
          ymax = 55) %>%
  st_union() %>%
  as_Spatial() %>%
  ms_simplify() %>%
  st_as_sf()

library(ggplot2)
sf_study_area %>%
  ggplot() +
  geom_sf()

NEA_simplified_FR <- sf_study_area %>%
  as_Spatial()

usethis::use_data(cds)
usethis::use_data(cds_krig)
usethis::use_data(Europe)
usethis::use_data(NEA, overwrite = T)
usethis::use_data(NEA_isobath100)
usethis::use_data(NEA_isobath200)
usethis::use_data(lbrt93_proj)
usethis::use_data(NEA_simplified_FR)



# Data part ---------------------------------------------------------------

source("C:/Users/mgenu.RATUFA/Documents/Projet/Build_example_for_package/res/01_build_study_area_with_effort_obs/")

# create a shapefile for vignette example
load("C:/Users/mgenu.RATUFA/Documents/Projet/simul_ppp/res/effort_observation.RData")
effort_example <- effort
observation_example <- observation
shape_example <- spdf
usethis::use_data(shape_example, overwrite = T)

# create observation table
usethis::use_data(observation_example, overwrite = T)

# create effort table
effort_example$CenterTime <- lubridate::date(effort_example$CenterTime)
usethis::use_data(effort_example, overwrite = T)

# create MOLMOL data list for DSM
load("C:/Users/mgenu.RATUFA/Documents/Projet/MOLMOL_data_vignette_geffaeR/res/DSM_pack_MOLMOL/DSM_pack_MOLMOL.RData")
usethis::use_data(DSM_pack_MOLMOL, overwrite = T)


# pkgdown ----------------------------------------------------------------

# run one at the beginning
# usethis::use_pkgdown()
pkgdown::build_site()
pkgdown::build_home()
pkgdown::build_reference()
pkgdown::build_article(name = "vignettes/Package_example.Rmd")
usethis::use_github_action("pkgdown")


# add files in Rbuildignore
#--------------------------
usethis::use_build_ignore(c("index.Rmd","index.md","logo.png","readme.md","read","script_creation_package.R"))
