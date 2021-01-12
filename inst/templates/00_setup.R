# clean environment
rm(list = ls())


#CRAN packages vector
cran_packages <- c(
  "devtools",
  "foreign",
  "rgdal",
  "rlist",
  "tidyverse",
  "magrittr",
  "ggspatial",
  "sf",
  "rgeos",
  "cli",
  "crayon",
  "janitor",
  "stringr",
  "lubridate"
)

github_packages <- c(
  "mathUtil",
  "geffaeR"
)

github_repository <- c(
  "MathieuGenu/mathUtils",
  "MathieuGenu/geffaeR"
)

# c_p_n_i : cran packages not installed
c_p_n_i <- cran_packages[!(cran_packages %in% installed.packages())]
g_p_n_i <- which(!(github_packages %in% installed.packages()))


# installation of packages
lapply(c_p_n_i, install.packages, dependencies = TRUE)
lapply(github_repository[g_p_n_i], devtools::install_github , dependencies = TRUE)

#install packages
lapply(c(cran_packages,github_packages), function(x){
  library(x, character.only = TRUE, quietly = TRUE)
})

rm(c_p_n_i, g_p_n_i, cran_packages, github_packages, github_repository)
