
# Source needed packages --------------------------------------------------
source("R/00_setup.R")

# Create mirror repository in "res/" ----------------------------------------
path_save <- "res/03_prepare_effort_obs"
if(!dir.exists(path_save)) {
  dir.create(path_save, recursive = T)
}


# Get standard_effort and standard_obs ------------------------------------
standard_effort <- read.dbf("res/02_change_colnames/standard_effort.dbf")
standard_obs <- read.dbf("res/02_change_colnames/standard_obs.dbf")



# Prepare species parameters ----------------------------------------------

# It's possible to put multiple taxon at a time in "sp" but give it one name in "sp_name"
sp <- list("Petit delphinine",
           "PHOPHO",
           "TURTRU",
           "ALCURI", 
           "SULBAS",
           "LARMIN",
           "LARRID",
           "Goeland noir",
           "Goeland gris",
           "Sternidae",
           "Procellariidae",
           "Hydrobatidae",
           "RISTRI", 
           "CATSKU", 
           c("SMAGUL","LARMEL","LARRID","LARMIN"),
           c("LARRID","LARMEL"),
           "LARGUL",
           "JELLY",
           "TRASH",
           c("Delphininae", "Delphinidae"),
           "Requin",
           c("THUSPP", "LARFIS", "XIPGLA"),
           "MOLMOL",
           "Dechet",
           "BUOY",
           "Mouette",
           "SMASHE"
           )

sp_name <- c("Petit delphinine",
             "PHOPHO", 
             "TURTRU",
             "ALCURI",
             "SULBAS", 
             "LARMIN", 
             "LARRID",
             "Goeland noir",
             "Goeland gris",
             "Sternidae",
             "Procellariidae",
             "Hydrobatidae",
             "RISTRI",
             "CATSKU", 
             "grp_smagul_larmel_larmin_larrid",
             "grp_larrid_larmel",
             "LARGUL",
             "JELLY",
             "TRASH",
             "delphini_nae_dae",
             "Requin",
             "THUSPP_LARFIS_XIPGLA",
             "MOLMOL",
             "Dechet",
             "BUOY",
             "Mouette",
             "Petits_puffins")

trunc <-  c(0.5,
            0.4,
            0.6,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            0.5,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA,
            NA)

# save species
save(sp, sp_name, file = paste(path_save, "species_info.RData", sep="/"))




# Prepare effort data -----------------------------------------------------

covariable <- NULL

# projection of shapefile
shape <- "data/shape/Shape_study.shp"
study_area <- sf::st_read(shape)


# projection for future analysis
New_projection <-  lbrt93_proj # in package lbrt93_proj saved

# assciation of region with its surface
Block_Area <- data.frame(Block = c("ATL_N1"), # col strate_sec de standard_effort
                         Area = c(14947)) # area of each element of strate_sc of standard_effort


effort_output <- prepare_data_effort(effort_base = standard_effort, 
                                     shape = shape,
                                     optimal = T,
                                     block_area = Block_Area, 
                                     New_projection = New_projection,
                                     covariable = covariable)

legdata <- effort_output$legdata
segdata <- effort_output$segdata



# Prepare observation data ------------------------------------------------

list_prepare_obs_by_sp <- list()
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

for (s in 1:length(sp)) {
  
  if(!is.na(trunc[s])) {
    truncation = trunc[s]
  } else {
    truncation <- NULL
  }
  
  
  observation_output <- prepare_data_obs(sp = sp[[s]], 
                                         obs_base = standard_obs, 
                                         truncation = truncation,
                                         legdata, 
                                         segdata, 
                                         shape, 
                                         projection = projection)
  
  list_prepare_obs_by_sp[[s]] <- observation_output
  names(list_prepare_obs_by_sp)[s] <- paste(sp_name[[s]],"obs_output",sep="_")
  
}


# Save list of prepared observation and effort ----------------------------

save(list_prepare_obs_by_sp, file = paste(path_save, "list_prepare_obs_by_sp.RData", sep="/"))
save(effort_output, file = paste(path_save, "effort_output.RData", sep="/"))

rm(list=ls())