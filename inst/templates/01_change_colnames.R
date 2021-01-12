source("R/00_setup.R")


# Create mirror repository in "res/" ----------------------------------------
path_save <- "res/01_change_colnames"
if(!dir.exists(path_save)) {
  dir.create(path_save, recursive = T)
}


# Get effort and observation in "data/" -----------------------------------
effort <- read.dbf("data/effort/##########.dbf")
observation <- read.dbf("data/observation/##########.dbf")
  


# Prepare effort and observation  -----------------------------------------

# delete all columns with NA and useless 
observation <-  observation %>% 
  remove_empty(which = "cols") %>% 
  select(!taxon_eng:cd_taxsup)


# Change colnames of effort and observation with change_varName() ---------
standard_effort <- change_effort_varName(effort)
standard_obs <- change_obs_varName(observation)



# Verify if all observation are in effort data ----------------------------
verif_effort_obs(c("legId"), standard_obs, standard_effort)


# Save standard_effort and standard_obs in "res/" -------------------------
write.dbf(standard_effort, paste(path_save,"standard_effort",sep="/"))
write.dbf(standard_obs, paste(path_save,"standard_obs",sep="/"))

rm(effort, observation)