# Source needed packages --------------------------------------------------
source('R/00_setup.R')
set.seed(666)

# Create mirror repository in "res/" ----------------------------------------
path_save <- "res/CDS/01_adjust_plot_detection"
if(!dir.exists(path_save)) {
  dir.create(path_save, recursive = T)
}


# get sp, effort and observation prepared ---------------------------------
path_res <- "res/02_prepare_effort_and_obs"
load(file.path(path_res,"effort_output.RData"))
load(file.path(path_res,"list_prepare_obs_by_sp.RData"))
load(file.path(path_res,"species_info.RData"))

# need to load standard_obs
standard_obs <- read.dbf("res/01_change_colnames/standard_obs.dbf")



# Adjust detection plot with plot_detection() -----------------------------

list_effort_w_obs <- list()
list_plot_detection_by_session <- list()


for(s in 1:length(sp)) {
  
  # hacking method to have cds estimation by session
  list_prepare_obs_by_sp[[s]]$distdata$Region.Label <- paste(list_prepare_obs_by_sp[[s]]$distdata$Region.Label,
                                                             list_prepare_obs_by_sp[[s]]$distdata$session,
                                                             sep = "_")
  list_prepare_obs_by_sp[[s]]$obsdata$Region.Label <- paste(list_prepare_obs_by_sp[[s]]$obsdata$Region.Label,
                                                            list_prepare_obs_by_sp[[s]]$obsdata$session, 
                                                            sep = "_")
  
  # segregate seabirds
  obs_base_sp <- standard_obs %>% 
    select(taxon, group, family, species) %>% 
    filter_all(any_vars(. %in% sp[[s]]))
  
  # strip-transect for non marine mammal
  if ((unique(obs_base_sp$taxon) %in% c("Seabird","Oiseau marin","Coastal Bird","Oiseau cotier",
                                        "Land Bird","Oiseau terrestre")) | 
      unique(obs_base_sp$species %in% c("JELLY","TRASH","Dechet","BUOY"))) {
    exclude <- TRUE
  } else {
    exclude <- FALSE
  }
  
  if (exclude) {
    list_prepare_obs_by_sp[[s]]$distdata <- list_prepare_obs_by_sp[[s]]$distdata %>% 
      filter(distance <= 0.2 | is.na(distance))
  }
  
  
  # Add observation on effort
  effort_w_obs <- ajout_obs(countdata_leg = list_prepare_obs_by_sp[[s]]$countdata_leg,
                            countdata_seg = list_prepare_obs_by_sp[[s]]$countdata_seg,
                            legdata = effort_output$legdata,
                            segdata = effort_output$segdata)
  list_effort_w_obs[[s]] <- effort_w_obs
  names(list_effort_w_obs)[s] <- paste("effort_w_obs",sp_name[s],sep="")
  
  
  # isolate distdata
  distdata <- list_prepare_obs_by_sp[[s]]$distdata
  
  # get parameter depending on exclude
  list_g_det <- list()
  
  if(!exclude){
    bin = seq(0.0,list_prepare_obs_by_sp[[s]]$trunc, 0.05)
    key = "halfnorm"
    upper = list_prepare_obs_by_sp[[s]]$trunc
  } else {
    bin = NULL
    key = NULL
    upper = NULL
  }
  
  # adjust function
  g_det_all <- plot_detection(distdata = distdata,
                              bin = bin,
                              key = key,
                              upper = upper,
                              is_seabird = exclude)
  
  if(!exclude){
    g_det_all$graph <- g_det_all$graph + 
      ggtitle(paste(sp[[s]],"all session")) 
  }
  
  
  list_g_det[[1]] <- g_det_all
  
  names(list_g_det)[1] <- "all_session"
  
  list_plot_detection_by_session[[s]] <- list_g_det
  names(list_plot_detection_by_session)[s] <- paste("plot_detection",sp_name[s],sep="")
}



# Save effort_w_obs and list_plot_detection_by_session --------------------
rm(list=setdiff(ls(), c("list_plot_detection_by_session","path_save", "sp_name","list_effort_w_obs")))
save("list_plot_detection_by_session", 
     file = file.path(path_save,"list_plot_detection_by_session.RData"))
save("list_effort_w_obs", file=file.path(path_save,"list_effort_w_obs.RData"))
rm(list_plot_detection_by_session,path_save, list_effort_w_obs, sp_name)



