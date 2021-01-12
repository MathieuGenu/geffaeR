
# Source needed packages and function -------------------------------------
source("R/00_setup.R")
source("R/CDS/function/function_summarize_summary_ds.R")

# Create mirror repository in "res/" --------------------------------------
path_save <- "res/CDS/02_get_detection_info_and_graph"
if(!dir.exists(path_save)) {
  dir.create(path_save, recursive = T)
}


# Get necessary data ------------------------------------------------------

# load species name and list
path_res <- "res/CDS/01_adjust_plot_detection/"
load(file.path("res/02_prepare_effort_and_obs/species_info.RData"))

# load dectection plot adjustment
load(file.path(path_res,"list_plot_detection_by_session.RData"))



# Initiate esw table ------------------------------------------------------
tab_esw_cv <- data.frame(matrix(rep(NA,5), ncol=5))
colnames(tab_esw_cv) <- c("species","esw","esw_cv_min","esw_cv_max","CV_%")


for(s in (1:length(list_plot_detection_by_session))) {
  
  # stock esw and cv
  if(!is.null(list_plot_detection_by_session[[s]]$all_session$esw)) {
    temp_esw_tab <- c(sp_name[s],
                      list_plot_detection_by_session[[s]]$all_session$esw,
                      list_plot_detection_by_session[[s]]$all_session$esw_cv)
    if(s==1) {
      tab_esw_cv[1,] <- temp_esw_tab
    } else {
      tab_esw_cv <- rbind(tab_esw_cv, temp_esw_tab)
      row.names(tab_esw_cv) <- NULL
    }
  }
  
  ds_model <- list_plot_detection_by_session[[s]]$all_session$distFit
  tab <- sum_ds_table(ds_model)
  species <- rep(sp_name[s],nrow(tab))
  tab <- cbind(species, tab)
  
  if(s == 1) {
    summary_all_sp <- tab
  } else {
    summary_all_sp <- rbind(summary_all_sp, tab)
  }
  
  # export plot detection graph
  if(!is.null(list_plot_detection_by_session[[s]]$all_session$esw)) {
    
    list_plot_detection_by_session[[s]]$all_session$graph +
      scale_x_continuous(breaks = seq(0,1,0.1))
    
    ggsave(filename = paste(path_save,paste("detection_plot_for",sp_name[s],".png",sep="_"),sep="/"),
           width = 8, height = 4.5, dpi = 300)
  }
  
}


# Save tables -------------------------------------------------------------

write.table(summary_all_sp, paste(path_save,"summary_all_sp.csv",sep="/"),
            sep=",", dec=".", col.names = T, row.names = F)
write.table(summary_all_sp, paste(path_save,"summary_all_sp_sep_semicolon.csv",sep="/"),
            sep=";",dec=",", col.names = T, row.names = F)
write.table(tab_esw_cv, paste(path_save,"esw_cv_sep_semicolon.csv",sep="/"),
            sep=";",dec=",", col.names = T, row.names = F)

rm(list=ls())
