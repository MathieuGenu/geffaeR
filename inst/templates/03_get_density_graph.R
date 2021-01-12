
# Source needed packages and function -------------------------------------
source("R/00_setup.R")

# Create mirror repository in "res/" --------------------------------------
path_save <- "res/CDS/03_get_density_graph"
if(!(dir.exists(path_save))) {
  dir.create(path_save, recursive = T)
}



# Get necessary data ------------------------------------------------------

# get csv tab of distance compilation (summary_all_sp.csv)
df_sum_ds <- read.table("res/CDS/02_get_detection_info_and_graph/summary_all_sp_sep_semicolon.csv",
                        sep=";", h=T, dec=",")

# Load species name and list
path_res_species <- "res/02_prepare_effort_and_obs/"
load(paste(path_res_species,"species_info.RData",sep="/"))

# add session and year column
df_sum_ds <- df_sum_ds %>% 
  mutate(session = str_extract(Region, '[:digit:]$'),
         Année = ifelse(session > 4, "2020","2019")) %>% 
  mutate(session = factor(session, levels = c("1","5","2","6","3","7","4","8")))


# graph
# x_label <- c("1 (hiver)","2 (printemps)","3 (ete)","4 (automne)")
# x_label <- unique(df_sum_ds$Region)
# x_label <- c("1","2","3","4","5","6","7","8")
for (s in (1:length(sp_name))) {
  
  if (sp_name[s] != "TRASH") {
    y_label <- "Animal density (ind.km²)"
  } else {
    y_label <- "density (item.km²)"
  }
  
  select_sp_df <- df_sum_ds %>% 
    filter(species == sp_name[s])
  
  ggplot(data = select_sp_df , aes(x = session, y = animal_density)) +
    geom_bar(stat="identity", color = "black", aes(fill = Année)) +
    geom_errorbar(data = select_sp_df, aes(ymin = animal_density-animal_density*cv_animal_density,
                                           ymax = animal_density+animal_density*cv_animal_density),
                  width = 0.5) +
    theme_minimal() +
    theme(axis.text = element_text(size = 10)) +
    scale_x_discrete(xlab("Session")) + 
    scale_y_continuous(ylab(y_label))
  ggsave(filename = paste(path_save,"/density_for_",sp_name[s],".png",sep=""),
         width = 8, height = 4.5, dpi = 300)
    
}


