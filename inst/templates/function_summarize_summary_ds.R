# summary of distance function (package Distance)
sum_ds_table <- function(sum_ds) {
  
  nrow_sum <- nrow(sum_ds$dht$clusters$summary)
  
  # cluster summary
  clust_part <- sum_ds$dht$clusters$summary[,c(1,2,3,4,5,7,9)]
  if (nrow_sum >1) {
    clust_part <- clust_part[-nrow_sum,]
  }
  colnames(clust_part) <- c("Region", "area_Surface", "surface_Sample", "effort_Sampled", "sightings", "Encounter_rate", "cv_Encounter_rate")
  clust_part[,c(3,4,6,7)] <- round(clust_part[,c(3,4,6,7)], 6)
  
  # expected
  expected_part <- sum_ds$dht$Expected.S
  expected_part$cv.Expected.S <- expected_part$se.Expected.S/expected_part$Expected.S
  expected_part <- expected_part[-nrow_sum, -c(1,3)]
  colnames(expected_part) <- c("mean_group", "cv_mean_group")
  expected_part <- round(expected_part, 6)
  
  # abondance
  abond_part <- sum_ds$dht$individuals$N[,c(2,5,6)]
  if (nrow_sum >1) {
    abond_part <- abond_part[-nrow_sum,]
  }
  colnames(abond_part) <- c("animal_abundance", "Abon_min_(IC95%)", "Abon_max_(IC95%)")
  abond_part <- round(abond_part)
  
  # density
  density_part <- sum_ds$dht$individuals$D[,c(2,4)]
  if (nrow_sum >1) {
    density_part <- density_part[-nrow_sum,]
  }
  colnames(density_part) <- c("animal_density", "cv_animal_density")
  density_part <- round(density_part, 6)
  
  # colbind all
  summary_ds <- cbind(clust_part, expected_part, density_part, abond_part)
  
  return(summary_ds)
  
}

