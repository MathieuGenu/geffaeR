#' Preparation of observation data
#'
#' This function transform raw observation data into multiple sub data.frame for next analysis of
#' other functions of the package.
#'
#' @param sp Name of taxon, group, family or species.
#' @param obs_base Observation database.
#' @param legdata "Legdata" data.frame, built with \code{\link{prepare_data_effort}}.
#' @param segdata "Segdata" data.frame, built with \code{\link{prepare_data_effort}}.
#' @inheritParams prepare_data_effort
#' @param projection Projection of \code{shape} object.
#' @param truncation Default = \code{NULL}. disstance of truncation of the sampling (in km)
#' @param remove_sp Species to remove of the data set.
#' @param unit_km Default = \code{FALSE}. Is the unit of distance in the data set is in km ?
#' @return This function return a list containing :
#'         \enumerate{
#'           \item distdata : Data.frame with distance information, which aim to be an input for
#'           \code{\link[Distance]{ds}}.
#'           \item obsdata : Data.frame containing information at observation scale.
#'           \item countdata_leg : A data.frame that merge leg scale effort informations
#'           with number of sightings and number of individuals (N and Y).
#'           \item countdata_seg : A data.frame that merge segment scale effort informations
#'           with number of sightings and number of individuals (N and Y).
#'         }
#' @examples
#'
#'
#' @export

prepare_data_obs <- function(sp, obs_base, legdata, segdata, shape, shape_layer, projection,
                             truncation = NULL, remove_sp = NULL,
                             unit_km = FALSE) {

  # polygons sampling
  if(any("character" %in% is(shape))){
    poly_NC <- readOGR(dsn = paste(shape), layer = shape_layer, verbose = F) # NC
  } else {
    poly_NC <- shape
  }


  # load Observation base
  raw_obs <- obs_base
  raw_obs <- subset(raw_obs, segId %in% unique(segdata$Seg))

  # remove center observation
  raw_obs <- subset(raw_obs, side != "CENTER")

  # # ne prendre que l'espece choisie
  # if(group) {
  #   sp_data <- subset(raw_obs, group %in% sp)
  # } else if(taxon) {
  #   sp_data <- subset(raw_obs, taxon %in% sp)
  # } else if(family) {
  #   sp_data <- subset(raw_obs, family %in% sp)
  # } else {
  #   sp_data <- subset(raw_obs, species %in% sp)
  # }

  # ne prendre que l'espece choisie ###TEST###
  # cas pour plusieurs espèces en même temps

  sp_data <- raw_obs[0,]

  if (any(sp %in% unique(raw_obs$group))) {
    match_group <- sp[which(sp %in% unique(raw_obs$group))]
    sp_data_group <- subset(raw_obs, group %in% match_group)
    sp_data <- rbind(sp_data, sp_data_group)
  }
  if (any(sp %in% unique(raw_obs$taxon))) {
    match_taxon <- sp[which(sp %in% unique(raw_obs$taxon))]
    sp_data_taxon <- subset(raw_obs, taxon %in% match_taxon)
    sp_data <- rbind(sp_data, sp_data_taxon)
  }
  if (any(sp %in% unique(raw_obs$family))) {
    match_family <- sp[which(sp %in% unique(raw_obs$family))]
    sp_data_family <- subset(raw_obs, family %in% match_family)
    sp_data <- rbind(sp_data, sp_data_family)
  }
  if (any(sp %in% unique(raw_obs$species))) {
    match_species <- sp[which(sp %in% unique(raw_obs$species))]
    sp_data_species <- subset(raw_obs, species %in% match_species)
    sp_data <- rbind(sp_data, sp_data_species)
  }


  # qu'est ce ?
  if(!is.null(remove_sp)) { sp_data <- subset(sp_data, species != remove_sp) }

  # ne garder que les obs dans la bande pour les oiseaux
  if(all(sp_data$taxon == "Oiseau marin")) {
    cli_alert_info("Keeping only observations within the 200 m around the transect")
    sp_data <- subset(sp_data, decAngle %in% c(1, 3))
  }

  ## distance en km
  if(!unit_km){
    sp_data$distance <- sp_data$perpDist/1000
  } else {
    sp_data$distance <- sp_data$perpDist
  }

  # troncation
  if(!is.null(truncation)) {
    wa <- truncation
  } else {
    #tronque a 5%, arrondie a la classe superieure
    if(unit_km == F){
      corr_km <- 1/1000
    } else {
      corr_km <- 1
    }
    pas <- 50 * corr_km
    wa <- as.numeric(with(sp_data, quantile(distance,  prob = 0.95)))
    wa <- pas*floor(wa/pas) + ifelse(wa %% pas == 0, 0, pas)
  }

  sp_data <- subset(sp_data, distance <= wa)

  ## countdata seg et leg pour rajouter les observation à legdata et segdata avec la fonction ajout_obs()
  if("session" %in% colnames(sp_data)) {
    countdata_seg <- as.data.frame(sp_data[, c("transect", "legId", "segId","podSize","session")] %>%
                                     group_by(transect, legId, segId) %>%
                                     summarize(n = n(), count = sum(podSize)))
  } else {
    countdata_seg <- as.data.frame(sp_data[, c("transect", "legId", "segId","podSize")] %>%
                                     group_by(transect, legId, segId) %>%
                                     summarize(n = n(), count = sum(podSize)))
  }

  colnames(countdata_seg) <- c("Transect.Label", "Sample.Label", "Seg","n", "y")

  countdata_leg <- as.data.frame(countdata_seg %>%
    group_by(Transect.Label, Sample.Label) %>%
    summarize(n_detected = sum(n), n_ind = sum(y)))

  # creation des tableaux necessaires a l'analyse ----
  #--------------------------------------------------#

  ## table : distdata with truncation
  if("session" %in% colnames(sp_data)) {
    distdata <- sp_data[, c("transect", "strate", "legId", "segId", "lon", "lat", "lon",
                            "lat", "podSize", "distance", "observerId","session")]
    distdata$strate <- paste(sp_data$subRegion, sp_data$strate)
    names(distdata) <- c("Transect.Label", "Region.Label", "Sample.Label", "Seg", "X", "Y",
                         "longitude", "latitude", "size", "distance","observerId","session")
  } else {
    distdata <- sp_data[, c("transect", "strate", "legId", "segId", "lon", "lat", "lon",
                            "lat", "podSize", "distance", "observerId")]
    distdata$strate <- paste(sp_data$subRegion, sp_data$strate)
    names(distdata) <- c("Transect.Label", "Region.Label", "Sample.Label", "Seg", "X", "Y",
                         "longitude", "latitude", "size", "distance","observerId")
  }


  distdata$object <- as.numeric(row.names(distdata))
  distdata$detected <- 1
  distdata <- subset(distdata, Sample.Label %in% legdata$Sample.Label)

  # re-projeter les obs dans le système de projection indiqué dans projection
  distdata_xy <- distdata[, c("longitude", "latitude")]
  coordinates(distdata_xy) <- ~ longitude + latitude
  distdata_xy@proj4string <- CRS(projection)
  distdata_xy <- spTransform(distdata_xy, CRS(as.character(poly_NC@proj4string)))
  distdata[, c("X", "Y")] <- distdata_xy@coords

  if("session" %in% colnames(sp_data)) {
    distdata <- left_join(dplyr::select(legdata, -survey, -left, -right),
                          dplyr::select(distdata, -Transect.Label, -Region.Label,-session),
                          by="Sample.Label")
  } else {
    distdata <- left_join(dplyr::select(legdata, -survey, -left, -right),
                          dplyr::select(distdata, -Transect.Label, -Region.Label),
                          by="Sample.Label")
  }


  distdata$detected[is.na(distdata$detected)] <- 0 # à verifier

  ### table : obsdata
  if("session" %in% colnames(sp_data)) {
    obsdata <- subset(distdata, detected == 1)[, c("distance", "size", "Transect.Label", "Region.Label",
                                                   "Seg", "Sample.Label","observerId","session")]
  } else {
    obsdata <- subset(distdata, detected == 1)[, c("distance", "size", "Transect.Label", "Region.Label",
                                                   "Seg", "Sample.Label","observerId")]
  }
  distObject <- distdata$object
  obsdata$object <- distObject[!is.na(distObject)]

  ## rassembler
  return(list(distdata = distdata,
              obsdata = obsdata,
              countdata_leg = countdata_leg,
              countdata_seg = countdata_seg,
              trunc = wa))
}
