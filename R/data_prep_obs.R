#' \encoding{Préparation des données d'observation.}
#'
#' \encoding{Cette fonction prépare les donnees d'observation pour qu'elles soient compatibles
#' pour les autres fonctions de ce package.}
#'
#' @param species \encoding{Nom du groupe, de la famille, du taxon d'interêt.}
#' @param obs_base \encoding{nom de la base de données contenant les observations de \code{species}.}
#' @param DataDir \encoding{Direction du répertoire où se trouvent le shape file et la base d'observation
#' Cela suppose donc que le shape file et la base sont dans le meme fichier.}
#' @param legdata \encoding{data.frame "Legdata" au prealable construit avec la fonction
#'        \code{\link{prepare_data_effort}}.}
#' @param segdata \encoding{data.frame "Segdata" au prealable construit avec la fonction
#'        \code{\link{prepare_data_effort}}.}
#' @inheritParams prepare_data_effort
#' @return Cette fonction renvoie une liste contenant :
#'         \enumerate{
#'           \item distdata : Un data.frame contenant les infos sur chaque .....
#'           \item obsdata : Un data.frame contenant les infos sur chaque observation.
#'           \item \encoding{countdata_leg : Un data.frame contenant le nombre de detection
#'           et d'individus total pour les legs
#'           où il y a eu observation pour l'espece ou group en question.}
#'           \item \encoding{countdata_seg : Un data.frame contenant le nombre de detection
#'           et d'individus total pour les segment
#'           ou il y a eu observation pour l'espece ou group en question.}
#'         }
#' @examples
#'
#'
#' @export

prepare_data_obs <- function(sp, obs_base, DataDir, legdata, segdata, shape, shape_layer, projection = projection,
                             group = FALSE, family = FALSE, taxon = FALSE, truncation = NULL, remove_sp = NULL,
                             optimal = TRUE,bird = FALSE, unit_km = FALSE) {

  # polygons sampling
  poly_NC <- readOGR(dsn = paste(DataDir, shape, sep = "/"), layer = shape_layer, verbose = F)


  # load Observation base
  raw_obs <- obs_base
  raw_obs <- subset(raw_obs, segId %in% unique(segdata$Seg))

  # ne prendre que l'espece choisie
  if(group) {
    sp_data <- subset(raw_obs, group_ %in% sp)
    } else if(taxon) {
      sp_data <- subset(raw_obs, taxon %in% sp)
    } else if(family) {
      sp_data <- subset(raw_obs, family %in% sp)
    } else {
        sp_data <- subset(raw_obs, species %in% sp)
    }
  # qu'est ce ?
  if(!is.null(remove_sp)) { sp_data <- subset(sp_data, species != remove_sp) }

  # ne garder que les obs dans la bande pour les oiseaux
  if(all(sp_data$taxon == "Oiseau marin")) {
    writeLines("Keeping only observations within the 200 m around the transect")
    sp_data <- subset(sp_data, decAngle %in% c(1, 3))
  }

  ## distance en km
  if(!unit_km){
    sp_data$distance <- sp_data$PerpDist/1000
  } else {
    sp_data$distance <- sp_data$PerpDist
  }

  # troncation
  if(!is.null(truncation)) {
    wa <- truncation
  } else {
    #tronque a 5%, arrondie a la classe superieure
    pas <- 50
    wa <- as.numeric(with(sp_data, quantile(PerpDist,  prob = 0.95)))
    wa <- pas*floor(wa/pas) + ifelse(wa%%pas == 0, 0, pas)
  }

  sp_data <- subset(sp_data, distance <= wa)

  ## countdata seg et leg pour rajouter les observation à legdata et segdata avec la fonction ajout_obs()
  if("session" %in% colnames(sp_data)) {
    countdata_seg <- as.data.frame(sp_data[, c("transect", "IdLeg", "segId","podSize","session")] %>%
                                     group_by(transect, IdLeg, segId) %>%
                                     summarize(n = n(), count = sum(podSize)))
  } else {
    countdata_seg <- as.data.frame(sp_data[, c("transect", "IdLeg", "segId","podSize")] %>%
                                     group_by(transect, IdLeg, segId) %>%
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
    distdata <- sp_data[, c("transect", "strate", "IdLeg", "segId", "lon", "lat", "lon",
                            "lat", "podSize", "distance", "observerId","session")]
    distdata$strate <- paste(sp_data$subRegion, sp_data$strate)
    names(distdata) <- c("Transect.Label", "Region.Label", "Sample.Label", "Seg", "X", "Y",
                         "longitude", "latitude", "size", "distance","observerId","session")
  } else {
    distdata <- sp_data[, c("transect", "strate", "IdLeg", "segId", "lon", "lat", "lon",
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
    distdata <- left_join(dplyr::select(legdata, -survey, -left_, -right_, -session),
                          dplyr::select(distdata, -Transect.Label, -Region.Label),
                          by="Sample.Label")
  } else {
    distdata <- left_join(dplyr::select(legdata, -survey, -left_, -right_),
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
              countdata_seg = countdata_seg))
}
