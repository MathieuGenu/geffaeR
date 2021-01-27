#' Preparation of effort data.
#'
#' Transform raw effort data into multiple sub data.frame for next analysis of
#' other functions of the package.
#'
#' @param effort_base \code{data.frame} with effort data.
#' @param covariable \code{vector} of covariable names to keep in output of the function.
#' @param block_area Data.frame with 2 colnames :
#'                     \enumerate{
#'                       \item Block.
#'                       \item Area.
#'                     }
#' @param shape Shapefile of the study area. It can be either a \code{SpatialPolygonsDataFrame} class object,
#'  in this case it is not necessary to give shape_layer argument. Or it can be the name of the shape
#'  object with its extension ".shp" (ex : "data/studyAreaShapefile.shp").
#' @param shape_layer Layer of the shapefile if shape is a \code{character string}.
#' @param New_projection New projection of longitude and latitude of columns
#'        (POINT_X et POINT_Y) in \code{Proj4String} format,
#'        see : \code{\link[sp]{CRS}} for more infos.
#' @param optimal Argument which allows to keep data sampled in optimal conditions.
#'        Defaults settings are all data are kept. In case of optimal = T, indexes
#'        \code{"c("GG", "GM", "MG", "EG", "GE", "EE", "ME", "EM", "MM")"} are kept.
#' @param col2keep Columns to keep from effort_base in output of the function.
#' @return This function return a list containing:
#'         \enumerate{
#'           \item legdata : \code{data.frame} with infos at leg scale.
#'           \item segdata : \code{data.frame} with infos at segment scale.
#'         }
#' @import dplyr sp rgdal
#' @importMethodsFrom raster as.character
#' @importFrom lubridate ymd
#' @examples
#'
#' @export


prepare_data_effort <- function(effort_base,
                                covariable = NULL,
                                block_area,
                                shape,
                                shape_layer,
                                New_projection,
                                optimal = FALSE,
                                col2keep = NULL) {

  effort <- effort_base

  # verifier si les colonnes sont bien dans le DF effort
  col_name_neces <- c("lon","lat","seaState","subjective",
                      "survey","strateSec","transect","legId","segLength",
                      "segId","left","right","CenterTime")


  if(!all(col_name_neces %in% colnames(effort))){
    var_alone <- setdiff(col_name_neces, colnames(effort))

    stop(paste("Les variables : ",var_alone, "ne sont pas dans le tableau effort.", sep="\n",collapse = ", "),
         paste("utiliser la fonction change_effort_varName","\n",sep=""))
  }

  ### Block et Surface
  # bon nom pour block et area
  if(!is.data.frame(block_area)) {
    stop("block_area n'est pas un objet de type data.frame.")
  } else if(!all(names(block_area) %in% c("Block","Area"))) {
    stop(paste('block_area doit contenir les colonnes : "Block" et "Area".',sep=''))
  } else {
    block_area <- block_area
  }
  # correspondance entre strateSec et block area
  if(!any(block_area$Block %in% effort$strateSec)){
    stop(cat("la variable Block de block_area ne correspond pas aux valeurs de strateSec dans la table effort :",
               "- soit le strateSec créé par change_varName n'est pas du bon format",
               "- soit les valeurs de Block ne sont pas dans dans la variable strateSec de la table effort",sep="\n"))
  }

  # polygons sampling
  if(any("character" %in% is(shape))){
    poly_NC <- readOGR(dsn = paste(shape), layer = shape_layer, verbose = F) # NC
  } else {
    poly_NC <- shape
  }


  # Covariable
  allvar = covariable

  # Grille de la zone en 2008
  effort$POINT_X <- effort$lon
  effort$POINT_Y <- effort$lat

  # Reprojeter dans système de projection  renseingé en argument
  effort_xy <- effort[, c("POINT_X", "POINT_Y")]
  coordinates(effort_xy) <- ~ POINT_X + POINT_Y
  effort_xy@proj4string <- CRS(as.character(poly_NC@proj4string))
  effort_xy <- spTransform(effort_xy, CRS(New_projection))
  effort[, c("POINT_X", "POINT_Y")] <- effort_xy@coords

  # selection effort et obs en bonnes conditions
  if(optimal==T) {
    effort <- subset(effort, seaState <= 3 & subjective %in% c("GG", "GM", "MG", "EG", "GE", "EE", "ME", "EM", "MM"))
  }


  # creation des tableaux necessaires a l'analyse  ---------------

  #-- legdata --#
  #-------------#
  if("session" %in% colnames(effort)) {
    legdata <- effort %>%
      group_by(survey, strateSec, transect, legId, left, right, session) %>%
      summarize(Effort = sum(segLength),
                seaState = unique(seaState),
                subjective = unique(subjective))
  } else {
    legdata <- effort %>%
      group_by(survey, strateSec, transect, legId, left, right) %>%
      summarize(Effort = sum(segLength),
                seaState = unique(seaState),
                subjective = unique(subjective))
  }

  legdata <- as.data.frame(legdata)

  names(legdata)[which(names(legdata) %in% c("strateSec", "transect", "legId"))] <- c("Region.Label", "Transect.Label",
                                                                                       "Sample.Label")

  # merge col2keep
  if(!is.null(col2keep)){

    col2keep <- col2keep[!(col2keep %in% colnames(legdata))]

    legdata_col_joined <- legdata %>%
      left_join(effort_base %>%
                  select(.dots = c(col2keep,"legId")) %>%
                  `colnames<-`(c(col2keep,"legId")),
                by = c("Sample.Label"="legId"))

    if(nrow(legdata) < nrow(legdata_col_joined)){
      stop("col2keep are not unique at leg scale")
    }

    legdata <- legdata_col_joined

  }

  # Assigner area à legdata en fonction du nom du block commun avec block_area
  legdata$Area <- sapply(legdata$Region.Label, function(id) {
    as.numeric(block_area$Area[which(block_area$Block == id)])
  })



  #-- segdata --#
  #-------------#
  if("session" %in% colnames(effort)) {
    segdata <- data.frame(effort[, c("CenterTime", "survey", "transect", "legId", "segId",
                                     "segLength", "POINT_X",
                                     "POINT_Y", "lon", "lat", "strateSec",
                                     "seaState", "subjective","session", allvar, col2keep)
                                 ])
  } else {
    segdata <- data.frame(effort[, c("CenterTime", "survey", "transect", "legId", "segId",
                                     "segLength", "POINT_X",
                                     "POINT_Y", "lon", "lat", "strateSec",
                                     "seaState", "subjective", allvar, col2keep)
                                 ])
  }


  segdata$CenterTime <- ymd(segdata$CenterTime)
  names(segdata)[1:11] <- c("date", "survey", "Transect.Label", "Sample.Label", "Seg", "Effort", "X", "Y",
                            "longitude", "latitude", "Region.Label")

  # Assigner area à segdata en fonction du nom du block commun avec block_area
  segdata$Area <- sapply(segdata$Region.Label, function(id) {block_area$Area[which(block_area$Block == id)]})


  # renvoyer les outputs
  return(list(legdata = legdata, segdata = segdata))

}
