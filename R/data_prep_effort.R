#' \encoding{Préparation des données d'effort.}
#'
#' \encoding{Cette fonction prépare les données d'effort pour qu'elles soient compatibles
#' pour les autres fonctions de ce package.}
#'
#' @param effort_base \encoding{data.frame contenant les données d'effort.}
#' @param covariable \encoding{Vecteur contenant le nom des covariables à conserver en sortie de la fonction.}
#' @param block_area \encoding{data.frame contenant une colonne :
#'                     \enumerate{
#'                       \item Block.
#'                       \item Area.
#'                     }
#'                   }
#' @param shape \encoding{Nom du shape file de la zone d'étude.}
#' @param shape_layer \encoding{Couche du shape file d'interêt.}
#' @param New_projection \encoding{Nouvelle projection des longitude et latitude des
#'        colonnes (POINT_X et POINT_Y) format Proj4String,
#'        voir : \code{\link[sp]{CRS}} pour plus d'infos.}
#' @param optimal \encoding{Paramètre permettant de selectionner (ou non) les données prelevées dans les
#'        conditions otpimales. Par défaut toutes les données sont selectionnées. Dans le cas où seules les
#'        données en conditions optimales sont selectionnees, ce sont les indices
#'        \code{c("GG", "GM", "MG", "EG", "GE", "EE", "ME", "EM", "MM")} qui sont conservés.}
#' @return Cette fonction renvoie une liste  contenant :
#'         \enumerate{
#'           \item Legdata : Un data.frame contenant les infos sur chaque leg.
#'           \item Segdata : Un data.frame contenant les infos sur chaque segment.
#'         }
#' @examples
#'
#'
#' @export


prepare_data_effort <- function(effort_base, covariable = NULL, block_area, shape, shape_layer,
                                New_projection, optimal = FALSE) {

  effort <- effort_base

  # verifier si les colonnes sont bien dans le DF effort
  col_name_neces <- c("lon","lat","seaState","subjective",
                      "survey","strate_sec","transect","IdLeg","Shape_Leng",
                      "segId","left_","right_")
  if(!all(col_name_neces %in% colnames(effort))){
    var_alone <- col_name_neces[!(col_name_neces %in% colnames(effort))]
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
  # correspondance entre strate_sec et block area
  if(!any(block_area$Block %in% effort$strate_sec)){
    stop(cat("la variable Block de block_area ne correspond pas aux valeurs de strate_sec dans la table effort :",
               "- soit le strate_sec créé par change_varName n'est pas du bon format",
               "- soit les valeurs de Block ne sont pas dans dans la variable strate_sec de la table effort",sep="\n"))
  }

  # polygons sampling
  poly_NC <- readOGR(dsn = paste(DataDir, shape, sep = "/"), layer = shape_layer, verbose = F) # NC


  # Covariable
  allvar = covariable

  # Grille de la zone en 2008
  effort$POINT_X <- effort$lon
  effort$POINT_Y <- effort$lat

  # Reprojeter dans système de projection  renseingé en argument
  effort_xy <- effort[, c("POINT_X", "POINT_Y")]
  coordinates(effort_xy) <- ~ POINT_X + POINT_Y
  effort_xy@proj4string <- CRS(New_projection)
  effort_xy <- spTransform(effort_xy, CRS(as.character(poly_NC@proj4string)))
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
        group_by(survey, strate_sec, transect, IdLeg, left_, right_, session) %>%
        summarize(Effort = sum(Shape_Leng),
                  seaState = unique(seaState),
                  subjective = unique(subjective))
    } else {
      legdata <- effort %>%
        group_by(survey, strate_sec, transect, IdLeg, left_, right_) %>%
        summarize(Effort = sum(Shape_Leng),
                  seaState = unique(seaState),
                  subjective = unique(subjective))
    }

  legdata <- as.data.frame(legdata)
  names(legdata)[which(names(legdata) %in% c("strate_sec", "transect", "IdLeg"))] <- c("Region.Label", "Transect.Label",
                                                                                       "Sample.Label")

  # Assigner area à legdata en fonction du nom du block commun avec block_area
  legdata$Area <- sapply(legdata$Region.Label, function(id) {as.numeric(block_area$Area[which(block_area$Block == id)])})



  #-- segdata --#
  #-------------#
  if("session" %in% colnames(effort)) {
    segdata <- data.frame(effort[, c("CenterTime", "survey", "transect", "IdLeg", "segId", "Shape_Leng", "POINT_X",
                                     "POINT_Y", "lon", "lat", "strate_sec",
                                     "seaState", "subjective","session", allvar)
                                 ])
  } else {
    segdata <- data.frame(effort[, c("CenterTime", "survey", "transect", "IdLeg", "segId", "Shape_Leng", "POINT_X",
                                     "POINT_Y", "lon", "lat", "strate_sec",
                                     "seaState", "subjective", allvar)
                                 ])
  }


  segdata$CenterTime <- lubridate::ymd(segdata$CenterTime)
  names(segdata)[1:11] <- c("date", "survey", "Transect.Label", "Sample.Label", "Seg", "Effort", "X", "Y",
                            "longitude", "latitude", "Region.Label")

  # Assigner area à segdata en fonction du nom du block commun avec block_area
  segdata$Area <- sapply(segdata$Region.Label, function(id) {block_area$Area[which(block_area$Block == id)]})


  # renvoyer les outputs
  return(list(legdata = legdata, segdata = segdata))

}
