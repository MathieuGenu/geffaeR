#' \encoding{Standardisation des noms de colonnes}
#'
#' \encoding{Cette fonction permet d'attribuer au data.frame effort les bons noms de colonnes pour
#' la suite des analyse à effectuer.}
#'
#' @param effort_base \encoding{data.frame contenant les données d'effort pour lesquelles les variables sont à standardiser.}
#'
#' @return \encoding{Cette fonction renvoie un data.frame effort avec les noms de variables au bon format
#'         pour la suite de l'analyse.}
#'
#' @note \encoding{Lorsque strate_sec n'est pas trouvé à la fin de la fonction, celle-ci créée un strate_sec
#'       correspondant à la moitié de transect plus survey séparés par un "_".
#'       Example, pour une ligne donnée, si \code{transect = O2/306} et \code{survey = FR_N_C} on aura alors
#'       \code{strate_sec = FR_N_c_02}.}
#'
#' @examples
#'
#'
#' @export



change_effort_varName <- function(effort_base){

  col_name_neces <- c("lon","lat","seaState","subjective",
                      "survey","strate_sec","transect","IdLeg","Shape_Leng","segId","left_","right_")

  if("session_" %in% colnames(effort_base)){
    effort_base$session_ <- as.factor(effort_base$session_)
  }

  if(all(col_name_neces %in% colnames(effort_base))){
    return(effort_base)
  } else {
    var_alone <- col_name_neces[!(col_name_neces %in% colnames(effort_base)) ]

    # which var is not right
    pos_var <- which(col_name_neces %in% var_alone)
    # lon
    if("strate_sec" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("STRATE_SEC","subRegion_strate")] <- "strate_sec"
    }
    if("lon" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("point_X","Point_X","pointX","PointX",
                                                         "POINT_X","POINTX")] <- "lon"
    }
    if("lat" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("point_Y","Point_Y","pointY","PointY",
                                                         "POINT_Y","POINTY")] <- "lat"
    }
    if("seaState" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("Sea_State","SeaState","sea_state","SEA_STATE",
                                                         "Beaufort","BEAUFORT","beaufort")] <- "seaState"
    }
    if("subjective" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("Subjective","SUBJECTIVE")] <- "subjective"
    }
    if("survey" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("Survey","Campaign","CAMPAIGN")] <- "survey"
    }
    if("transect" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("TRANSECT","Transect")] <- "transect"
    }
    if("IdLeg" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("IDLEG","Id_Leg","Id_leg","Sample.Label")] <- "IdLeg"
    }
    if("Shape_Leng" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("Length","length","LengthKm",'LengthKM')] <- "Shape_Leng"
    }
    if("segId" %in% var_alone){
      # cas particulier peut y avoir plusieurs segId dans la table effort (5k et 10k) -> prendre celui sans 0
      effort_match <- effort_base[,colnames(effort_base) %in% c("IdSeg","Id_Seg","IDSEG","Seg","seg","SegID",
                                                      "segId5k","segId10k","segId5K","segId10K"), drop = FALSE]
      colnames(effort_base)[colnames(effort_base) %in% names(effort_match)[colSums(effort_match) > 1]] <- "segId"
    }
    if("left_" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("LEFT_REAR")] <- "left_"
    }
    if("right_" %in% var_alone){
      colnames(effort_base)[colnames(effort_base) %in% c("RIGHT_REAR")] <- "right_"
    }
    if(!all(is.na(as.numeric(as.character(effort_base$IdLeg))))) {
      if("session_" %in% colnames(effort_base)) {
        effort_base$IdLeg <- paste(effort_base$IdLeg,effort_base$session_,sep="_")
      }
    }
    if(!all(is.na(as.numeric(as.character(effort_base$segId))))) {
      if("session_" %in% colnames(effort_base)) {
        effort_base$segId <- paste(effort_base$segId,effort_base$session_,sep="_")
        colnames(effort_base)[colnames(effort_base) %in% c("session_")] <- "session"
      }
    }

    # tester si strate_sec est fabriqué, sinon le fabriquer avec "survey" et "transect" si ils existent...
    if(!("strate_sec" %in% colnames(effort_base)) & all(c("survey","transect") %in% colnames(effort_base))){
      split <- strsplit(effort_base$transect,"/")
      half_transect <- sapply( split, "[", 1)
      effort_base$strate_sec <- paste(effort_base$survey,half_transect,sep="_")
    }
    # tester si c'est OK maintenant sino renvoyer un message et dire quel équivalent est absent
    if(!all(col_name_neces %in% colnames(effort_base))){
      var_still_alone <- col_name_neces[!(col_name_neces %in% colnames(effort_base)) ]
      message(paste(c("Pour la base effort, la fonction ne trouve pas d'équivalent pour : ", var_still_alone), collapse="\n "))
    } else {
      return(effort_base = effort_base)
    }
  }
}

#' \encoding{Standardisation des noms de colonnes de la base observation}
#'
#' \encoding{Cette fonction permet d'attribuer au data.frame observation les bons noms de colonnes pour
#' la suite des analyses à effectuer.}
#'
#' @param obs_base \encoding{data.frame contenant les données d'observation pour lesquelles les variables
#'        sont à standardiser.}
#'
#' @return \encoding{Cette fonction renvoie un data.frame effort avec les noms de variables au bon format
#'         pour la suite de l'analyse.}
#'
#' @examples
#'
#'
#' @export



change_obs_varName <- function(obs_base) {


  col_name_neces <- c("strate","subRegion","transect","taxon","group_",
                      "family","species","podSize","decAngle","lat","lon",
                      "IdLeg","PerpDist","segId")
  if("session_" %in% colnames(obs_base)) {
    obs_base$session_ <- as.factor(obs_base$session_)
  }
  if("observer" %in% colnames(obs_base)) {
    obs_base$observerId <- as.factor(obs_base$observer)
    obs_base$observer <- NULL
  }

  if(all(col_name_neces %in% colnames(obs_base))) {
    return(obs_base)
  } else {
    var_alone <- col_name_neces[!(col_name_neces %in% colnames(obs_base)) ]

    # which var is not right
    pos_var <- which(col_name_neces %in% var_alone)

    # regarder ce qu'il reste si on enleve les majuscules
    varName_to_lower <- tolower(colnames(obs_base))
    varTrue <- varName_to_lower %in% var_alone
    colnames(obs_base)[varTrue] <- varName_to_lower[varTrue]

    var_alone_after_lower <- col_name_neces[!(col_name_neces %in% colnames(obs_base)) ]

    # lon
    if("lon" %in% var_alone_after_lower) {
      colnames(obs_base)[colnames(obs_base) %in% c("point_X","Point_X","pointX","PointX","POINT_X","POINTX")] <- "lon"
    }
    if("lat" %in% var_alone_after_lower) {
      colnames(obs_base)[colnames(obs_base) %in% c("point_Y","Point_Y","pointY","PointY","POINT_Y","POINTY")] <- "lat"
    }
    if("subRegion" %in% var_alone_after_lower) {
      colnames(obs_base)[colnames(obs_base) %in% c("SECTEUR","Secteur")] <- "subRegion"
    }
    if("PerpDist" %in% var_alone_after_lower) {
      colnames(obs_base)[colnames(obs_base) %in% c("Distance","distance")] <- "PerpDist"
    }
    if("IdLeg" %in% var_alone_after_lower) {
      colnames(obs_base)[colnames(obs_base) %in% c("IDLEG","Id_Leg","Id_leg","Sample.Label")] <- "IdLeg"
    }
    if("segId" %in% var_alone_after_lower) {
      obs_match <- obs_base[,colnames(obs_base) %in% c("IdSeg","Id_Seg","IDSEG","Seg","seg","SegID",
                                                       "segId5k","segId10k","segId5K","segId10K"),
                            drop = F]
      colnames(obs_base)[colnames(obs_base) %in% names(obs_match)[colSums(obs_match) > 1]] <- "segId"
    }
    if("decAngle" %in% var_alone_after_lower) {
      colnames(obs_base)[colnames(obs_base) %in% c("DEC_ANGLE","Dec_Angle","dec_angle","dec_Angle")] <- "decAngle"
    }
    if("podSize" %in% var_alone_after_lower) {
      colnames(obs_base)[colnames(obs_base) %in% c("POD_SIZE","Pod_Size","pod_size")] <- "podSize"
    }
    if("taxon" %in% var_alone_after_lower) {
      colnames(obs_base)[colnames(obs_base) %in% c("TAXON_Fr")] <- "taxon"
    }
    if(!all(is.na(as.numeric(as.character(obs_base$IdLeg))))) {
      if("session_" %in% colnames(obs_base)) {
        obs_base$IdLeg <- paste(obs_base$IdLeg,obs_base$session_,sep="_")
      }
    }
    if(!all(is.na(as.numeric(as.character(obs_base$segId))))) {
      if("session_" %in% colnames(obs_base)) {
        obs_base$segId <- paste(obs_base$segId,obs_base$session_,sep="_")
        colnames(obs_base)[colnames(obs_base) %in% c("session_")] <- "session"
      }
    }
    # tester si c'est OK maintenant sino renvoyer un message et dire quel équivalent est absent
    if(!all(col_name_neces %in% colnames(obs_base))) {
      var_still_alone <- col_name_neces[!(col_name_neces %in% colnames(obs_base)) ]
      message(paste(c("Pour la base observation, la fonction ne trouve pas d'équivalent pour : ", var_still_alone), collapse="\n "))
    } else {
      return(obs_base = obs_base)
    }
  }
}
