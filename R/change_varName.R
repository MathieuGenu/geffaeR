#' Standardisation of effort data column names.
#'
#' Attribute standard names to column names of effort data.frame
#' for the next steps of the analysis.
#'
#' @param effort_base Data.frame containing column names to standardize.
#'
#' @return This function return the same data.frame as in input with the standardized column name.
#'
#' @note When "strate_sec" is not found at the end of the function, it creates a strate_sec
#'       column automatically corresponding to the merge of "subRegion" and "strate". There are merged by "_".
#'       Example, for a row, if \code{strate = N1} and \code{subRegion = ATL} we will get
#'       \code{strateSec = ATL_N1}.
#' @export
#' @import assertthat

change_effort_varName <- function(effort_base){

  assert_that(is.data.frame(effort_base))

  col_name_neces <- c("lon","lat","seaState","subjective","subRegion","strate", "survey",
                      "strateSec","transect","legId","segLength","segId","left","right", "CenterTime")

  if(all(col_name_neces %in% colnames(effort_base))){

    return(effort_base)

  } else {

    # Put everything tolowercase and without "_" and check if there are matches with colnames ----
    lower_colnames <- tolower(colnames(effort_base))
    lower_no_under_colnames <- gsub("_", "", lower_colnames)
    lower_neces <- tolower(col_name_neces)

    for (n in 1:length(lower_neces)) {
      if (lower_neces[n] %in% lower_no_under_colnames) {
        colnames(effort_base)[lower_no_under_colnames %in% lower_neces[n]] <- col_name_neces[n]
      }
    }

    # Are there all the needed columns ? ----
    missing_needed_col <- which(!(col_name_neces %in% colnames(effort_base)))

    if (!("lon" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("pointx","longitude")] <- "lon"
    }
    if (!("lat" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("pointy","latitude")] <- "lat"
    }
    if (!("seaState" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("seastate","beaufort")] <- "seaState"
    }
    if (!("subjective" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("Subjective")] <- "subjective"
    }
    if (!("subRegion" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("subregion")] <- "subRegion"
    }
    if (!("strate" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("strate")] <- "strate"
    }
    if (!("survey" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("survey","campaign")] <- "survey"
    }
    if (!("strateSec" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("stratesec","subregstra")] <- "strateSec"
    }
    if (!("transect" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("transect")] <- "transect"
    }
    if (!("legId" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("legid","idleg","samplelabel")] <- "legId"
    }
    if (!("segLength" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("seglength","length","lengthkm",
                                                           "segleng10k","segleng10km","segleng5km",
                                                           "segleng5k", "seglength5k","seglength5km",
                                                           "seglength5","seglength10k","seglength10km",
                                                           "seglength10")] <- "segLength"
    }
    if (!("segId" %in% colnames(effort_base))) {
      # particularity for segId
      # there are potentially 2 columns segId (10km and 5km)
      # check if there is only one segId
      #  if yes -> replace
      #  if not -> check of one of them has only 0 in it
      #    if yes -> keep one with no 0
      #    if not -> error message, delete one
      potential_segId <- which(lower_no_under_colnames %in% c("idseg","seg","segid",
                                                              "segid5km","segid10km",
                                                              "segid5k","segid10k"))
      if (length(potential_segId) == 1) {
        colnames(effort_base)[potential_segId] <- "segId"
      } else {

        # TO DELETE, user should select its segId by itself #
        length_unique_col <- sapply(effort_base[potential_segId], function(x){length(unique(x))})
        if (length(length_unique_col[length_unique_col > 3]) == 1) {
          colnames(effort_base)[colnames(effort_base) %in% names(length_unique_col[length_unique_col])] <- "segId"
        } else {
        #####################################################
          stop("Can't determine a segId column because there are several segId columns(5k, 10k or more...)\n
             choose between one of them.")
        }
      }
    }

    # Replace facultatif colnames ----
    if (!("left" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("left")] <- "left"
    }
    if (!("right" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("right")] <- "right"
    }
    if (!("session" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("session")] <- "session"
    }

    # Add session in seg_id if both exist ----
    if (all(c("segId","session") %in% colnames(effort_base))) {
      effort_base$segId <- paste(effort_base$segId,effort_base$session,sep="_")
    }

    # Add session in legId if both exist ----
    if (all(c("legId","session") %in% colnames(effort_base))) {
      effort_base$legId <- paste(effort_base$legId,effort_base$session,sep="_")
    }

    # Build strate_sec if it doesn't exist ----
    if (!("strateSec" %in% colnames(effort_base))) {
      effort_base$strateSec <- paste(effort_base$subRegion,effort_base$strate,sep="_")
    }

    # Error message if needed columns are missing ----
    if (!all(col_name_neces %in% colnames(effort_base))) {
      var_still_alone <- col_name_neces[!(col_name_neces %in% colnames(effort_base))]
      stop(paste(c("For effort data, can't find equivalent name for : ", var_still_alone), collapse="\n"))
    }


    return(effort_base)
  }

}

#' Standardisation of observation data column names.
#'
#' Attribute standard names to column names of observation data.frame
#' for the next steps of the analysis.
#'
#' @param obs_base Data.frame containing column names to standardize.
#'
#' @return This function return the same data.frame as in input with the standardized column name.
#' @export
#' @importFrom assertthat assert_that



change_obs_varName <- function(obs_base) {

  assert_that(is.data.frame(obs_base))

  col_name_neces <- c("strate","subRegion","transect","taxon","group",
                      "family","species","podSize","decAngle","lat","lon",
                      "legId","perpDist","segId", "survey")

  if(all(col_name_neces %in% colnames(obs_base))){

    return(obs_base)

  } else {

    # Put everything tolowercase and without "_" and check if there are matches with colnames ----
    lower_colnames <- tolower(colnames(obs_base))
    lower_no_under_colnames <- gsub("_", "", lower_colnames)
    lower_neces <- tolower(col_name_neces)

    for (n in 1:length(lower_neces)) {
      if (lower_neces[n] %in% lower_no_under_colnames) {
        colnames(obs_base)[lower_no_under_colnames %in% lower_neces[n]] <- col_name_neces[n]
      }
    }

    # Are there all the needed columns ? ----
    missing_needed_col <- which(!(col_name_neces %in% colnames(obs_base)))

    if (!("strate" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("strate")] <- "strate"
    }
    if (!("subRegion" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("subregion","secteur")] <- "subRegion"
    }
    if (!("transect" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("transect")] <- "transect"
    }
    if (!("taxon" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("taxon","taxonfr")] <- "taxon"
    }
    if (!("group" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("group","groupe","groupefr")] <- "group"
    }
    if (!("family" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("family","famille","famillefr")] <- "family"
    }
    if (!("species" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("species")] <- "species"
    }
    if (!("podSize" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("podsize","size")] <- "podSize"
    }
    if (!("decAngle" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("decangle")] <- "decAngle"
    }
    if (!("lon" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("pointx","longitude")] <- "lon"
    }
    if (!("lat" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("pointy","latitude")] <- "lat"
    }
    if (!("strateSec" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("stratesec","subregstra")] <- "strateSec"
    }
    if (!("legId" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("legid","idleg","samplelabel")] <- "legId"
    }
    if (!("perpDist" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("perpdist","distance")] <- "perpDist"
    }
    if (!("segId" %in% colnames(obs_base))) {
      # particularity for segId
      # there are potentially 2 columns segId (10km and 5km)
      # check if there is only one segId
      #  if yes -> replace
      #  if not -> check of one of them has only 0 in it
      #    if yes -> keep one with no 0
      #    if not -> error message, delete one
      potential_segId <- which(lower_no_under_colnames %in% c("idseg","seg","segid",
                                                              "segid5km","segid10km",
                                                              "segid5k","segid10k"))
      if (length(potential_segId) == 1) {
        colnames(obs_base)[potential_segId] <- "segId"
      } else {
        length_unique_col <- sapply(obs_base[potential_segId], function(x){length(unique(x))})
        if (length(length_unique_col[length_unique_col > 3]) == 1) {
          colnames(obs_base)[colnames(obs_base) %in% names(length_unique_col[length_unique_col])] <- "segId"
        } else {
          stop("Can't determine a segId column because there are several segId columns(5k, 10k or more...)\n
             choose between one of them.")
        }
      }
    }
    if (!("survey" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("survey")] <- "survey"
    }

    # Replace facultatif colnames ----
    if (!("left" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("left")] <- "left"
    }
    if (!("right" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("right")] <- "right"
    }
    if (!("session" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("session")] <- "session"
    }

    # Transform session in factor ----
    if ("session" %in% colnames(obs_base)) {
      obs_base$session <- as.factor(obs_base$session)
    }

    # Add session in seg_id if both exist ----
    if (all(c("segId","session") %in% colnames(obs_base))) {
      obs_base$segId <- paste(obs_base$segId,obs_base$session,sep="_")
    }

    # Add session in legId if both exist ----
    if (all(c("legId","session") %in% colnames(obs_base))) {
      obs_base$legId <- paste(obs_base$legId,obs_base$session,sep="_")
    }

    # Build strate_sec if it doesn't exist ----
    if (!("strateSec" %in% colnames(obs_base))) {
      obs_base$strateSec <- paste(obs_base$subRegion,obs_base$strate,sep="_")
    }
    # Transform observer into observerId ----
    if("observer" %in% colnames(obs_base)) {
      obs_base$observerId <- as.factor(obs_base$observer)
      obs_base$observer <- NULL
    }
    # Error message if needed columns are missing ----
    if (!all(col_name_neces %in% colnames(obs_base))) {
      var_still_alone <- col_name_neces[!(col_name_neces %in% colnames(obs_base))]
      stop(paste(c("For observation data, can't find equivalent name for : ", var_still_alone), collapse="\n "))
    }
    return(obs_base)
  }
}
