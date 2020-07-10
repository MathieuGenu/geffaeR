#' Adding observation on effort.
#'
#' Adding the number of individuals and the group size on prepared effort data
#' (e.g. on legdata and segdata).
#'
#' @inheritParams prepare_data_obs
#' @param countdata_leg data.frame containing the detection number and the total number of
#' individuals for legs for which observation have beebn made.
#' @param coundata_seg data.frame containing the detection number and the total number of
#' individuals for segments for which observation have beebn made.
#' @return This function return a list containing :
#'         \enumerate{
#'           \item legdata_obs : data.frame corresponding to legdata with observation information.
#'           \item segdata_obs : data.frame corresponding to segdata with observation information.
#'         }
#' @examples
#'
#'
#' @export


ajout_obs <- function(countdata_leg, legdata, countdata_seg, segdata){

  ### merge Legdata avec obs pour n_detected et n_ind ###
  legdata_obs <- merge(legdata, countdata_leg, by = c("Transect.Label","Sample.Label"),
                       all.x = TRUE)
  legdata_obs$n_detected[is.na(legdata_obs$n_detected)] <- 0
  legdata_obs$n_ind[is.na(legdata_obs$n_ind)] <- 0

  ### merge Segdata avec obs pour n_detected et n_ind ###
  segdata_obs <- merge(segdata, countdata_seg, by = c("Transect.Label","Seg","Sample.Label"),
                       all.x = TRUE)
  segdata_obs$n[is.na(segdata_obs$n)] <- 0
  segdata_obs$y[is.na(segdata_obs$y)] <- 0

  return(list(segdata_obs = segdata_obs, legdata_obs = legdata_obs))
}


