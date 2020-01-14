#' \encoding{Ajout des Observations sur les données effort.}
#'
#' \encoding{Cette fonction permet d'ajouter le nombre d'individus/groupe et le nombre d'observation totale
#' aux données d'effort preparées (i.e. a legdata et segdata).}
#'
#' @inheritParams prepare_data_obs
#' @param countdata_leg \encoding{Un data.frame contenant le nombre de détection et d'individus total pour les legs
#'           où il y a eu observation pour l'espèce ou group en question.}
#' @param coundata_seg \encoding{Un data.frame contenant le nombre de détection et d'individus total pour les segments
#'           où il y a eu observation pour l'espèce ou group en question.}
#'
#' @return Cette fonction renvoit une liste contenant :
#'         \enumerate{
#'           \item legdata_obs : Un data.frame contenant Legdata avec des infos d'observations.
#'           \item segdata_obs : Un data.frame contenant Segdata avec des infos d'observations.
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


