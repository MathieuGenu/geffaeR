#' @export

cv_dsm <- function(dsm.fc, predata, inside = NULL) {
  if(is.null(inside)) { preddata.varprop <- predata[which(predata[, "inside"] == 1), ] }
  else { preddata.varprop <- predata[which(predata$zone %in% inside), ] }

  preddata.list <- split(preddata.varprop, 1:nrow(preddata.varprop))
  dsm.xy.varprop <- dsm.var.gam(dsm.fc,
                                pred.data = preddata.list,
                                off.set = preddata.varprop$Area
  )
  # colone pour la variance dans predata
  preddata.varprop$var <- as.numeric(dsm.xy.varprop$pred.var)
  # colone pour le CV dans predata
  preddata.varprop$cv <- with(preddata.varprop, round (100*sqrt(var)/pred, 2))
  # results
  return(list(var.dsm = dsm.xy.varprop, preddata.varprop = preddata.varprop))
}
