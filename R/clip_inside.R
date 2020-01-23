#' @export


clip_inside <- function(predata, poly_df) {
  spdf_predata <- predata
  coordinates(spdf_predata) <- ~ longitude + latitude
  spdf_predata@proj4string <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  if(as.character(spdf_predata@proj4string) != as.character(poly_df@proj4string)) {
    spdf_predata <- spTransform(spdf_predata, poly_df@proj4string)
  }
  if(length(poly_df) > 1) {
    dd <- do.call("cbind",
                  lapply(1:nrow(poly_df@data),
                         function(i) {
                           x_inside <- rep(0, nrow(predata))
                           x_inside[which(!is.na(over(spdf_predata, poly_df[i, ])$Id))] <- 1
                           return(x_inside)
                         }
                  )
    )
    dd <- apply(dd, 1, sum)
    return(dd)
  }
  else {
    x_inside <- rep(0, nrow(predata))
    x_inside[which(!is.na(over(spdf_predata, poly_df)$Id))] <- 1
    return(x_inside)
  }
}
