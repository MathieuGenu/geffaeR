#' Prediction map of fitted dsm model.
#'
#' @param predata data.frame of prediction.
#' @param grid gridata in \code{sf} format or \code{crs} value. It have to be usable by \code{\link[sf]{st_crs}}.
#' @param var Column name as character used as value to be plotted.
#' @param facet_param Column name as character used as facetting parameter.
#' @param title character string to put as title on the map.
#'
#' @return \code{ggplot} object.
#' @export
#'
#' @import dplyr sf sp ggplot2
#' @importFrom viridis scale_fill_viridis
#' @examples
pred_map_dsm <- function(predata, grid, var, facet_param = NULL, title = NULL) {

  sf_predict_reprojected <- predata %>%
    st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
    st_transform(st_crs(grid))

  predict_reprojected <- sf_predict_reprojected %>%
    as_Spatial() %>%
    as.data.frame() %>%
    rename(lon = coords.x1, lat = coords.x2 )

  no_sf <- sf_predict_reprojected %>%
    as_Spatial() %>%
    as.data.frame() %>%
    rename(lon = coords.x1, lat = coords.x2 )

  width <- no_sf %>%
    group_by(lon) %>%
    summarize(n=n()) %>%
    arrange(lon) %>%
    mutate(diff = c(NA,diff(lon, lag=1))) %>%
    summarize(max = max(diff, na.rm = T)) %>%
    as.numeric()

  height <- no_sf %>%
    group_by(lat) %>%
    summarize(n=n()) %>%
    arrange(lat) %>%
    mutate(diff = c(NA,diff(lat, lag=1))) %>%
    summarize(max = max(diff, na.rm = T)) %>%
    as.numeric()


  bbox_study_area <- st_bbox(sf_predict_reprojected)

  g_map_pred <- ggplot() +
    geom_sf(data = sf_predict_reprojected, colour = alpha(colour = "black", alpha = 0)) +
    geom_tile(data = no_sf, mapping = aes(x = lon, y = lat, fill = get(var), width = width, height = height)) +
    geom_sf(data = st_as_sf(NEA_simplified_FR), colour = "grey50") +
    coord_sf(xlim = bbox_study_area[c(1,3)],
             ylim = bbox_study_area[c(2,4)]) +
    viridis::scale_fill_viridis() +
    labs(fill = "Density (ind/km²)") +
    theme_bw()

  if(!is.null(facet_param)) {
    g_map_pred <- g_map_pred + facet_wrap(get(facet_param)~.)
  }

  if(!is.null(title)) {
    g_map_pred <- g_map_pred + ggtitle(title)
  }


  return(g_map_pred)
}
