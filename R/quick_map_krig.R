#' Visualize kriging prediction.
#'
#' Map projection of \code{\link[geffaeR]{predict_kriging}} output.
#'
#'
#' @param df_sf Object of class "sf", output of \code{\link[geffaeR]{predict_kriging}}.
#' @param varname Choose which column of df to plot between "mean_pred" or "se_pred".
#' @param segdata_obs Optionnal. Allows to plot observations on the map,
#' default is NULL. Output of \code{\link[geffaeR]{ajout_obs}}.
#' @import ggplot2
#' @importFrom ggspatial annotation_north_arrow annotation_scale
#' @importFrom sf st_as_sf
#' @importFrom viridisLite viridis
#' @export

quick_map_krig <- function(df_sf, varname, segdata_obs = NULL) {

  df_dt <- df_sf %>%
    as_Spatial() %>%
    as.data.frame() %>%
    rename(lon = coords.x1, lat = coords.x2)

  # kriging graph
  NEA_sf <- st_as_sf(NEA_simplified_FR) %>%
    st_transform(st_crs(df_sf))

  # get graph window values
  range_x <- range(df_dt[,"lon"])
  range_y <- range(df_dt[,"lat"])

  x_lim  <- range_x + c((0.05 * diff(range_x))*c(-1,1))
  y_lim  <- range_y + c((0.05 * diff(range_y))*c(-1,1))


  g_graph <- ggplot() +
    geom_tile(data = df_dt,
              aes(x = lon, y = lat, fill = get(varname))) +
    geom_sf(data = NEA_sf) +
    coord_sf(xlim = x_lim, ylim = y_lim, expand = FALSE) +
    scale_fill_gradientn(colours = viridisLite::viridis(256),
                         name = varname) +
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "tr",
                           which_north = "true",
                           pad_x = unit(0.1, "cm"),
                           pad_y = unit(0.05, "cm"),
                           style = north_arrow_fancy_orienteering) +
    theme(legend.position = "bottom",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 8, angle = 30),
          panel.grid = element_line(colour = "transparent"),
          plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
          legend.box = "vertical",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(size = "ln(nb_ind+1)")

  if(!is.null(segdata_obs)){
    g_graph <- g_graph +
      geom_point(data = segdata_obs,
                 aes(x = X, y = Y, size = log1p(y)), col = alpha("red",0.3))
  }

  return(g_graph)

}
