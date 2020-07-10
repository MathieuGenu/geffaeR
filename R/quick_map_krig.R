#' @export

quick_map_krig <- function(df, df_proj, varname, lat, lon, NEA, lbrt93_proj, segdata_obs) {

  # kriging graph
  NEA_lbrt93 <- spTransform(NEA, lbrt93_proj)
  europe_sf <- st_as_sf(NEA_lbrt93)

  # get graph window values
  range_x <- range(df[,lon])
  range_y <- range(df[,lat])

  x_lim  <- range_x + c((0.05 * diff(range_x))*c(-1,1))
  y_lim  <- range_y + c((0.05 * diff(range_y))*c(-1,1))


  g_graph <- ggplot() +
    geom_tile(data = df_proj,
              aes(x = df_proj[,lon], y = df_proj[,lat], fill = df_proj[,varname])) +
    geom_point(data = segdata_obs,
               aes(x = X, y = Y, size = log1p(y)), col = alpha("red",0.3)) +
    geom_sf(data = europe_sf) +
    coord_sf(xlim = x_lim, ylim = y_lim, expand = FALSE) +
    scale_fill_gradientn(colours = viridisLite::viridis(256),
                         name = varname
    ) +
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
          axis.title.y = element_blank()
    )+
    labs(size = "ln(nb_ind+1)")

}
