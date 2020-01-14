#' @export

env_graph <- function(predata, varname = "Depth") {
  g_env_graph <- ggplot() +
    geom_tile(data = predata, aes_string(x = "X", y = "Y", fill = varname)) +
    theme_map() +
    scale_fill_gradientn(name = varname,
                         colours = viridisLite::viridis(256)
    ) +
    guides(fill = guide_colorbar(barwidth = 8, barheight = 0.6,
                                 title.position = "top",
                                 title.hjust = 0.5),
           size = guide_legend(title.position = "left")
    ) +
    xlab("Eastings") + ylab("Northings") +
    theme(legend.position = "top",
          plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 16)
    )
  return(g_env_graph)
}
