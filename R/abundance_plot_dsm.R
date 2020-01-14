# graphique des abondances

#' @export

abundance_plot <- function(predata, stat = "mean", inside = NULL, abundance = TRUE, saturation = NULL) {

  if(stat == "mean") {
    varname <- "pred" ; titre <- "Densité (ind/km^2)"
    if(any(names(predata) == "Area")) { predata$pred <- predata$pred/predata$Area }
  }
  else{ varname <- "cv" ; titre <- "CV (%)" }

  ### enlever les predictions en dehors d'une boite
  if(is.null(inside)) { predata[which(predata[, "inside"] == 0), varname] <- NA }
  else { predata[-which(predata[, "zone"] %in% inside), varname] <- NA }

  if(!is.null(saturation)) {
    if(!is.numeric(saturation)) {
      writeLines("Attention valeur non numérique fournie: par défaut saturation au quantile 99%")
      saturation <- quantile(predata$pred[!is.na(predata$pred)], prob = 0.99)
    }
    else{
      if(saturation >= 1) {
        writeLines("Attention la valeur fournie n'est pas un quantile: par défaut saturation au quantile 99%")
        saturation <- quantile(predata$pred[!is.na(predata$pred)], prob = 0.99)
      }
      else{
        saturation <- quantile(predata$pred[!is.na(predata$pred)], prob = saturation)
      }
    }
    predata$pred <- ifelse(predata$pred > saturation, saturation, predata$pred)
  }

  predata$cv <- 100*round(predata$se/predata$pred, 3)
  predata$cv <- ifelse(predata$cv > 100, 100, predata$cv)

  ### enlever les predictions en dehors d'une boite
  if(is.null(inside)) { predata[which(predata[, "inside"] == 0), varname] <- NA }
  else { predata[-which(predata[, "zone"] %in% inside), varname] <- NA }   # A garder ou pas ?

  g_abundance_plot <- ggplot() +
    geom_tile(data = predata, aes_string(x = "X", y = "Y", fill = varname)) +
    coord_equal() + theme_map() +
    scale_fill_gradientn(name = titre,
                         colours = viridisLite::viridis(256)
    ) +
    guides(fill = guide_colorbar(barwidth = 8, barheight = 0.6,
                                 title.position = "top",
                                 title.hjust = 0.5),
           size = guide_legend(title.position = "left")
    ) +
    ylab("Northings") + xlab("Eastings") +
    theme(legend.position = "top",
          plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 16)
    )
  return(g_abundance_plot)
}
