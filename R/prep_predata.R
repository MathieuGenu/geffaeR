#' \encoding{Création de predata et imputation de données manquantes}
#'
#' \encoding{Cette fonction permet de créer predata, le data.frame contenant les données environnementales dans la
#' grille de la zone d'étude. Cette fonction permet également l'imputation de données manquantes pour les data.frames
#' segdata et predata}
#'
#' @inheritParams prepare_data_obs
#' @param gridfile_name \encoding{Grille de la zone d'étude avec toutes les variables d'intérêt qui permettront la fabrication
#'        de predata.}
#' @param proj_grid \encoding{Systême de projection du fichier de gride.}
#' @param varenviro \encoding{Vecteur contenant les variables environnementales dynamiques (i.e chlorophylles, SST).}
#' @param do_log_enviro \encoding{vecteur précisant parmi les variables environnementales (varenviro), lesquelles
#'        sont à transformer en logarithme néperien + 1 (\code{\link[base]{log1p}}).}
#' @param varphysio \encoding{Vecteur contenant les variables environnementales non-dynamiques (i.e depth, dist200m).}
#' @param do_log_physio \encoding{vecteur précisant parmi les variables physiologiques (varphysio), lesquelles
#'        sont à transformer en logarithme néperien + 1 (\code{\link[base]{log1p}}).}
#' @param imputation \encoding{Méthode permettant l'imputation de données manquantes pour segdata et predata.
#'        Parmi les méthodes il y a :
#'        \itemize{
#'          \item{"PCA" : }{\encoding{par méthode ACP. méthode issue du package \pkg{missMDA}}.}
#'          \item{"Amelia" : }{\encoding{par méthode Amelia. méthode issue du pakage \pkg{Amelia}}.}
#'        }
#'        Par défaut la méthode est "PCA".}
#' @param shape_pred \encoding{shapefile délimitant la zone où la prédiction est effectuée.}
#' @param layer_pred \encoding{Layer du shpaefile donnant les limite de la zone de prédiction.}
#' @param saturate_predata \encoding{Booléen. Si \code{TRUE}, la fonction saturate est appliquée à segdata sur
#'        toutes les colonnes varenviro et varphysio. La fonction saturate va exclure les valeurs extrêmes de ces
#'        colonnes en conservant uniquement les quantiles à 95% et 5%.}
#' @param saturate_segdata \encoding{Booléen. Si \code{TRUE}, la fonction saturate est appliquée à predata sur
#'        toutes les colonnes varenviro et varphysio. La fonction saturate va exclure les valeurs extrêmes de ces
#'        colonnes en conservant uniquement les quantiles à 95% et 5%.}
#' @return Cette fonction renvoit une liste contenant :
#'         \enumerate{
#'           \item predata : Un data.frame contenant predata
#'           \item segdata : Un data.frame contenant segdata
#'           \item pca_pred : Une sortie de la fonction \code{\link[FactoMineR]{PCA} sur le data.frame segdata}
#'           \item pca_seg : Une sortie de la fonction \code{\link[FactoMineR]{PCA} sur le data.frame predata}
#'         }
#'
#' @examples
#'
#' @export


prep_predata <- function(segdata,
                         gridfile_name,
                         varenviro, do_log_enviro,
                         varphysio, do_log_physio,
                         imputation = "Amelia",
                         shape_pred,layer_pred,
                         saturate_predata = F, saturate_segdata = F,
                         inbox_poly = T) {

  # check si toutes les covariables sont dans segdata
  if(!all(varphysio %in% colnames(segdata))) {
    miss_physio <- varphysio[!(varphysio %in% colnames(segdata))]
    stop(paste(miss_physio, "n'est pas/ ne sont pas, dans segdata",collapse = "\n"))
  }
  if(!all(varenviro %in% colnames(segdata))) {
    miss_enviro <- varenviro[!(varenviro %in% colnames(segdata))]
    stop(paste(miss_enviro, "n'est pas/ ne sont pas, dans segdata",collapse = "\n"))
  }

  # dire à l'utilisateur qu'il faut bien que les données manquantes soient au format NA et pas autre chose
  message("Les données manquantes dans les tableaux grid et segdata doivent impérativement être sous la forme de NA
          et rien d'autre (pas de 0)")

  segdata_old <- segdata
  segdata <- NULL

  ## prediction
  pred.poly <- readOGR(dsn = paste(shape_pred, sep = "/"), layer = layer_pred)

  ### Covariable(s)
  # grille de la zone d'étude pour 2017
  grid <-  read.dbf(paste(gridfile_name, sep = "/"), as.is = TRUE)
  grid <- grid[which(duplicated(grid) == FALSE), ]

  # colname of coord in grid in good format
  if(all(c("lat","lon") %in% colnames(grid))) {
    colnames(grid)[colnames(grid) %in% 'lat'] <- "LATITUDE"
    colnames(grid)[colnames(grid) %in% 'lon'] <- "LONGITUDE"
  }

  ## MA ## Ne conserver que les valeurs dans la zone de prédiction
  grid_xy <- grid[, c("LONGITUDE", "LATITUDE")]
  coordinates(grid_xy) <- ~ LONGITUDE + LATITUDE
  grid_xy@proj4string <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  if(proj4string(grid_xy) != proj4string(pred.poly)) {
    grid_xy <- spTransform(grid_xy, proj4string(pred.poly))
  }
  xy_inside <- which(!is.na(over(grid_xy, pred.poly)[,"Id"]))
  grid$inbox <- NA
  grid$inbox[xy_inside] <- 1
  grid$inbox[!xy_inside] <- 0
  if(inbox_poly == T) {
    grid <- subset(grid,inbox==1)
  }

  # # !!! cas particulier pour ANT_GUY !!! #
  # ## SL ## renommer champ correctement, annuler fev pour ne faire prédiction que sur oct
  # grid$CHL_month     <- grid$CHL201710M
  # grid$CHL_monthClim <- grid$CHLm10Clim
  # grid$SST_month     <- grid$SST201710M
  # grid$SST_monthClim <- grid$SSTm10Clim


  ### Covariable
  allvar <- c(varphysio, varenviro)

  if(!all(allvar %in% names(grid))) {
    stop("Les covariables suivantes n'apparaissent pas dans la grille de prediction (grid) :",
         allvar[which(allvar %in% names(grid) == FALSE)], sep = " ")
  }

  ########## création des tableaux nécessaires à l'analyse  ##########
  # recuperer les donnees dans la table d'effort
  segdata <- data.frame(segdata_old[, c("date", "survey", "Transect.Label", "Sample.Label", "Seg", "Effort",
                                    "X", "Y", "longitude", "latitude", "Region.Label","seaState","Area","subjective",
                                    allvar)]
  )

  # à enlever
  names(segdata)[1:10] <- c("date", "survey", "Transect.Label", "Sample.Label","Seg", "Effort", "X",
                            "Y", "longitude", "latitude")


  ## 0 is missing data for varenviro (should be done here to compute accurately the quantiles)
  for(j in varenviro) {
    segdata[, j] <- ifelse(segdata[, j] == 0, NA, segdata[, j])
  }

  # saturer les variables
  # calculer les seuils pour la grille de prediction
  segdata_threshold <- apply(segdata[, allvar], 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  saturate <- function(x, threshold, upper = TRUE) {
    y <- x[which(!is.na(x))]
    if(upper) {
      y[y > threshold] <- threshold
    } else {
      y[y < threshold] <- threshold
    }
    x[which(!is.na(x))] <- y
    return(x)
  }
  if(saturate_segdata == T) {
    for (j in allvar) {
      segdata[, j] <- saturate(segdata[, j], threshold = segdata_threshold[1, j], upper = FALSE)
      segdata[, j] <- saturate(segdata[, j], threshold = segdata_threshold[2, j], upper = TRUE)
    }
  }

  ## impute missing values
  mi_segdata <- segdata[, c("longitude", "latitude", allvar)]

  # echelle log pour la chl et SSTsd si remmoa 2017 --> à changer
  if(!is.null(do_log_enviro)) {
    for(j in do_log_enviro) {
      mi_segdata[, j] <- log1p(mi_segdata[, j])
    }
  }
  if(!is.null(do_log_enviro)) {
    for(j in do_log_physio) {
      mi_segdata[, j] <- log1p(mi_segdata[, j])
    }
  }

  seg_mipat <- NULL

  if(any(apply(mi_segdata, 2, function(j) { any(is.na(j)) }))) {
    ### missingness patterns
    seg_mipat <- summary(VIM::aggr(mi_segdata, sortVar = TRUE, plot = F))
    cat("il y a",sum(seg_mipat$missings$Count),"cases sans valeurs dans segdata","\n","\n", sep=" ")
    if(any(seg_mipat$missings$Count == nrow(mi_segdata))) {
      writeLines(paste("Missing data imputation not done on data for segments because some variables
                       only have missing values.\n Please check",
                       seg_mipat$missings$Variable[which(seg_mipat$missings$Count == nrow(mi_segdata))], sep = ' '))
      seg_pca <- NULL
    } else {
      if(imputation == "PCA") {
        ## missMDA
        segdata[, c("longitude", "latitude", allvar)] <-
          imputePCA(mi_segdata,ncp = estim_ncpPCA(mi_segdata,ncp.max = ncol(mi_segdata),
                                                  method.cv = "Kfold",
                                                  nbsim = 10
          )$ncp)$completeObs
      } else {
        ## amelia
        amelia_segdata <- amelia(mi_segdata, m = 10)
        mi_segdata <- array(NA, dim = c(nrow(mi_segdata), ncol(mi_segdata), 10))
        for(k in 1:10) { mi_segdata[ , , k] <- as.matrix(amelia_segdata$imputations[[k]]) }
        segdata[, c("longitude", "latitude", allvar)] <- as.data.frame(apply(mi_segdata, c(1, 2), mean))
      }

      # enlever echelle log
      if(!is.null(do_log_enviro)) {
        for(j in do_log_enviro) {
          segdata[, j] <- ifelse(exp(segdata[, j]) < 1, 0, exp(segdata[, j]) - 1)
        }
      }
      if(!is.null(do_log_physio)) {
        for(j in do_log_physio) {
          segdata[, j] <- ifelse(exp(segdata[, j]) < 1, 0, exp(segdata[, j]) - 1)
        }
      }

      # re-seuiler pour être cohérent
      if(saturate_segdata == T){
        for (j in allvar) {
          segdata[, j] <- saturate(segdata[, j], threshold = segdata_threshold[1, j], upper = FALSE)
          segdata[, j] <- saturate(segdata[, j], threshold = segdata_threshold[2, j], upper = TRUE)
        }
      }

      ### PCA analyses
      seg_pca <- PCA(segdata[, c("longitude", "latitude", allvar)], graph = FALSE)
    }
  } else {
    amelia_segdata <- NULL
    ### PCA analyses
    seg_pca <- PCA(segdata[, allvar], graph = FALSE)
  }

  ### table : predata
  predata <- data.frame(grid[, c("LONGITUDE", "LATITUDE", "POINT_X", "POINT_Y", "Area",
                                 allvar)])
  names(predata)[1:4] <- c("longitude", "latitude", "X", "Y")
  predata_xy <- predata[, c("longitude", "latitude")]
  coordinates(predata_xy) <- ~ longitude + latitude
  predata_xy@proj4string <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  ## 0 is missing data for varenviro
  for(j in varenviro) {
    predata[, j] <- ifelse(predata[, j] == 0, NA, predata[, j])
  }

  # saturer les valeurs extremes si demandé
  if(saturate_predata == T) {
    predata_threshold <- apply(predata[, allvar], 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    for (j in allvar) {
      predata[, j] <- saturate(predata[, j], threshold = predata_threshold[1, j], upper = FALSE)
      predata[, j] <- saturate(predata[, j], threshold = predata_threshold[2, j], upper = TRUE)
    }
  }

  ## impute missing values
  mi_predata <- predata[, c("longitude", "latitude", allvar)]

  # echelle log pour la chl et SSTsd
  if(!is.null(do_log_enviro)) {
    for(j in do_log_enviro) {
      mi_predata[, j] <- log1p(mi_predata[, j])
    }
  }
  if(!is.null(do_log_physio)) {
    for(j in do_log_physio) {
      mi_predata[, j] <- log1p(mi_predata[, j])
    }
  }

  pred_mipat <- NULL

  if(any(apply(mi_predata, 2, function(j) { any(is.na(j)) }))) {
    ### missingness patterns
    pred_mipat <- summary(VIM::aggr(mi_predata, sortVar = TRUE, plot = F))
    cat("il y a",sum(pred_mipat$missings$Count),"cases sans valeurs dans predata","\n","\n", sep=" ")
    if(any(pred_mipat$missings$Count == nrow(mi_predata))) {
      writeLines(paste("Missing data imputation not done on data for predictions because some
                       variables only have missing values.\n Please check",
                       pred_mipat$missings$Variable[which(pred_mipat$missings$Count == nrow(mi_predata))], sep = ' '))
      ### PCA analyses
      pred_pca <- NULL
    } else {
      if(imputation == "PCA") {
        ## missMDA
        predata[, c("longitude", "latitude", allvar)] <-
          imputePCA(mi_predata,ncp = estim_ncpPCA(mi_predata,ncp.max = ncol(mi_predata),
                                                  method.cv = "Kfold",
                                                  nbsim = 10)$ncp)$completeObs
      } else {
        ## amelia
        amelia_predata <- amelia(mi_predata, m = 10)
        mi_predata <- array(NA, dim = c(nrow(mi_predata), ncol(mi_predata), 10))
        for(k in 1:10) {
          mi_predata[ , , k] <- as.matrix(amelia_predata$imputations[[k]])
        }
        predata[, c("longitude", "latitude", allvar)] <- as.data.frame(apply(mi_predata, c(1, 2), mean))
      }

      # enlever echelle log
      if(!is.null(do_log_enviro)){
        for(j in do_log_enviro) {
          predata[, j] <- ifelse(exp(predata[, j]) < 1, 0, exp(predata[, j]) - 1)
        }
      }
      if(!is.null(do_log_physio)){
        for(j in do_log_physio) {
          predata[, j] <- ifelse(exp(predata[, j]) < 1, 0, exp(predata[, j]) - 1)
        }
      }
      # re-seuiler pour être cohérent
      if(saturate_predata == T) {
        for (j in allvar) {
          predata[, j] <- saturate(predata[, j], threshold = segdata_threshold[1, j], upper = FALSE)
          predata[, j] <- saturate(predata[, j], threshold = segdata_threshold[2, j], upper = TRUE)
        }
      }

      ### PCA analyses
      pred_pca <- PCA(predata[, c("longitude", "latitude", allvar)], graph = FALSE)
    }
  } else {
    amelia_predata <- NULL
    ### PCA analyses
    pred_pca <- PCA(predata[, allvar], graph = FALSE)
  }

  ## rassembler
  return(list(predata = predata,
              segdata = segdata,
              pca_seg = seg_pca,
              pca_pred = pred_pca,
              seg_mipat = seg_mipat,
              pred_mipat = pred_mipat))
}


