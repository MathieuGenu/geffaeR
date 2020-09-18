#' Creation of predata and data impuation of missing data.
#'
#' Creation of predata, a data.frame containing environmental parameters in a grid of the study area.
#' Also it allows to impute missing data for predata and segdata.
#'
#' @inheritParams prepare_data_obs
#' @param gridfile_name Grid of the study area with all the environmental variables which
#' allows to build predata.
#' @param varenviro Vector containing all dynamic environmental variables (e.g. chlorophyl, SST).
#' @param do_log_enviro Vector precising which variables among varenviro that are necessary
#' to transform to neperian logarithm + 1 (\code{\link[base]{log1p}}).
#' @param varphysio Vector containing non-dynamic environmental variables (i.e depth, dist200m).
#' @param do_log_physio Vector precising which variables among varphysio that are necessary
#' to transform to neperian logarithm + 1 (\code{\link[base]{log1p}}).
#' @param imputation Data imputation method to missing values of segdata and predata. 2 methods allowed :
#'        \itemize{
#'          \item "PCA" : PCA method from \pkg{missMDA}.
#'          \item "Amelia" : Amelia method from \pkg{Amelia}.
#'        }
#'        Default method is PCA.
#' @param saturate_predata Boolean. If \code{TRUE}, saturate function is applied on all varenviro
#'  an varphysio columns of segdata. Saturate function excludes extreme valuesand keep values between
#'  quantiles 95 percent and 5 percent.
#' @param saturate_segdata Boolean. If \code{TRUE}, saturate function is applied on all varenviro
#'  an varphysio columns of predata. Saturate function excludes extreme valuesand keep values between
#'  quantiles 95 percent and 5 percent.
#' @param col2keep \code{character string} corresponding to the columns wanted to appear in output in predata.
#' @param inbox_poly When \code{TRUE}, keep only part of the grid that match with the study area.
#' @return This function return a list containing :
#'         \enumerate{
#'           \item predata : data.frame of predata.
#'           \item segdata : data.frame of segdata.
#'           \item pca_pred : Output of \code{\link[FactoMineR]{PCA}} function on segdata.
#'           \item pca_seg : Output of \code{\link[FactoMineR]{PCA}} function on predata.
#'         }
#' @import dplyr
#' @importFrom Amelia amelia
#' @importFrom FactoMineR PCA
#' @importFrom foreign read.dbf
#' @importFrom missMDA estim_ncpPCA imputePCA
#' @importFrom sf st_read st_as_sf st_cast st_crs st_intersects st_transform
#' @importFrom VIM aggr
#' @examples
#'
#' @export


prep_predata <- function(segdata,
                         gridfile_name,
                         varenviro, do_log_enviro,
                         varphysio, do_log_physio,
                         imputation = "Amelia",
                         shape,
                         saturate_predata = F, saturate_segdata = F,
                         inbox_poly = T,
                         col2keep = NULL) {

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
  # shape can be character dir or object directly
  if(any("character" %in% is(shape))){
    # pred.poly <- readOGR(dsn = paste(shape), layer = layer, verbose = F) # NC
    pred.poly <- st_read(paste(shape)) %>%
      mutate(Id = 1:n())
  } else {
    pred.poly <- shape %>% st_as_sf()
  }

  ### Covariable(s)
  # grille de la zone d'étude
  # grid can be character dir or object directly
  if(any("character" %in% is(gridfile_name))){
    grid <-  read.dbf(paste(gridfile_name, sep = "/"), as.is = TRUE)
  } else {
    grid <- gridfile_name
  }

  grid <- grid[which(duplicated(grid) == FALSE), ]

  # check if col2keep are in colnames of grid
  if(!all(col2keep %in% colnames(grid))){
    miss_col2keep <- col2keep[!(col2keep %in% colnames(grid))]
    stop(paste(miss_col2keep, "n'est pas/ ne sont pas, dans grid",collapse = "\n"))
  }

  # colname of coord in grid in good format
  if(all(c("lat","lon") %in% colnames(grid))) {
    colnames(grid)[colnames(grid) %in% 'lat'] <- "LATITUDE"
    colnames(grid)[colnames(grid) %in% 'lon'] <- "LONGITUDE"
  }

  ## MA ## Ne conserver que les valeurs dans la zone de prédiction
  grid_xy <- grid[, c("LONGITUDE", "LATITUDE")] %>%
    st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
    st_cast("POINT")
  # coordinates(grid_xy) <- ~ LONGITUDE + LATITUDE
  # grid_xy@proj4string <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  if(st_crs(grid_xy) != st_crs(pred.poly)) {
    grid_xy <- grid_sf %>% st_transform(crs = st_crs(pred.poly))
  }
  # if(proj4string(grid_xy) != proj4string(pred.poly)) {
  #   grid_xy <- spTransform(grid_xy, proj4string(pred.poly))
  # }

  xy_inside <- st_intersects(grid_xy, pred.poly, sparse = FALSE)#which(!is.na(over(grid_xy, pred.poly)[,"Id"]))
  grid$inbox <- apply(xy_inside, 1, any)
  # grid$inbox <- NA
  # grid$inbox[xy_inside] <- 1
  # grid$inbox[!xy_inside] <- 0
  if(inbox_poly == T) {
    grid <- subset(grid, inbox == TRUE)
  }

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
  predata <- data.frame(grid[, c("LONGITUDE", "LATITUDE", "POINT_X", "POINT_Y", "Area", col2keep,
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


