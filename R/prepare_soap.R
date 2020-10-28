
autocruncher <- function(bnd,knots,nmax=200,k=10, xname="x", yname="y") {
  ## setup soap film  smooth - nmax is number of grid cells for longest side
  ## it's important that grid cells are square!

  autocrunch.knots <- function(G,knots,x0,y0,dx,dy){
    ## finds indices of knot locations in solution grid
    ## the knot x,y locations are given in the `knots' argument.
    badk <- NULL
    nk <- length(knots$x)
    ki <- rep(0,nk)
    nx <- ncol(G);ny <- nrow(G)
    if (nk==0) return(NULL)
    for (k in 1:nk) {
      i <- round((knots$x[k]-x0)/dx)+1
      j <- round((knots$y[k]-y0)/dy)+1
      if (i>1&&i<=nx&&j>1&&j<=ny) {
        ki[k] <- G[j,i]
        if (ki[k] <= 0) {
          badk <- c(badk, k)
          #str <- paste("knot",k,"is on or outside boundary")
          #stop(str)
        }
      }
    } ## all knots done
    #ki ## ki[k] indexes kth knot in solution grid
    badk
  } ## end crunch.knots


  # boundary names must be x and y!
  bnd <- lapply(bnd, function(x, xname, yname){
    if(!all(names(x) == c("x", "y"))){
      x$x <- x[[xname]]
      x[[xname]] <- NULL
      x$y <- x[[yname]]
      x[[yname]] <- NULL
    }
    x
  }, xname=xname, yname=yname)
  knots$x <- knots[[xname]]
  knots[[xname]] <- NULL
  knots$y <- knots[[yname]]
  knots[[yname]] <- NULL

  ## check boundary...
  if (!inherits(bnd,"list")) stop("bnd must be a list.")
  n.loops <- length(bnd)
  if (n.loops!=length(k)) {
    if (length(k)==1) k <- rep(k,n.loops)
    else stop("lengths of k and bnd are not compatible.")
  }
  bnd <- mgcv:::process.boundary(bnd) ## add distances and close any open loops

  ## create grid on which to solve Laplace equation
  ## Obtain grid limits from boundary 'bnd'....
  x0 <- min(bnd[[1]]$x);x1 <- max(bnd[[1]]$x)
  y0 <- min(bnd[[1]]$y);y1 <- max(bnd[[1]]$y)
  if (length(bnd)>1) for (i in 2:length(bnd)) {
    x0 <- min(c(x0,bnd[[i]]$x)); x1 <- max(c(x1,bnd[[i]]$x))
    y0 <- min(c(y0,bnd[[i]]$y)); y1 <- max(c(y1,bnd[[i]]$y))
  } ## now got the grid limits, can set it up

  if (x1-x0>y1-y0) { ## x is longest side
    dy <- dx <- (x1-x0) /(nmax-1)
    nx <- nmax
    ny <- ceiling((y1-y0)/dy)+1
  } else { ## y is longest side
    dy <- dx <- (y1-y0) /(nmax-1)
    ny <- nmax
    nx <- ceiling((x1-x0)/dy)+1
  }
  ## so grid is now nx by ny, cell size is dx by dy (but dx=dy)
  ## x0, y0 is "lower left" cell centre

  ## Create grid index G
  bnc <- mgcv:::bnd2C(bnd) ## convert boundary to form required in C code

  G <- matrix(0,ny,nx)
  nb <- rep(0,bnc$n.loop)

  oo <- .C(mgcv:::C_boundary,G=as.integer(G), d=as.double(G), dto=as.double(G), x0=as.double(x0),
           y0 = as.double(y0), dx=as.double(dx), dy = as.double(dy),
           nx=as.integer(nx),as.integer(ny), x=as.double(bnc$x),y=as.double(bnc$y),
           breakCode=as.double(bnc$breakCode),n=as.integer(bnc$n),nb=as.integer(nb))

  ret <- list(G=matrix(oo$G,ny,nx),nb=oo$nb,d=oo$d[oo$d >= 0],x0=x0,y0=y0,dx=dx,dy=dy,bnd=bnd)
  rm(oo)
  ## Now create the PDE coefficient matrix
  #  n.inside <- sum(ret$G > - nx*ny)
  #  xx <- rep(0,5*n.inside)
  #  o1 <- .C(C_pde_coeffs,as.integer(ret$G),xx=as.double(xx),ii=as.integer(xx),jj=as.integer(xx),
  #            n=as.integer(0),as.integer(nx),as.integer(ny),as.double(dx),as.double(dy))
  #  ind <- 1:o1$n
  #  X <- sparseMatrix(i=o1$ii[ind]+1,j=o1$jj[ind]+1,x=o1$xx[ind])
  #  er <- expand(lu(X))
  #  ret$Q <- er$Q;ret$U <- er$U;ret$L <- er$L;ret$P <- er$P
  #  ret$ng <- n.inside ## the number of cells to solve for
  #  rm(er);rm(X)
  ## ... so the sparse LU decomposition of X can be used to solve PDE.
  ## X = PLUQ where P and Q are permuation matrices.

  ## now obtain location of knots in solution ...
  ret <- autocrunch.knots(ret$G,knots,x0,y0,dx,dy)

  ret
} ## end of setup.soap



#' Prepare soap for dsm modelling
#'
#' Create necessary output for dsm modelling using soap smoother.
#'
#' @param data data.frame containing at least 3 columns. Among those columns it must have \code{longitude} and
#'  \code{latitude}.
#' @param contour Land contour of the study area where the soap will be implemented (optionnal).
#'  Must be a Spatial or sf object containing polygons defining land contour (default is NULL).
#' @param polygon_pred Polygon defining the area where the soap will be implemented.
#' Either a Spatial or sf object containing polygon.
#' @param ratio_simplify Correspond to the 'keep' argument of \code{\link[rmapshaper]{ms_simplify}} (default 1e-4).
#' @param N Number of knots per row and column (default 23).
#'
#' @return \enumerate{
#'           \item sf_cropped_NEA : \code{data.frame} (given if \code{polygon_pred} is \code{NULL}).
#'           \item bnd : \code{list} of \code{list} delimiting soap area. Compatible with \code{\link[mgcv]{gam}}.
#'           \item knots : \code{data.frame} containing longitude and latitude of konts.
#'           \item df_cropped_data : \code{data.frame} corresponding to \code{data} without the points outside of \code{bnd}.
#'           \item gg_cropped_data : \code{gg} object that synthetised information on the soap process which has been made.
#'           \item sf_ocean : \code{sf} object containging the contour of the study area where knots are placed.
#'         }
#' @export
#'
#' @examples
#' @references This function uses \code{autocruncher}, a function created by Simon Wood and improved by David L. Miller
#' (\link[gtihub link]{https://github.com/dill/soap_checker/blob/master/autocrunch.R}).
prepare_soap <- function(data,
                         contour = NULL,
                         polygon_pred = NULL,
                         ratio_simplify = 0.00010,
                         N = 23) {

  if(!all(c("longitude","latitude") %in% colnames(data))){
    stop("data must contain 'longitude' and 'latitude' columns")
  }

  if(ncol(data) < 3) {
    stop("data must contain at least one column in addition of columns 'longitude' and 'latitude'")
  }

  cat("longitude and latitude must be projected in WGS84\n")

  # découper NEA à l'échelle de la zone d'étude
  x <- range(data$longitude) * c(0.999,1.001)
  y <- range(data$latitude) * c(0.999,1.001)

  bbox <- c(xmin = x[1], xmax = x[2], ymin = y[1], ymax = y[2])

  if(is.null(polygon_pred)) {

    # force NEA to be a particular file
    if(!is.null(contour)) {

      # contour should be sf or Spatial
      if(!any(is(contour) %in% c("Spatial","sf"))) {
        stop("contour should be either a 'Spatial' or a 'sf' object")
      }

      if(any(is(contour) != "sf")) {
        contour %<>%
          st_as_sf()
      }

      # contour should be crs = 4326
      if(st_crs(contour)[[1]] != 4326 | is.na(st_crs(contour)[[1]])) {
        contour %<>%
          st_transform(crs = 4326)
      }

      land <- contour

    } else {

      land <- NEA %>%
        st_as_sf()

    }


    sf_cropped_NEA <- land %>%
      st_crop(st_bbox(bbox)) %>%
      ms_simplify(keep = ratio_simplify) %>%
      st_as_sf() %>%
      st_union()


    # get / create polygon where to put prediction points ---------------------

    rect_coords <- as.data.frame(rbind(c(x[1],y[1]), c(x[1],y[2]), c(x[2],y[2]), c(x[2],y[1]), c(x[1],y[1])))
    colnames(rect_coords) <- c("x","y")

    sf_polygon <- rect_coords %>%
      st_as_sf(coords = c("x", "y"), crs = 4326) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_cast("POLYGON")

    # diff entre rect et cropped_NEA
    sf_ocean <- st_difference(sf_polygon, sf_cropped_NEA)
    ocean <- sf_ocean %>% st_as_sf(crs = 4326) %>%
      as_Spatial()

  }


  if(!is.null(polygon_pred)){

    if(!any(is(polygon_pred) == "sf")) {
      stop("polygon_pred should be a sf object")
    }

    sf_ocean <- polygon_pred %>%
      st_buffer(0) %>%
      ms_simplify(keep = ratio_simplify) %>%
      st_union() %>%
      st_transform(4326)

  }

  ocean <- sf_ocean %>%
    as_Spatial()

  ocean_unlisted <- ocean %>%
    broom::tidy() %>%
    dplyr::select(long,lat,piece)

  names(ocean_unlisted) <- c("longitude", "latitude", "piece")
  borderlist <- split(ocean_unlisted, ocean_unlisted$piece)
  names(borderlist)
  border.aut <-lapply(borderlist,`[`,c(1,2))
  nr <-seq(1,max(as.numeric(names(borderlist))))
  border.aut <-lapply(sort(as.character(nr)), function(n)as.list.data.frame(border.aut[[n]]))
  lapply(border.aut, function(x) {nrow(as.data.frame(x))})

  # N <- 10
  gx <- seq(min(x) + 0.05*(max(x)-min(x)), max(x), len = N)
  gy <- seq(min(y) + 0.05*(max(y)-min(y)), max(y), len = N)
  gp <- expand.grid(gx, gy)
  names(gp) <- c("longitude","latitude")

  knots <- gp[with(gp, mgcv::inSide(bnd = border.aut, longitude, latitude)), ]
  names(knots) <-c("longitude", "latitude")

  # delete knots to close to border (Function of D. Miller)
  crunch_ind <- autocruncher(border.aut, knots, k=56, xname="longitude", yname="latitude")

  # crop data and convert in data.frame
  data <- data %>%
    st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
    as_Spatial()

  cropped_data <- raster::crop(data, ocean)
  df_cropped_data <- data.frame(cropped_data)
  df_cropped_data %<>%
    rename(longitude = coords.x1,
           latitude = coords.x2)

  # get lines of data that are deleted during crop
  sf_data <- data %>%
    st_as_sf(coords = c("longitude","latitude"), crs = 4326)

  sf_deleted_points <- st_difference(sf_data, sf_ocean)

  # plot ocean, knots, df_cropped_data
  gg_cropped_data <- ggplot() +
    geom_sf(data = sf_ocean) +
    geom_point(data = knots, aes(x = longitude, y = latitude), colour = alpha('red',1)) +
    geom_point(data = knots[crunch_ind, ], aes(x = longitude, y = latitude), colour = alpha("green",0.5)) +
    geom_sf(st_as_sf(df_cropped_data,
                     coords = c("longitude","latitude"),
                     crs = 4326),
            mapping = aes(),
            colour = alpha("blue", 0.5),
            fill = alpha("blue", 0.5),
            shape = 22) +
    geom_sf(data = sf_deleted_points, mapping = aes(), colour = "black") +
    ggtitle(paste0("There were ",nrow(data)-nrow(df_cropped_data)," obs deleted during the crop")) +
    theme_light()


  if(!is.null(polygon_pred)){

    return(list(bnd = border.aut,
                knots = knots[-crunch_ind,],
                df_cropped_data = df_cropped_data,
                gg_cropped_data = gg_cropped_data,
                sf_ocean = sf_ocean))

  } else {

    return(list(sf_cropped_NEA = sf_cropped_NEA,
                bnd = border.aut,
                knots = knots[-crunch_ind,],
                df_cropped_data = df_cropped_data,
                gg_cropped_data = gg_cropped_data,
                sf_ocean = sf_ocean))

  }

}


