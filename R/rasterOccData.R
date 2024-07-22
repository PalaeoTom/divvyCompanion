#' Rasterise occurrence data and create seed matrix for cookies2
#'
#' @param occData A data.frame or matrix containing taxon names, Cartesian coordinates, and any associated variables. All columns should be labelled.
#' @param res A numeric value specifying the resolution (in metres) of the output raster.
#' @param xyCoords1 Either `NULL` (the default) or a numeric vector specifying the first Cartesian coordinate(s) of the centroid(s) of the radially defined region(s) you wish to sample using [cookies2()]. Vector should be same length as that submitted for `xyCoords2`. Coordinates should use coordinate reference system specified using `occData.crs`.
#' @param xyCoords2 Either `NULL` (the default) or a numeric vector specifying the second Cartesian coordinate(s) of the centroid(s) of the radially defined region(s) you wish to sample using [cookies2()]. Vector should be same length as that submitted for `xyCoords1`. Coordinates should use coordinate reference system specified using `occData.crs`.
#' @param occData.crs A character string specifying the coordinate reference system (CRS) of either the coordinates specified in columns `xyCoords` of `occData` (if `xyCoords1` and `xyCoords2` are null) or the coordinates submitted as numeric vectors to `xyCoords1` and `xyCoords2` (if you are defining a seed matrix). Default is `"EPSG:4326"`.
#' @param raster.crs A character string specifying the coordinate reference system (CRS) to be used to rasterise `occData` (if `xyCoords1` and `xyCoords2` are null) or the coordinates submitted as numeric vectors to `xyCoords1` and `xyCoords2` (if you are defining a seed matrix). Default is `"EPSG:8857"`.
#' @param xyCoords A character vector with two elements, specifying the names of columns in `occData` containing Cartesian coordinates or the names to be used for the coordinates submitted as numeric vectors to `xyCoords1` and `xyCoords2` (if you are defining a seed matrix). Default is `xyCoords = c('paleolng','paleolat')`.
#' @param xyCell A character vector with two elements, specifying the names that should be given to the columns in the output containing the grid cell centroid coordinates. Note, occurrences (i.e., rows in `occData` or the seed points specified as numeric vectors via `xyCoords1` and `xyCoords2`) that fall within the same grid cell will have identical values for columns `xyCell` in the output. Default is `xyCell = c('cellX','cellY')`.
#' @param uniqID A character string, specifying the name to be given to the column in the output containing the unique identifying number for each cell. Note, occurrences (i.e., rows in `occData` or the seed points specified as numeric vectors via `xyCoords1` and `xyCoords2`) that fall within the same grid cell will have identical values for column `uniqID` in the output. Default is `uniqID = "cell`. Cannot match other column names provided or `"id"` (the latter is applied when generating a seed matrix).
#'
#' @return A data frame or matrix (matching the class of `occData`).
#'
#' If `xyCoords1` and `xyCoords2` are null, a rasterised version of `occData` will be returned with *r* rows and at least *c* columns, where *r* is the number of rows in `occData` and *c* is the number of columns in `occData`. If not already present in `occData`, columns named `xyCell` and `uniqID` will be added. If they are present, they will be overwritten.
#'
#' If `xyCoords1` and `xyCoords2` are numeric vectors, a seed matrix will be returned with *l* rows and at least *c* columns, where *l* is the length of vector `xyCoords1` (this should match length of `xyCoords2`) and *c* is the number of columns in `occData`. If not already present in `occData`, columns named `xyCoords`, `xyCell`, `uniqID`, and `"id"` (a column denoting the rows as manually defined seed points) will be added.
#'
#' @export
#'
#' @description This is a wrapper function that employs [terra] functions to rasterise occurrence data and define rasterised seed points for use with [cookies2()] (via argument `seeding`).
#'
#' If defining a seed matrix, it is *highly recommended* that you first rasterise your occurrence data using this function, then submit the rasterised occurrence data as `occData`, *keeping all other arguments the same*, when defining your seed points **as this will ensure the coordinates of your seeds are assigned to the same columns as the coordinates of your occurrences** .
#'
#' You can use the following formats to define coordinate reference systems: WKT, PROJ.4 (e.g., `crs = +proj=longlat +datum=WGS84`), or an EPSG code (e.g., `crs = "EPSG:4326"`). But note that the PROJ.4 notation has been deprecated, and you can only use it with the WGS84/NAD83 and NAD27 datums. Other datums are silently ignored.
#'
#' @examples
#' # First, load terra
#' library(terra)
#'
#' # Work in Equal Earth project coordinates
#' prj <- 'EPSG:8857'
#'
#' # Generate occurrence data
#' n <- 100
#' set.seed(5)
#'
#' # 100 sets of x and y coordinates
#' x <- runif(n, 0, 50)
#' y <- runif(n, 0, 50)
#'
#' # Equal Earth is in metres
#' # Convert from km
#' x <- x * 1000
#' y <- y * 1000
#'
#' # Combine into data frame and label columns
#' pts <- data.frame(x, y)
#' colnames(pts) <- c("x", "y")
#'
#' # Rasterise (raster crs is 'EPSG:8857' by default)
#' # using 5km grid cells
#' raster <- rasterOccData(occData = pts, res = 5000,
#' xyCoords = c("x", "y"), occData.crs = prj,
#' xyCell = c("cellX", "cellY"), uniqID = "cell")
#'
#' # How many viable radially constrained regions do we have if we use
#' # 10km radii, require a minimum of 2 grid cells,
#' # and demand no overlap in sites? Let's find out.
#' standard.seeding <- cookies2(dat = raster,
#' rarefaction = "none", seeding = NULL, uniqID = "cell",
#' xy = c("cellX", "cellY"), nSite = 2, r = 10000,
#' oThreshold = 0, oType = "sites", output = "locs")
#' length(standard.seeding)
#'
#' # Now let's use it to generate a seed matrix for use with [cookies2()]
#' # Set coordinates for our seeds
#' # in km
#' seed.x <- c(10000, 15000, 20000, 25000, 30000, 35000, 40000)
#' seed.y <- c(10000, 15000, 20000, 25000, 30000, 35000, 40000)
#'
#' seed.matrix <- rasterOccData(occData = raster, xyCoords1 = seed.x,
#' xyCoords2 = seed.y, res = 5000, xyCoords = c("x", "y"),
#' occData.crs = prj,
#' xyCell = c("cellX", "cellY"), uniqID = "cell")
#'
#' # Now let's see how many of these seeds
#' # produce usable 10km radius regions
#' manual.seeding <- cookies2(dat = raster, rarefaction = "none",
#' seeding = seed.matrix, uniqID = "cell",
#' xy = c("cellX", "cellY"), nSite = 2, r = 10000, oThreshold = 0,
#' oType = "sites", output = "locs")
#' length(manual.seeding)
rasterOccData <- function(occData, res, xyCoords1 = NULL, xyCoords2 = NULL, occData.crs = 'EPSG:4326', raster.crs = 'EPSG:8857', xyCoords = c('paleolng','paleolat'), xyCell = c('cellX','cellY'), uniqID = "cell"){
  if(is.null(xyCoords1) && is.null(xyCoords2)){
    ## define grid
    rPrj <- terra::project(x = terra::rast(), y = raster.crs, res = res)
    terra::values(rPrj) <- 1:terra::ncell(rPrj)
    ## convert occurrence data to spatVector
    llOccs <- terra::vect(occData, geom = xyCoords, crs = occData.crs)
    ## change to new coordinate reference system (if)
    if(!occData.crs == raster.crs){
      prjOccs <- terra::project(llOccs, raster.crs)
    } else {
      prjOccs <- llOccs
    }
    ## get cell numbers for each occurrence
    occData[, uniqID] <- terra::cells(rPrj, prjOccs)[,"cell"]
    ## get coordinates of each cell
    occData[, xyCell] <- terra::xyFromCell(rPrj, occData$cell)
    ## return rasterised occurrence data
    return(occData)
  } else {
    if(is.numeric(xyCoords1) && is.vector(xyCoords1) && is.numeric(xyCoords2) && is.vector(xyCoords2)){
      ## define resolution and coordinate system of raster
      rPrj <- terra::project(x = terra::rast(), y = raster.crs, res = res)
      terra::values(rPrj) <- 1:terra::ncell(rPrj)
      ## create mock matrix
      seedMatrix <- as.data.frame(matrix("seed", ncol = ncol(occData), nrow = length(xyCoords1)))
      colnames(seedMatrix) <- colnames(occData)
      seedMatrix[,xyCoords[1]] <- xyCoords1
      seedMatrix[,xyCoords[2]] <- xyCoords2
      ## convert seedMatrix to spatVector
      llOccs <- terra::vect(seedMatrix, geom = xyCoords, crs = occData.crs)
      ## change to new coordinate reference system (if different)
      if(!occData.crs == raster.crs){
        prjOccs <- terra::project(llOccs, raster.crs)
      } else {
        prjOccs <- llOccs
      }
      ## get cell numbers
      seedMatrix[, uniqID] <- terra::cells(rPrj, prjOccs)[,"cell"]
      ## get coordinates for each cell
      seedMatrix[, xyCell] <- terra::xyFromCell(rPrj, seedMatrix$cell)
      ## create id column
      seedMatrix$id <- paste0("seed", 1:nrow(seedMatrix))
      ## return seedMatrix
      return(seedMatrix)
    } else {
      stop("check xyCoords1 and 2 inputs")
    }
  }
}
