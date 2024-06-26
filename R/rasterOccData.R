#' Rasterise occurrence data and create seed matrix for cookies2
#'
#' @param occData D
#' @param res D
#' @param xyCoords1 D
#' @param xyCoords2 D
#' @param occData.crs D
#' @param raster.crs D
#' @param xyCoords D
#' @param xyCell D
#'
#' @return D
#' @export
#'
#' @examples
#' print("example")
rasteriseOccData <- function(occData, res, xyCoords1 = NULL, xyCoords2 = NULL, occData.crs = 'EPSG:4326', raster.crs = 'EPSG:8857', xyCoords = c('paleolng','paleolat'), xyCell = c('cellX','cellY')){
  if(is.null(xyCoords1) && is.null(xyCoords2)){
    ## initialise
    rWorld <- terra::rast()
    ## define resolution and coordinate system of raster
    rPrj <- terra::project(x = rWorld, y = raster.crs, res = res)
    values(rPrj) <- 1:terra::ncell(rPrj)
    ## convert occurrence data to spatVector
    llOccs <- terra::vect(occData, geom = xyCoords, crs = occData.crs)
    ## change to new coordinate reference system
    prjOccs <- terra::project(llOccs, raster.crs)
    ## get cell numbers for each occurrence
    occData$cell <- terra::cells(rPrj, prjOccs)[,'cell']
    ## get coordinates of each cell
    occData[, xyCell] <- terra::xyFromCell(rPrj, occData$cell)
    ## return rasterised occurrence data
    return(occData)
  } else {
    if(is.numeric(xyCoords1) && is.vector(xyCoords1) && is.numeric(xyCoords2) && is.vector(xyCoords2)){
      ## initialise
      rWorld <- terra::rast()
      ## define resolution and coordinate system of raster
      rPrj <- terra::project(x = rWorld, y = raster.crs, res = res)
      values(rPrj) <- 1:terra::ncell(rPrj)
      ## create mock matrix
      seedMatrix <- as.data.frame(matrix("seed", ncol = ncol(occData), nrow = length(xyCoords1)))
      colnames(seedMatrix) <- colnames(occData)
      seedMatrix[,xyCoords[1]] <- xyCoords1
      seedMatrix[,xyCoords[2]] <- xyCoords2
      ## convert seedMatrix to spatVector
      llOccs <- terra::vect(seedMatrix, geom = xyCoords, crs = occData.crs)
      ## change to desired coordinate system
      prjOccs <- terra::project(llOccs, raster.crs)
      ## get cell numbers
      seedMatrix$cell <- terra::cells(rPrj, prjOccs)[,'cell']
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
