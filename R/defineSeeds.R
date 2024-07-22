#' Define centroids of radially constrained regions to be subsampled by [cookies()]
#'
#' @param grid A SpatVector defining the gridding template to be used. Note, if terra::values(grid) = NA for all cells, the grid cells defined by the template will be numbered, starting from 1, and the resulting data frame of seeds will use this numbering scheme. It is **highly recommended** that users label their gridding template SpatVectors before generating their seed matrices.
#' @param dat A data.frame with at least 2 columns. These 2 columns must contain the Cartesian coordinates of the seeds to be used by [cookies()]. Their names must be submitted via argument 'xy' as a character vector.
#' @param xy A character vector with two elements, specifying the names of the columns in 'dat' containing the Cartesian coordinates of the seeds.
#'
#' @return A data frame.
#'
#' @export
#'
#' @description This function defines rasterised seed points for use with [cookies2()] (via argument `seeding`).
#'
#' It is **highly recommended** that you first rasterise your occurrence data, and then use the same grid cell template you used to do so to define your seed matrix **as this will ensure your seeds and point occurrences are rasterised using the same grid of cells**.
#'
#' @examples
#' # First, load terra
#' library(terra)
#' # Work in Equal Earth project coordinates
#' prj <- 'EPSG:8857'
#' # Generate occurrence data
#' n <- 100
#' set.seed(5)
#' # 100 sets of x and y coordinates
#' x <- runif(n, 0, 50)
#' y <- runif(n, 0, 50)
#' # Units in equal earth are in meters
#' # So if we consider x and y as given in km
#' x <- x * 1000
#' y <- y * 1000
#' # Set name for coords
#' xy.name <- c("x", "y")
#' # Combine into data frame and label columns
#' pts <- data.frame(x, y)
#' colnames(pts) <- xy.name
#' # Finally, convert to SpatVector
#' ptsPrj <- terra::vect(pts, geom = xy.name,
#'                      crs = prj)
#' # Now let's define grid template with 5km grid cells
#' grid <- terra::project(x = terra::rast(),
#'                        y = prj, res = 5000)
#' # Number grid cells
#' terra::values(grid) <- 1:terra::ncell(grid)
#' # Now assign occurrences in toy dataset to grid cells
#' # First the cell number
#' pts[,"cell"] <- terra::cells(grid, ptsPrj)[,"cell"]
#' # Then the cell coordinates
#' pts[,c("cellX", "cellY")] <- terra::xyFromCell(grid,
#'                                                pts[,"cell"])
#' # Toy dataset ready
#' # Set coordinates for our seeds
#' # in km
#' seedX <- c(10000, 15000, 20000, 25000, 30000, 35000, 40000)
#' seedY <- c(10000, 15000, 20000, 25000, 30000, 35000, 40000)
#' # Concatenate into data.frame
#' dat <- data.frame(seedX, seedY)
#' colnames(dat) <- xy.name
#' # Use defineSeeds to create seed matrix
#' seedMatrix <- defineSeeds(grid = grid,
#'                          dat = dat, xy = xy.name)
#' # Use cookies to assess viability of 10km radius regions
#' # about seeds. Only viable regions returned.
#' subsamples <- cookies(dat = pts, xy = c("cellX", "cellY"),
#'                      seeding = seedMatrix,
#'                      r = 10000, rarefy = FALSE,
#'                      output = "full")
#' length(subsamples)
#' # Just two of the seven seed points are viable!
defineSeeds <- function(grid, dat, xy = c("x","y")){
  if(inherits(grid, "SpatRaster")){
    if(all(xy %in% colnames(dat))){
      ## get crs EPSG code
      prj <- paste(unlist(terra::crs(grid, describe = T)[2:3]), collapse = ":")
      if(all(is.na(suppressWarnings(terra::values(grid))))){
        terra::values(grid) <- 1:terra::ncell(grid)
      }
      ## prepare output
      seedMatrix <- dat
      ## convert dat to spatVector
      sv <- terra::vect(dat, geom = xy, crs = prj)
      ## add cell numbers
      seedMatrix[,"cell"] <- terra::cells(grid, sv)[,"cell"]
      ## get coordinates for each cell
      seedMatrix[,c("cellX","cellY")] <- terra::xyFromCell(grid, seedMatrix$cell)
      ## create id column
      seedMatrix[,"id"] <- paste0("seed", 1:nrow(seedMatrix))
      ## return seedMatrix
      return(seedMatrix)
    } else {
      stop('column names specified in "xy" are not all present in column names of "dat"')
    }
  } else {
    stop('"grid" is not a SpatRaster')
  }
}
