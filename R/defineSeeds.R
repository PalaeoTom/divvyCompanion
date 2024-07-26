#' Define centroids of radially constrained regions to be subsampled by [cookies()]
#'
#' @param dat A data.frame with at least 2 columns. These 2 columns must contain the Cartesian coordinates of the seeds to be used by [cookies()]. Their names must be submitted via argument 'xy' as a character vector.
#' @param xy A character vector with two elements, specifying the names of the columns in 'dat' containing the Cartesian coordinates of the seeds.
#'
#' @return A data frame.
#'
#' @export
#'
#' @description This function defines seed points for use with [cookies2()] (via argument `seeding`).
#'
#' @examples
#' # First, load terra
#' library(terra)
#' # Generate occurrence data
#' n <- 100
#' set.seed(5)
#' # 100 sets of x and y coordinates in meters
#' x <- runif(n, 0, 50000)
#' y <- runif(n, 0, 50000)
#' # Set name for coords
#' xy.name <- c("x", "y")
#' # Combine into data frame and label columns
#' pts <- data.frame(x, y)
#' colnames(pts) <- xy.name
#' # Toy dataset ready
#' # Set coordinates for our seeds
#' seedX <- c(10000, 15000, 20000, 25000, 30000, 35000, 40000)
#' seedY <- c(10000, 15000, 20000, 25000, 30000, 35000, 40000)
#' # Concatenate into data.frame with
#' # same coordinate names
#' dat <- data.frame(seedX, seedY)
#' colnames(dat) <- xy.name
#' # Define seedMatrix
#' seedMatrix <- defineSeeds(dat = dat, xy = xy.name)
#' # Get subsamples
#' subsamples <- cookies(dat = pts, xy = xy.name, r = 10000,
#' seeding = seedMatrix, output = "full")
defineSeeds <- function(dat, xy = c("x","y")){
    if(all(xy %in% colnames(dat))){
      ## prepare output
      seedMatrix <- dat
      ## create id column
      seedMatrix[,"id"] <- paste0("seed", 1:nrow(seedMatrix))
      ## assign seedMatrix subclass
      attr(seedMatrix, "class") <- c("data.frame", "seedMatrix")
      ## return seedMatrix
      return(seedMatrix)
    } else {
      stop('column names specified in "xy" are not all present in column names of "dat"')
    }
}
