#' Derive radially constrained spatial subsamples from occurrence data
#'
#' @param dat A data.frame or matrix containing taxon names, Cartesian coordinates, and any associated variables.
#' @param xy 	A vector of two elements, specifying the name of columns in `dat` containing Cartesian coordinates, e.g. longitude and latitude. Coordinates for any shared sampling sites should be identical, and where sites are raster cells, coordinates are usually expected to be cell centroids.
#' @param r A numeric value specifying the radius (in metres) to use to define radially constrained regions.
#' @param rarefy Whether sites within radially constrained regions should be subsampled (`rarefy = TRUE`, default) or not (`rarefy = FALSE`). If `TRUE`, `nSite` sites will be drawn from `iter` randomly selected radially constrained regions.
#' @param weight Whether sites within the radially constrained region should be drawn at random (`weight = FALSE`, default) or with probability inversely proportional to the square of their distance from the centre of the region (`weight = TRUE`) during subsampling. Only applicable when `rarefy = TRUE`.
#' @param seeding Either `NULL` (the default) or a `data.frame` with additional class `seedMatrix` generated using the [defineSeeds()] or `cookies` (with `output = "seeds"`) functions, specifying the location of the centroids to be used to define radially constrained regions. If `NULL`, each unique populated site within `dat` will be assessed for viability as the center of a radially constrained region.
#' @param iter A numeric value specifying the number of subsamples to be drawn from each radially constrained region if `rarefy = "divvySites"`, `"weightedDivvySites"`, `"sites"`, `"occs"`, or `"sitesThenOccs"`. Default is `100`.
#' @param nSite A numeric value specifying the minimum number of unique populated sites that must fall within a radially constrained region for it to be considered viable. Also, specifies the number of sites that are randomly drawn *without* replacement from each radially constrained region if `rarefy = "sites"` or `"sitesThenOccs"`. Default is `3`.
#' @param crs A character string specifying a coordinate reference system (CRS). You can use the following formats to define coordinate reference systems: WKT, PROJ.4 (e.g., `crs = +proj=longlat +datum=WGS84`), or an EPSG code (e.g., `crs = "EPSG:4326"`). But note that the PROJ.4 notation has been deprecated, and you can only use it with the WGS84/NAD83 and NAD27 datums. Other datums are silently ignored. Default is `"EPSG:8857"`.
#' @param output A character string, either `"full"`,`"locs"`, or `"seeds"`. Specifies whether the returned data should be two columns of subsample site coordinates (`output = "locs"`), the subset of rows from dat associated with those coordinates (`output = "full"`), or a data.frame of the seeds used (`output = "seeds"`).
#'
#' @return A list of length iter.
#' - If `rarefy = TRUE`, each element is a data.frame or matrix (matching the class of `dat`) with *l* rows, where *l* is the number of occurrences associated with `nSite` sites within a viable radially constrained region. The radially constrained region of origin for each subsample is chosen at random.
#' - If `rarefy = FALSE`, each element is a data.frame or matrix (matching the class of `dat`) containing *all* occurrences associated with *all* sites within *each* viable radially constrained region.
#
#' If `output = 'locs'` (default), only the coordinates of the `dat` rows associated with radially constrained regions/subsamples are returned. If `output = 'full'`, all `dat` columns are returned for the rows associated with the radially constrained region/subsample. If `output = 'seeds'`, only the coordinates of the viable seeds are returned.
#'
#' @export
#'
#' @description This function is an expanded version of [divvy::cookies()].
#'
#' In brief, the function takes a single location as a starting (seed) point, circumscribes a circular buffer of `r` metres around it, and then identifies which observations fall within this circle. This process is repeated for each unique site in the dataset (unless these seed points are manually specified). The radially constrained regions (RCRs) that result are then assessed against user-specified criteria; if they contain insufficient data (measured by the number of unique populated sites they contain) and/or overlap with other RCRs by more than is deemed acceptable, they are dropped from the analysis. Those that are retained are then either returned in full (i.e., all occurrences that fall within them are returned) or subsampled.
#'
#' This version of the function differs from the original version in two key ways:
#' 1. Allows users to explore viability of specific regions for spatial subsampling through manual seeding.
#' 2. Allows users to assess how changing values of `r`, `nSite`, and `oThreshold` affects the number of viable RCRs through the new option to return them without subsampling.
#' 3. It allows users to return the viable seed points as a data.frame with additional class `seedMatrix`, which can be used to manually seed (via argument `seeding`) future runs.
#'
#' The function also now returns the coordinates of the seed point in each subsample when `output = "full"`, which can be used to specify the location of radially constrained region from which the subsample is drawn as a random effect in mixed models.
#'
#' # Manual seeding and raw radially constrained region outputs #
#' These functionalities were added to give the function greater exploratory utility. Now users can repeatedly sample the same regions using different radii (`r`) and site minima (`nSite`) to determine optimum values for data sampling.
#'
#' # Subsampling procedures #
#' For explanations of how the subsampling procedures used when `rarefy = TRUE` and `weight = TRUE/FALSE`, see [divvy::cookies()] documentation.
#'
#' # Returning seeds #
#' This functionality was added to further empower users to explore their data.
#'
#' @examples
#' # First, load terra
#' library(terra)
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
#' # Try spatial subsampling
#' # 10km radius, minimum 3 sites, no overlap
#' subsamples <- cookies(pts, xy = xy.name,
#' r = 10000, rarefy = TRUE, output = "full")
cookies <- function(dat, xy, r, rarefy = TRUE, weight = FALSE, seeding = NULL, iter = 100, nSite = 3,
                     crs = "EPSG:8857", output = "locs"){
  coords <- divvy::uniqify(dat, xy)
  coords$id <- paste0("loc", 1:nrow(coords))
  if(is.null(seeding)){
    allPools <- findSeeds(coords, dat, "id", xy, r, nSite, crs)
  } else {
    if(any(class(seeding) == "seedMatrix")){
      allPools <- findSeeds(coords, dat, "id", xy, r, nSite, crs, seeding)
    } else {
      stop("'seeding' is not a 'seedMatrix' object")
    }
  }
  if (length(allPools) < 1) {
    stop("not enough close sites for any subsample or all cookies exceed overlap threshold. Please adjust nSite or oThreshold.")
  }
  seeds <- names(allPools)
  if(output == "seeds"){
    if(!is.null(seeding)){
      return(seeding)
    } else {
      seedMatrix <- coords[which(coords[,"id"] %in% seeds),-which(colnames(coords)=="id")]
      seedMatrix[,"id"] <- paste0("seed", 1:nrow(seedMatrix))
      attr(seedMatrix, "class") <- c("data.frame", "seedMatrix")
      return(seedMatrix)
    }
  } else {
    if(rarefy){
      subsamples <- replicate(iter, cookie(dat, seeds, "id", xy, nSite, allPools, weight, seeding, coords, crs, output, rarefy), simplify = FALSE)
      return(subsamples)
    } else {
      if(!rarefy){
        ## get sites subsamples if rarefying not by divvy
        subsamples <- cookie(dat, seeds, "id", xy, nSite, allPools, weight, seeding, coords, crs, output, rarefy)
        return(subsamples)
      } else {
        stop("`rarefy needs to be TRUE or FALSE")
      }
    }
  }
}
