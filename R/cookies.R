#' Derive radially constrained spatial subsamples from occurrence data
#'
#' @param dat A data.frame or matrix containing taxon names, Cartesian coordinates, and any associated variables.
#' @param xy 	A vector of two elements, specifying the name of columns in `dat` containing Cartesian coordinates, e.g. longitude and latitude. Coordinates for any shared sampling sites should be identical, and where sites are raster cells, coordinates are usually expected to be cell centroids.
#' @param r A numeric value specifying the radius (in metres) to use to define radially constrained regions.
#' @param rarefy Whether sites within radially constrained regions should be subsampled (`rarefy = TRUE`, default) or not (`rarefy = FALSE`). If `TRUE`, `nSite` sites will be drawn from `iter` randomly selected radially constrained regions.
#' @param weight Whether sites within the radially constrained region should be drawn at random (`weight = FALSE`, default) or with probability inversely proportional to the square of their distance from the centre of the region (`weight = TRUE`) during subsampling. Only applicable when `rarefy = TRUE`.
#' @param seeding Either `NULL` (the default) or a matrix generated using [defineSeeds()] function, specifying the location of the centroids to be used to define radially constrained regions. If `NULL`, each unique populated site within `dat` will be assessed for viability as the center of a radially constrained region.
#' @param iter A numeric value specifying the number of subsamples to be drawn from each radially constrained region if `rarefy = "divvySites"`, `"weightedDivvySites"`, `"sites"`, `"occs"`, or `"sitesThenOccs"`. Default is `100`.
#' @param nSite A numeric value specifying the minimum number of unique populated sites that must fall within a radially constrained region for it to be considered viable. Also, specifies the number of sites that are randomly drawn *without* replacement from each radially constrained region if `rarefy = "sites"` or `"sitesThenOccs"`. Default is `3`.
#' @param oThreshold A numeric value between `0` and `1`, specifying the acceptable proportion of overlap in sites or area between viable radially constrained regions. Default is 0 (i.e., all subsamples returned will be spatially independent). Set to 1 to have function behave like [divvy::cookies()].
#' @param oType A character string, either `"sites"` (the default) or `"area"`. If `oThreshold < 1`, `oType` specifies whether the degree of overlap between radially constrained regions should be quantified using sites (`oType = "sites"`) or area (`oType = "area"`).
#' @param oPruningMode A character string, either `"maxOccs"` (the default) or `"minOverlap"`. If `oThreshold < 1`, `oPruningMode` specifies whether the overlapping radially constrained regions with the least associated occurrence data (`oPruningMode = "maxOccs"`) or the most overlap with other regions (`oPruningMode = "minOverlap"`) should be dropped.
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
#' This version of the function differs from the original version in three key ways:
#' 1. Offers users ability to limit or eliminate overlap between spatial subsamples.
#' 2. Allows users to explore viability of specific regions for spatial subsampling through manual seeding.
#' 3. Allows users to assess how changing values of `r`, `nSite`, and `oThreshold` affects the number of viable RCRs through the new option to return them without subsampling.
#'
#' It also allows users to wrap the output of the function in a list which includes the seed points used to derive the sampled RCRs.
#'
#' # Overlap elimination #
#' The ability to control the amount of overlap between RCRs was added so as to give users the ability to eliminate pseudoreplication of spatial regions and standardise the sampling of different geographic areas. The latter is useful when the function is applied to global fossil occurrence datasets, as it curbs the impact of some geographic areas having an order of magnitude more occurrence data associated with them than others (e.g., the oversampling of Western European and North American localities relative to other parts of the world). Currently, the degree of overlap between two RCRs is quantified either as the proportion of sites within an RCR that fall within another (`oType = "sites"`) or the the proportion of the geographic area of an RCR that is included with another (`oType = "area"`). The former is the recommended option and the default.
#'
#' The function offers two methods for iteratively selecting and removing RCRs until those that remain do not exceed the overlap threshold (`oThreshold`): `oPruningMode = "maxOccs"` specifies that the overlapping RCRs with the least associated occurrence data should be the first to be removed, whereas `oPruningMode = "minOverlap"` specifies that the RCRs that overlap the most with other regions should be dropped first. Ties are decided randomly. In practice, the former results in fewer, more richly radially constrained regions being retained, whereas the latter results in a greater number of less richly populated regions. The latter also has a habit of returning RCRs that are centered on the boundaries of densely populated geographic regions.
#'
#' # Manual seeding and raw radially constrained region outputs #
#' These functionalities were added to give the function greater exploratory utility. Now users can repeatedly sample the same regions using different radii (`r`) and site minima (`nSite`) to determine optimum values for data sampling.
#'
#' # Subsampling procedures #
#' For explanations of how the subsampling procedures used when `rarefy = TRUE` and `weight = TRUE/FALSE`, see [divvy::cookies()] documentation.
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
                     oThreshold = 0, oType = "sites", oPruningMode = "maxOccs",
                     crs = "EPSG:8857", output = "locs"){
  coords <- divvy::uniqify(dat, xy)
  coords$id <- paste0("loc", 1:nrow(coords))
  if(is.null(seeding)){
    allPools <- findSeeds(coords, dat, "id", xy, r, nSite, crs, oThreshold, oType, oPruningMode)
  } else {
    allPools <- findSeeds(coords, dat, "id", xy, r, nSite, crs, oThreshold, oType, oPruningMode, seeding)
  }
  if (length(allPools) < 1) {
    stop("not enough close sites for any subsample or all cookies exceed overlap threshold. Please adjust nSite or oThreshold.")
  }
  seeds <- names(allPools)
  if(output == "seeds"){
    if(!is.null(seeding)){
      return(seeding)
    } else {
      out <- coords[which(coords[,"id"] %in% seeds),-which(colnames(coords)=="id")]
      return(out)
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
