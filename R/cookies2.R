#' Derive radially constrained spatial subsamples from occurrence data
#'
#' @param dat A data.frame or matrix containing taxon names, Cartesian coordinates, and any associated variables.
#' @param xy 	A vector of two elements, specifying the name of columns in `dat` containing Cartesian coordinates, e.g. longitude and latitude. Coordinates for any shared sampling sites should be identical, and where sites are raster cells, coordinates are usually expected to be cell centroids.
#' @param uniqID A character string specifying the name of the column in `dat` containing unique site identifiers (e.g., cell numbers added via [rasterOccData()]). Default is `"cell"`.
#' @param r A numeric value specifying the radius (in metres) to use to define radially constrained regions.
#' @param seeding Either `NULL` (the default) or a matrix generated using [defineSeeds()] function, specifying the location of the centroids to be used to define radially constrained regions. If `NULL`, each unique populated site within `dat` will be assessed for viability as the center of a radially constrained region.
#' @param rarefaction A character string specifying the subsampling method to be used, if any. Options are `"divvySites"` (the default, equivalent to [divvy::cookies()] with `weight = FALSE`), `"weightedDivvySites"` (equivalent to [divvy::cookies()] with `weight = TRUE`), `"sites"`, `"occs"`, `"sitesThenOccs"`, and `"none"`. See below for further details.
#' @param iter A numeric value specifying the number of subsamples to be drawn from each radially constrained region if `rarefaction = "divvySites"`, `"weightedDivvySites"`, `"sites"`, `"occs"`, or `"sitesThenOccs"`. Default is `100`.
#' @param nSite A numeric value specifying the minimum number of unique populated sites that must fall within a radially constrained region for it to be considered viable. Also, specifies the number of sites that are randomly drawn *without* replacement from each radially constrained region if `rarefaction = "sites"` or `"sitesThenOccs"`. Default is `3`.
#' @param nOcc A numeric value specifying the number of occurrences to be randomly drawn *with* replacement from each radially constrained region if `rarefaction = "occs"` or `"sitesThenOccs"`. Default is `100`.
#' @param oThreshold A numeric value between `0` and `1`, specifying the acceptable proportion of overlap in sites or area between viable radially constrained regions. Default is 0 (i.e., all subsamples returned will be spatially independent). Set to 1 to have function behave like [divvy::cookies()].
#' @param oType A character string, either `"sites"` (the default) or `"area"`. If `oThreshold < 1`, `oType` specifies whether the degree of overlap between radially constrained regions should be quantified using sites (`oType = "sites"`) or area (`oType = "area"`).
#' @param oPruningMode A character string, either `"maxOccs"` (the default) or `"minOverlap"`. If `oThreshold < 1`, `oPruningMode` specifies whether the overlapping radially constrained regions with the least associated occurrence data (`oPruningMode = "maxOccs"`) or the most overlap with other regions (`oPruningMode = "minOverlap"`) should be dropped.
#' @param crs A character string specifying a coordinate reference system (CRS). You can use the following formats to define coordinate reference systems: WKT, PROJ.4 (e.g., `crs = +proj=longlat +datum=WGS84`), or an EPSG code (e.g., `crs = "EPSG:4326"`). But note that the PROJ.4 notation has been deprecated, and you can only use it with the WGS84/NAD83 and NAD27 datums. Other datums are silently ignored. Default is `"EPSG:8857"`.
#' @param output A character string, either `"full"`,`"locs"`, or `"seeds"`. Specifies whether the returned data should be two columns of subsample site coordinates (`output = "locs"`), the subset of rows from dat associated with those coordinates (`output = "full"`), or a data.frame of the seeds used (`output = "seeds"`).
#'
#' @return Output format changes depending on arguments used.
#' - `rarefaction = "none"`: a list of length *n*, where *n* is the number of radially constrained regions identified as viable. Each element is a data.frame or matrix (matching the class of `dat`) containing *k* rows, where *k* is the number of occurrences within dat that fall within a specific radially constrained region.
#' - `rarefaction = "divvySites"` or `"weightedDivvySites"`: a list of length iter. Each element is a data.frame or matrix (matching the class of `dat`) with *l* rows, where *l* is the number of occurrences associated with `nSite` sites within a viable radially constrained region. The radially constrained region of origin for each subsample is chosen at random. If `rarefaction = "weightedDivvySites"`, the first row in each returned subsample data.frame corresponds to the seed point (i.e., the centroid of the radially constrained region). If `rarefaction = "divvySites"`, rows are listed in the random order of which they were drawn.
#' - `rarefaction = "sites"`, `"occs"`, or `"sitesThenOccs"`: a list of length *n*, where *n* is the number of radially constrained regions identified as viable. Each element is a list of length `iter`, each element of which is a data.frame or matrix (matching the class of dat). The dimensions of these data.frames/matrices depend on the rarefaction method used:
#'    + if `rarefaction = "sites"`, each data.frame/matrix will have *l* rows, where *l* is the number of occurrences associated with `nSite` sites within a viable radially constrained region.
#'    + if `rarefaction = "occs"`, each data.frame/matrix will have `nOcc` rows, where `nOcc` is a random sampling of the occurrences associated with a viable radially constrained region.
#'    + if `rarefaction = "sitesThenOccs"`, each data.frame/matrix will have `nOcc` rows, where `nOcc` is a random sampling of the occurrences associated with `nSite` randomly selected sites within a viable radially constrained region.
#'
#' If `output = 'locs'` (default), only the coordinates of the `dat` rows associated with radially constrained regions/subsamples are returned. If `output = 'full'`, all `dat` columns are returned for the rows associated with the radially constrained region/subsample. If `output = 'seeds'`, only the coordinates of the viable seeds are returned.
#'
#' @export
#'
#' @description This function is an expanded version of [divvy::cookies()].
#'
#' In brief, the function takes a single location as a starting (seed) point, circumscribes a circular buffer of `r` metres around it, and then identifies which observations fall within this circle. This process is repeated for each unique site in the dataset (unless these seed points are manually specified). The radially constrained regions (RCRs) that result are then assessed against user-specified criteria; if they contain insufficient data (measured by the number of unique populated sites they contain) and/or overlap with other RCRs by more than is deemed acceptable, they are dropped from the analysis. Those that are retained are then either returned in full (i.e., all occurrences that fall within them are returned) or subsampled.
#'
#' This version of the function differs from the original version in four key ways:
#' 1. Offers users ability to limit or eliminate overlap between spatial subsamples.
#' 2. Allows users to explore viability of specific regions for spatial subsampling through manual seeding.
#' 3. Allows users to assess how changing values of `r`, `nSite`, and `oThreshold` affects the number of viable RCRs through the new option to return them without subsampling.
#' 4. Implements new subsampling procedures inspired by classical occurrence rarefaction.
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
#' For explanations of how the subsampling procedures used when `rarefaction = "divvySites"` or `"weightedDivvySites"`, see [divvy::cookies()] documentation.
#'
#' The three new subsampling procedures are described below. Each of them differ from the two subsampling procedures offered by [divvy::cookies()] in that they draw the `iter` subsamples from each viable RCR, rather than randomly drawing `iter` subsamples from `iter` randomly selected viable RCRs. When used in conjunction with a low threshold for overlap (i.e., `oThreshold = 0`), this limits oversampling of occurrence-rich geographic areas.
#' 1. `rarefaction = "sites"`: draws `nSite` sites randomly *without* replacement from each RCR. Each subsample contains all occurrences associated with the RCR.
#' 2. `rarefaction = "occs"`: draws `nOcc` occurrences randomly *with* replacement from each RCR. Each subsample contains `nOcc` occurrences drawn from all sites within the RCR.
#' 3. `rarefaction = "sitesThenOccs"`: draws `nSite` sites randomly *without* replacement from each RCR, then randomly draws `nOcc` occurrences *with* replacement from the occurrences associated with selected sites. Each subsample contains `nOcc` occurrences drawn from `nSite` sites within the RCR.
#'
#' @examples
#' # Two examples, first with non-rasterised data
#' # Generate occurrence data
#' n <- 100
#' set.seed(5)
#'
#' # 100 sets of x and y coordinates
#' x <- runif(n, 0, 50)
#' y <- runif(n, 0, 50)
#'
#' # Give each site a unique identifier
#' z <- seq(1, n, 1)
#'
#' # Combine into data frame and label columns
#' pts <- data.frame(x, y, z)
#' colnames(pts) <- c("x", "y", "z")
#'
#' # Derive subsamples, 10km radius regions
#' # Minimum 3 sites, no overlap
#' subsamples <- cookies2(dat = pts, uniqID = "z",
#' xy = c("x", "y"), iter = 5, nSite = 3, r = 10000)
cookies2 <- function(dat, xy, uniqID = "cell", r, seeding = NULL, rarefaction = "divvySites", iter = 100, nSite = 3, nOcc = 100,
                     oThreshold = 0, oType = "sites", oPruningMode = "maxOccs",
                     crs = "EPSG:8857", output = "locs"){
  ## If rarefaction = "weightedDivvy", set argument
  if(rarefaction == "weightedDivvySites"){
    weight <- T
  } else {
    weight <- F
  }
  coords <- divvy::uniqify(dat, xy)
  coords$id <- paste0("loc", 1:nrow(coords))
  if(is.null(seeding)){
    allPools <- findSeeds2(coords, dat, "id", xy, uniqID, r, nSite, crs, oThreshold, oType, oPruningMode)
  } else {
    allPools <- findSeeds2(coords, dat, "id", xy, uniqID, r, nSite, crs, oThreshold, oType, oPruningMode, seeding)
  }
  if (length(allPools) < 1) {
    stop("not enough close sites for any subsample or all cookies exceed overlap threshold. Please adjust nSite or oThreshold.")
  }
  seeds <- names(allPools)
  if(output == "seeds"){
    if(seeding){
      return(seeding)
    } else {
      out <- coords[which(coords[,"id"] %in% seeds),-which(colnames(coords)=="id")]
      return(out)
    }
  } else {
    if(rarefaction == "divvySites" || rarefaction == "weightedDivvySites"){
      subsamples <- replicate(iter, cookie2(dat, seeds, xy, nSite, allPools, weight, coords, crs, output, divvyRarefaction = T), simplify = FALSE)
      return(subsamples)
    } else {
      ## get sites subsamples if rarefying not by divvy
      subsamples <- cookie2(dat, seeds, xy, nSite, allPools, weight, coords, crs, output, divvyRarefaction = F)
      if(rarefaction == "sites" || rarefaction == "occs" || rarefaction == "sitesThenOccs"){
        ## rarefy by sites
        if(rarefaction == "sites"){
          rareSubs <- lapply(1:length(subsamples), function(x){
            ## Get list of unique cells
            cells <- unique(subsamples[[x]][,uniqID])
            ## Get iter samples
            sub2samples <- lapply(1:iter, function(all){
              ## Get sample of cells
              samp <- sample(cells, size = nSite, replace = F)
              ## extract data frame
              s2sample <- subsamples[[x]][which(subsamples[[x]][,uniqID] %in% samp),]
            })
          })
          ## prepare output
          return(rareSubs)
        }
        ## rarefy by occurrences
        if(rarefaction == "occs"){
          rareSubs <- lapply(1:length(subsamples), function(x){
            ## Get repOccs samples
            sub2samples <- lapply(1:iter, function(all){
              ## Get sample of occurrences
              samp <- sample(1:nrow(subsamples[[x]]), size = nOcc, replace = T)
              ## extract data frame
              s2sample <- subsamples[[x]][samp,]
            })
          })
          ## prepare output
          return(rareSubs)
        }
        ## rarefy by sites, then occurrences
        if(rarefaction == "sitesThenOccss"){
          rareSubs <- lapply(1:length(subsamples), function(x){
            ## Get list of unique cells
            cells <- unique(subsamples[[x]][,uniqID])
            ## Get repOccs samples
            sub2samples <- lapply(1:iter, function(all){
              ## Get sample of cells
              sampCells <- sample(cells, size = nSite, replace = F)
              ## Get occurrences with these cells
              sampCellOccs <- which(subsamples[[x]][,uniqID] %in% sampCells)
              ## Get sample of occurrences
              sampOccs <- sample(sampCellOccs, size = nOcc, replace = T)
              ## extract data frame
              s2sample <- subsamples[[x]][sampOccs,]
            })
          })
          return(rareSubs)
        }
      } else {
        if(rarefaction == "none"){
          return(subsamples)
        } else {
          stop("argument rarefaction is not 'divvySites', 'weightedDivvySites', 'sites', 'occs', 'sitesThenOccss', or 'none'")
        }
      }
    }
  }
}
