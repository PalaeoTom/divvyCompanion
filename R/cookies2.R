#' Derive radially constrained spatial subsamples from occurrence data
#'
#' @param dataMat A data.frame or matrix containing taxon names, Cartesian coordinates, and any associated variables.
#' @param xy 	A vector of two elements, specifying the name or numeric position of columns in `dataMat` containing Cartesian coordinates, e.g. longitude and latitude. Coordinates for any shared sampling sites should be identical, and where sites are raster cells, coordinates are usually expected to be cell centroids.
#' @param r A numeric value specifying the radius (in metres) to use to define radially constrained regions.
#' @param seeding Either `NULL` (the default) or a matrix generated using [rasterOccData()] function, specifying the location of the centroids to be used to define radially constrained regions. If `NULL`, each unique populated site within `dataMat` will be assessed for viability as the center of a radially constrained region.
#' @param rarefaction A character string specifying the subsampling method to be used, if any. Options are `"divvySites"` (the default, equivalent to [divvy::cookies()] with `weight = FALSE`), `"weightedDivvySites"` (equivalent to [divvy::cookies()] with `weight = TRUE`), `"sites"`, `"occs"`, `"sitesThenOccs"`, and `"none"`. See below for further details.
#' @param reps A numeric value specifying the number of subsamples to be drawn from each radially constrained region if `rarefaction = "divvySites"`, `"weightedDivvySites"`, `"sites"`, `"occs"`, or `"sitesThenOccs"`. Default is `100`.
#' @param nSites A numeric value specifying the minimum number of unique populated sites that must fall within a radially constrained region for it to be considered viable. Also, specifies the number of sites that are randomly drawn with replacement from each radially constrained region if `rarefaction = "sites"` or `"sitesThenOccs"`. Default is `3`.
#' @param nOccs A numeric value specifying the number of occurrences to be randomly drawn with replacement from each radially constrained region if `rarefaction = "occs"` or `"sitesThenOccs"`. Default is `100`.
#' @param oThreshold A numeric value between `0` and `1`, specifying the acceptable proportion of overlap in sites or area between viable radially constrained regions. Default is 0 (i.e., all subsamples returned will be spatially independent). Set to 1 to have function behave like [divvy::cookies()].
#' @param oType A character string, either `"sites"` (the default) or `"area"`. If `oThreshold < 1`, `oType` specifies whether the degree of overlap between radially constrained regions should be quantified using sites (`oType = "sites"`) or area (`oType = "area"`).
#' @param oPruningMode A character string, either `"maxOccs"` (the default) or `"minOverlap"`. If `oThreshold < 1`, `oPruningMode` specifies whether the overlapping radially constrained regions with the least associated occurrence data (`oPruningMode = "maxOccs"`) or the most overlap with other regions (`oPruningMode = "minOverlap"`) should be dropped. In practice, the former results in fewer, more richly radially constrained regions being retained, whereas the latter results in a greater number of less richly populated regions.
#' @param returnSeeds `TRUE` or `FALSE` (the default). If `TRUE`, the output of this function is nested within a list that includes a second element containing the coordinates of the centroids of each radially constrained region represented in the other element of the output.
#' @param crs A character string specifying a coordinate reference system (CRS). You can use the following formats to define coordinate reference systems: WKT, PROJ.4 (e.g., `crs = +proj=longlat +datum=WGS84`), or an EPSG code (e.g., `crs = "EPSG:4326"`). But note that the PROJ.4 notation has been deprecated, and you can only use it with the WGS84/NAD83 and NAD27 datums. Other datums are silently ignored. Default is `"EPSG:4326"`.
#' @param output A character string, either `"full"` or `"locs"`. Specifies whether the returned data should be two columns of subsample site coordinates (`output = "locs"`) or the subset of rows from dataMat associated with those coordinates (`output = "full"`).
#'
#' @return Output format changes depending on arguments used.
#' If `returnSeeds = FALSE` and:
#' - `rarefaction = "none"`: a list of length *n*, where *n* is the number of radially constrained regions identified as viable. Each element is a data.frame or matrix (matching the class of `dataMat`) containing *k* rows, where *k* is the number of occurrences within dataMat that fall within a specific radially constrained region.
#' - `rarefaction = "divvySites"` or `"weightedDivvySites"`: a list of length reps. Each element is a data.frame or matrix (matching the class of `dataMat`) with *l* rows, where *l* is the number of occurrences associated with `nSites` sites within a viable radially constrained region. The radially constrained region of origin for each subsample is chosen at random. If `rarefaction = "weightedDivvySites"`, the first row in each returned subsample data.frame corresponds to the seed point (i.e., the centroid of the radially constrained region). If `rarefaction = "divvySites"`, rows are listed in the random order of which they were drawn.
#' - `rarefaction = "sites"`, `"occs"`, or `"sitesThenOccs"`: a list of length *n*, where *n* is the number of radially constrained regions identified as viable. Each element is a list of length `reps`, each element of which is a data.frame or matrix (matching the class of dataMat). The dimensions of these data.frames/matrices depend on the rarefaction method used:
#'    + if `rarefaction = "sites"`, each data.frame/matrix will have *l* rows, where *l* is the number of occurrences associated with `nSites` sites within a viable radially constrained region.
#'    + if `rarefaction = "occs"`, each data.frame/matrix will have `nOccs` rows, where `nOccs` is a random sampling of the occurrences associated with a viable radially constrained region.
#'    + if `rarefaction = "sitesThenOccs"`, each data.frame/matrix will have `nOccs` rows, where `nOccs` is a random sampling of the occurrences associated with `nSites` randomly selected sites within a viable radially constrained region.
#'
#' If `returnSeeds = TRUE`: a list with two elements, where the second is main output of the function specified using argument `rarefaction`. The first element is a data.frame or matrix (matching the class of `dataMat`) containing the rows used as seed points for the radially constrained regions.
#'
#' If `output = 'locs'` (default), only the coordinates of the `dataMat` rows associated with radially constrained regions/subsamples are returned. If `output = 'full'`, all `dataMat` columns are returned for the rows associated with the radially constrained region/subsample.
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
#' 3. Allows users to assess how different values for `r`, `nSites`, and `oThreshold` affect number of viable RCRs through option to return them without subsampling.
#' 4. Implements new subsampling procedures inspired by classical occurrence rarefaction.
#'
#' This version of the function also allows users to wrap the output of the function in a list which includes the seed points used to derive the sampled RCRs.
#'
#'
#'
#' For explanations of how the subsampling procedures used when `rarefaction = "divvySites"` or `"weightedDivvySites"`, see [divvy::cookies()] documentation.
#'
#' The three new subsampling procedures are described below. Each of them differ from the two subsampling procedures offered by [divvy::cookies()] in that they draw the `reps` subsamples from each viable RCR, rather than randomly drawing `reps` subsamples from `reps` randomly selected viable RCRs. When used in conjunction with a low threshold for overlap (i.e., `oThreshold = 0`), this limits oversampling of occurrence-rich geographic areas.
#' 1. `rarefaction = "sites"`: draws `nSites` sites randomly with replacement from each RCR. Each subsample contains all occurrences associated with the RCR.
#' 2. `rarefaction = "occs"`: draws `nOccs` occurrences randomly with replacement from each RCR. Each subsample contains `nOccs` occurrences drawn from all sites within the RCR.
#' 3. `rarefaction = "sitesThenOccs"`: draws `nSites` sites randomly with replacement from each RCR, then draws `nOccs` occurrences from the occurrences associated with selected sites. Each subsample contains `nOccs` occurrences drawn from `nSites` sites within the RCR.
#'
#' @examples
#' print("example")
cookies2 <- function(dataMat, xy, r, seeding = NULL, rarefaction = "divvySites", reps = 100, nSites = 3, nOccs = 100,
                     oThreshold = 0, oType = "sites", oPruningMode = "maxOccs",
                     returnSeeds = F, crs = "EPSG:4326", output = "locs"){
  ## If rarefaction = "weightedDivvy", set argument
  if(rarefaction == "weightedDivvySites"){
    weight <- T
  } else {
    weight <- F
  }
  coords <- divvy::uniqify(dataMat, xy)
  coords$id <- paste0("loc", 1:nrow(coords))
  if(is.null(seeding)){
    allPools <- findSeeds2(coords, dataMat, "id", xy, r, nSites, crs, oThreshold, oType, oPruningMode)
  } else {
    allPools <- findSeeds2(coords, dataMat, "id", xy, r, nSites, crs, oThreshold, oType, oPruningMode, seeding)
  }
  if (length(allPools) < 1) {
    stop("not enough close sites for any subsample or all cookies exceed overlap oThreshold. Please adjust nSites or oThreshold.")
  }
  seeds <- names(allPools)
  if(rarefaction == "divvySites" || rarefaction == "weightedDivvySites"){
    subsamples <- replicate(reps, cookie2(dataMat, seeds, xy, nSites, allPools, weight, coords, crs, output, divvyRarefaction = T, returnSeeds), simplify = FALSE)
    if(returnSeeds){
      usedSeeds <- sapply(1:length(subsamples), function(x) subsamples[[x]][[1]])
      seed_out <- coords[which(coords$id %in% usedSeeds),]
      subsamples_out <- lapply(1:length(subsamples), function(x) subsamples[[x]][[2]])
      names(subsamples_out) <- usedSeeds
      return(list("seeds" = seed_out, "subsamples" = subsamples_out))
    } else {
      return(subsamples)
    }
  } else {
    ## get subsamples
    subsamples <- cookie2(dataMat, seeds, xy, nSites, allPools, weight, coords, crs, output, divvyRarefaction = F, returnSeeds)
    if(rarefaction == "sites" || rarefaction == "occs" || rarefaction == "sitesThenOccs"){
      ## rarefy by sites
      if(rarefaction == "sites"){
        rareSubs <- lapply(1:length(subsamples), function(x){
          ## Get list of unique cells
          cells <- unique(subsamples[[x]][,"cell"])
          ## Get reps samples
          sub2samples <- lapply(1:reps, function(all){
            ## Get sample of cells
            samp <- sample(cells, size = nSites, replace = T)
            ## extract data frame
            s2sample <- subsamples[[x]][which(subsamples[[x]][,"cell"] %in% samp),]
          })
        })
        ## prepare output
        if(returnSeeds){
          names(rareSubs) <- seeds
          seed_out <- coords[which(coords$id %in% seeds),]
          return(list("seeds" = seed_out, "subsamples" = subsamples))
        } else {
          return(rareSubs)
        }
      }
      ## rarefy by occurrences
      if(rarefaction == "occs"){
        rareSubs <- lapply(1:length(subsamples), function(x){
          ## Get repOccs samples
          sub2samples <- lapply(1:reps, function(all){
            ## Get sample of occurrences
            samp <- sample(1:nrow(subsamples[[x]]), size = nOccs, replace = T)
            ## extract data frame
            s2sample <- subsamples[[x]][samp,]
          })
        })
        ## prepare output
        if(returnSeeds){
          names(rareSubs) <- seeds
          seed_out <- coords[which(coords$id %in% seeds),]
          return(list("seeds" = seed_out, "subsamples" = subsamples))
        } else {
          return(rareSubs)
        }
      }
      ## rarefy by sites, then occurrences
      if(rarefaction == "sitesThenOccs"){
        rareSubs <- lapply(1:length(subsamples), function(x){
          ## Get list of unique cells
          cells <- unique(subsamples[[x]][,"cell"])
          ## Get repOccs samples
          sub2samples <- lapply(1:reps, function(all){
            ## Get sample of cells
            sampCells <- sample(cells, size = nSites, replace = T)
            ## Get occurrences with these cells
            sampCellOccs <- which(subsamples[[x]][,"cell"] %in% sampCells)
            ## Get sample of occurrences
            sampOccs <- sample(sampCellOccs, size = nOccs, replace = T)
            ## extract data frame
            s2sample <- subsamples[[x]][sampOccs,]
          })
        })
        ## prepare output
        if(returnSeeds){
          names(rareSubs) <- seeds
          seed_out <- coords[which(coords$id %in% seeds),]
          return(list("seeds" = seed_out, "subsamples" = subsamples))
        } else {
          return(rareSubs)
        }
      }
    } else {
      if(rarefaction == "none"){
        if(returnSeeds){
          names(subsamples) <- seeds
          seed_out <- coords[which(coords$id %in% seeds),]
          return(list("seeds" = seed_out, "subsamples" = subsamples))
        } else {
          return(subsamples)
        }
      } else {
        stop("argument rarefaction is not 'divvySites', 'weightedDivvySites', 'sites', 'occs', 'sitesThenOccs', or 'none'")
      }
    }
  }
}
