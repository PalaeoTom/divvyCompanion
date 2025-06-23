#' Derive radially constrained spatial subsamples from occurrence data
#'
#' @param dat A data.frame or matrix containing taxon names, Cartesian coordinates, and any associated variables.
#' @param xy 	A vector of two elements, specifying the name of columns in `dat` containing Cartesian coordinates, e.g. longitude and latitude. Coordinates for any shared sampling sites should be identical, and where sites are raster cells, coordinates are usually expected to be cell centroids.
#' @param uniqID A character string specifying the name of the column in `dat` containing unique site identifiers (e.g., cell numbers added via [rasterOccData()]). Default is `"cell"`.
#' @param r A numeric value specifying the radius (in metres) to use to define radially constrained regions.
#' @param standardiseCells Either `TRUE` (the default) or `FALSE`. If `TRUE`, `nCookie` samples containing `nSite` grid cells are drawn from each radially constrained region and returned. If `FALSE`, the raw radially constrained regions are returned.
#' @param exhaustClusters Either `TRUE` or `FALSE` (the default). See below for explanation.
#' @param nSite A numeric value specifying the minimum number of unique populated sites that must fall within a radially constrained region for it to be considered viable. Also, specifies the number of sites that are randomly drawn *without* replacement from each radially constrained region if `standardiseCells = TRUE`. Default is `3`.
#' @param nCookie A numeric value specifying the number of subsamples of grid cells to be drawn from each radially constrained region if `standardiseCells = TRUE`. Default is `100`.
#' @param nBatch A numeric value specifying the number of subsamples of radially constrained regions to be drawn from the clusters identified. Default is `100`.
#' @param crs A character string specifying a coordinate reference system (CRS). You can use the following formats to define coordinate reference systems: WKT, PROJ.4 (e.g., `crs = +proj=longlat +datum=WGS84`), or an EPSG code (e.g., `crs = "EPSG:4326"`). But note that the PROJ.4 notation has been deprecated, and you can only use it with the WGS84/NAD83 and NAD27 datums. Other datums are silently ignored. Default is `"EPSG:8857"`.
#' @param output A character string, either `"full"` or `"seeds"`. Specifies whether the returned data should be the subset of rows from dat associated with those coordinates (`output = "full"`) or a data.frame of the seeds used (`output = "seeds"`).
#' @param n.cores A numeric value specifying the number of cores to used. Note, this parallelisation method relies on forking, and so `n.cores` must equal `1` on Windows.
#'
#' @return Output format changes depending on arguments used.
#' - `standardiseCells = FALSE`: a list of length *n*, where *n* is the number of clusters identified as viable. Each element of this list is a list of length *c*, where *c* is the number of radially constrained regions identified as viable within each cluster. Each element of these nested lists is a data.frame or matrix (matching the class of `dat`) containing *k* rows, where *k* is the number of occurrences within dat that fall within a specific radially constrained region.
#' - `standardiseCells = TRUE`: a list of length *n*, where *n* is the number of radially constrained regions sampled. If `exhaustClusters = TRUE`, then clusters may be represented by more than one radially constrained regions. Otherwise, each cluster will be represented by a single radially constrained region. Each element of this list is a list of length `nCookie`, each element of which is a data.frame or matrix (matching the class of dat) containing all occurrences associated with the `nSite` grid cells sampled *without replacement* from the radially constrained region.
#
#' If `output = 'full'`, all `dat` columns are returned for the rows associated with the radially constrained region/subsample. If `output = 'seeds'`, only the coordinates of the viable seeds are returned.
#'
#' If no viable clusters are present in the data, output is `NULL`.
#'
#' @export
#'
#' @description This function is a modified version of [divvy::cookies()].
#'
#' In brief, the function takes a single location as a starting (seed) point, circumscribes a circular buffer of `r` metres around it, and then identifies which observations fall within this circle. This process is repeated for each unique site in the dataset (unless these seed points are manually specified). The radially constrained regions (RCRs) that result are then assessed against user-specified criteria; if they contain insufficient data (measured by the number of unique populated sites they contain), they are dropped from the analysis. These RCRs are then grouped into spatially distinct clusters; RCRs within clusters may overlap by a cell or more but each cluster is distinct. The RCRs are then either returned in full (i.e., all occurrences that fall within them are returned) or subsampled.
#'
#' **Differences between this function and divvy's cookies:**
#' 1. Grouping of overlapping RCRs into clusters which are themselves subsampled
#' This extra layer of iteration was added so as to give users the ability to eliminate pseudoreplication of spatial regions, and therefore standardise the sampling of different geographic areas, without dropping heaps of occurrence of data from the analysis. The latter is useful when the function is applied to global fossil occurrence datasets, as it curbs the impact of some geographic areas having an order of magnitude more occurrence data associated with them than others (e.g., the oversampling of Western European and North American localities relative to other parts of the world). Overlap is defined as the sharing of 1 or more grid cells between radially constrained regions. Note, some radially constrained regions within clusters do not overlap but are bridged by one or more RCRs that collectively share at least 1 grid cell with each non-overlapping RCR. If `standardiseCells = TRUE` and `exhaustClusters = TRUE`, multiple non-overlapping RCRs may be selected to represent a cluster within a single subsample.
#'
#' 2. Subsampling of `nSite` grid cells of each radially constrained region `nCookie` times.
#' This functions takes a different approach to subsampling radially constrained regions than [divvy::cookies()]. It differs in that it draws the `nCookie` subsamples from each viable RCR, rather than randomly drawing `nCookie` subsamples from `nCookie` randomly selected viable RCRs. In conjunction with the clustering approach implemented, this guarantees that no grid cell will be sampled twice in a subsample, which is a possibility when [divvy::cookies()] is applied to occurrence data with high density clusters of points.
#'
#' @examples
#' # Two examples, first with non-rasterised data
#' # Generate occurrence data
#' n <- 100
#' set.seed(5)
#'
#' # 100 sets of x and y coordinates
#' x <- runif(n, 0, 100)
#' y <- runif(n, 0, 100)
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
#' subsamples <- traybake(dat = pts, xy = c("x", "y"),
#' uniqID = "z", r = 25, standardiseCells = TRUE,
#' exhaustClusters = TRUE, nSite = 3,
#' nCookie = 100, nBatch = 100)
traybake <- function(dat, xy, uniqID = "cell", r, standardiseCells = T, exhaustClusters = F, nSite = 3, nCookie = 100, nBatch = 100, crs = "EPSG:8857", output = "full", n.cores = 1){
  ## Reduce to unique cells
  coords <- divvy::uniqify(dat, xy)
  ## Label unique cells
  coords$id <- paste0("loc", 1:nrow(coords))
  ## Find clusters
  clusters <- findClusters(coords, "id", xy, r, nSite, crs)
  ## If no clusters, return NULL
  if (length(clusters) < 1) {
    return(NULL)
  }
  ## Get clusters and seeds
  if(output == "seeds"){
    seeds <- lapply(1:length(clusters), function(x){
      coords[which(coords[,"id"] %in% names(clusters[[x]])), -which(colnames(coords)=="id")]
    })
    names(seeds) <- names(clusters)
    return(seeds)
  } else {
    ## Otherwise, need to get batches
    if(exhaustClusters){
      batches <- lapply(1:nBatch, function(all){
        batch <- lapply(1:length(clusters), function(y){
          ## Exhaust cluster
          cluster <- clusters[[y]]
          ## Define tracker
          tracker <- names(cluster)
          ## Create out
          out <- list()
          ## Draw 1
          while(TRUE){
            ## Sample
            s <- sample(1:length(cluster), size = 1)
            ## get seed
            seed <- names(cluster)[s]
            ## get RCR
            rcr <- cluster[[s]]
            ## Add to out
            out <- c(out,list(c(rcr,seed)))
            ## Drop all RCRs that overlap with RCR
            ## Get boolean
            bool <- sapply(1:length(cluster), function(z) all(!rcr %in% cluster[[z]]))
            ## Update cluster and tracker
            tracker <- tracker[bool]
            cluster <- cluster[bool]
            if(length(tracker)==0){
              break
            }
          }
          return(out)
        })
        batch <- unlist(batch, recursive = F)
        ## After batch, isolate seeds
        seeds <- sapply(1:length(batch), function(x) batch[[x]][length(batch[[x]])])
        ## Remove them
        batch <- lapply(1:length(batch), function(x) batch[[x]][-length(batch[[x]])])
        ## name batch
        names(batch) <- seeds
        return(batch)
      })
      ## Get single RCR per cluster
    } else {
      batches <- lapply(1:nBatch, function(all){
        batch <- lapply(1:length(clusters), function(y){
          s <- sample(1:length(clusters[[y]]), size = 1)
          seed <- names(clusters[[y]])[s]
          out <- c(clusters[[y]][[s]],seed)
        })
        ## After batch, isolate seeds
        seeds <- sapply(1:length(batch), function(x) batch[[x]][length(batch[[x]])])
        ## Remove them
        batch <- lapply(1:length(batch), function(x) batch[[x]][-length(batch[[x]])])
        ## name batch
        names(batch) <- seeds
        return(batch)
      })
    }
    ## Once batches are defined, use cookie2 to get occurrences associated, Do this for each
    raw <- lapply(1:length(batches), function(x){
      subsamples <- cookie2(dat, seeds = names(batches[[x]]), xy, nSite, batches[[x]], weight = F, coords, crs, output, divvyRarefaction = F)
      names(subsamples) <- names(batches[[x]])
      return(subsamples)
    })
    ## if standardiseCells = TRUE, rarefy
    if(standardiseCells){
      subsamples <- parallel::mclapply(1:length(raw), mc.cores = n.cores, function(y){
        rareSubs <- lapply(1:length(raw[[y]]), function(x){
        ## Get list of unique cells
        cells <- unique(raw[[y]][[x]][,uniqID])
        ## Get nCookie samples
        sub2samples <- lapply(1:nCookie, function(all){
          ## Get sample of cells
          samp <- sample(cells, size = nSite, replace = F)
          ## extract data frame
          s2sample <- raw[[y]][[x]][which(raw[[y]][[x]][,uniqID] %in% samp),]
          })
        })
        names(rareSubs) <- names(raw[[y]])
        return(rareSubs)
        })
    ## prepare output
    return(subsamples)
    } else {
    return(raw)
    }
  }
}
