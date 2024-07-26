#' Find radially constrained spatial subsamples that fulfil all specified criteria
#'
#' @noRd
findSeeds <- function(coords, dat, siteId, xy, r, nSite, crs = "EPSG:4326", seeding = NULL){
  if(is.null(seeding)){
    sites <- coords[, siteId]
    coordsSV <- terra::vect(coords, geom = xy, crs = crs)
    posPools <- lapply(sites, function(s) {
      seedRow <- which(coords[, siteId] == s)[1]
      sPool <- findPool(seedRow, coordsSV, sites, xy, r, crs)
      n <- length(sPool)
      if (n >= nSite)
        sPool
    })
    names(posPools) <- sites
  } else {
    seedDF <- as.data.frame(matrix("seeds", nrow = nrow(seeding), ncol = ncol(coords)))
    colnames(seedDF) <- colnames(coords)
    seedDF[,which(colnames(coords) %in% xy)] <- seeding[,which(colnames(seeding) %in% xy)]
    seedDF[,which(colnames(coords) %in% siteId)] <- seeding[,which(colnames(seeding) %in% siteId)]
    seedCoords <- rbind(coords, seedDF)
    sites <- seedCoords[,siteId]
    seeds <- seeding[,siteId]
    coordsSV <- terra::vect(seedCoords, geom = xy, crs = crs)
    posPools <- lapply(seeds, function(s) {
      seedRow <- which(seedCoords[, siteId] == s)[1]
      sPool <- findPool(seedRow, coordsSV, sites, xy, r, crs)
      sPool <- sPool[-which(sPool %in% seeds)]
      n <- length(sPool)
      if (n >= nSite)
        sPool
    })
    names(posPools) <- seeds
  }
  posPools <- Filter(Negate(is.null), posPools)
  return(posPools)
}
