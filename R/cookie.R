#' Modified version of divvy internal function "cookie". Final version for export.
#'
#' @noRd
cookie <- function(dat, seeds, siteId, xy, nSite, allPools, weight, seeding, coords, crs, output, rarefy) {
  if(rarefy){
    if(length(seeds) > 1) {
      seed <- sample(sample(seeds), 1)
    } else {
      seed <- seeds
    }
    pool <- allPools[seed][[1]]
    if(weight){
      datSV <- terra::vect(coords, geom = xy, crs = crs)
      pool <- pool[!pool == seed]
      poolBool <- coords[, "id"] %in% pool
      poolPts <- datSV[poolBool, ]
      seedRow <- which(coords[, "id"] == seed)[1]
      seedPt <- datSV[seedRow, ]
      gcdists <- terra::distance(poolPts, seedPt)
      wts <- sapply(gcdists, function(x) x^(-2))
      samplIds <- c(seed, sample(sample(pool), nSite-1, prob = wts, replace = FALSE))
    } else {
      samplIds <- sample(sample(pool), nSite, replace = FALSE)
    }
    coordRows <- match(samplIds, coords$id)
    coordLocs <- coords[coordRows, xy]
    if(output == "full"){
      x <- xy[1]
      y <- xy[2]
      sampPtStrg <- paste(coordLocs[, x], coordLocs[,y], sep = "/")
      datPtStrg <- paste(dat[, x], dat[, y], sep = "/")
      inSamp <- match(datPtStrg, sampPtStrg)
      out <- dat[!is.na(inSamp), ]
      # Need to account for manual seeding
      if(is.null(seeding)){
        seed.xy <- coords[which(coords[,siteId] %in% seed),xy]
        names(seed.xy) <- c(paste0("seed.",xy))
        out <- cbind(out,seed.xy)
      } else {
        seed.xy <- seeding[which(seeding[,"id"] %in% seed),xy]
        names(seed.xy) <- c(paste0("seed.",xy))
        out <- cbind(out,seed.xy)
      }
    } else {
      if (output == "locs") {
        out <- coordLocs
      } else {
        stop("output argument must be one of c('full', 'locs')")
      }
    }
    return(out)
  } else {
    output <- lapply(1:length(seeds), function(s){
      pool <- allPools[seeds[[s]]][[1]]
      coordRows <- match(pool, coords$id)
      coordLocs <- coords[coordRows, xy]
      if (output == "full") {
        x <- xy[1]
        y <- xy[2]
        sampPtStrg <- paste(coordLocs[, x], coordLocs[,y], sep = "/")
        datPtStrg <- paste(dat[, x], dat[, y], sep = "/")
        inSamp <- match(datPtStrg, sampPtStrg)
        out <- dat[!is.na(inSamp), ]
        # Need to account for manual seeding
        if(is.null(seeding)){
          seed.xy <- coords[which(coords[,siteId] %in% seeds[[s]]),xy]
          names(seed.xy) <- c(paste0("seed.",xy))
          out <- cbind(out,seed.xy)
        } else {
          seed.xy <- seeding[which(seeding[,"id"] %in% seeds[[s]]),xy]
          names(seed.xy) <- c(paste0("seed.",xy))
          out <- cbind(out,seed.xy)
        }
      } else {
        if (output == "locs") {
          out <- coordLocs
        } else {
          stop("output argument must be one of c('full', 'locs', 'seeds')")
        }
      }
      return(out)
    })
    return(output)
  }
}
