#' Modified version of divvy internal function "cookie"
#'
#' @param dat D
#' @param seeds D
#' @param xy D
#' @param nSite D
#' @param allPools D
#' @param weight D
#' @param coords D
#' @param crs D
#' @param output D
#' @param divvyRarefaction D
#' @param returnSeeds D
#'
#' @noRd
#'
#' @return D
#'
#' @examples
#' print("example")
cookie2 <- function(dat, seeds, xy, nSite, allPools, weight, coords, crs, output, divvyRarefaction, returnSeeds) {
  if(divvyRarefaction){
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
    } else {
      if (output == "locs") {
        out <- coordLocs
      } else {
        stop("output argument must be one of c('full', 'locs')")
      }
    }
    if(returnSeeds){
      return(list("seed" = seed, "subsample" = out))
    } else {
      return(out)
    }
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
      } else {
        if (output == "locs") {
          out <- coordLocs
        } else {
          stop("output argument must be one of c('full', 'locs')")
        }
      }
      return(out)
    })
    return(output)
  }
}
