#' Find radially constrained spatial subsamples that fulfil all specified criteria
#'
#' @noRd
findSeeds2 <- function(dat, rawData, siteId, xy, uniqID, r, nSite, crs = "EPSG:4326", oThreshold = 0, oType = "sites", oPruningMode = "maxOccs", seeding = NULL){
  if(is.null(seeding)){
    sites <- dat[, siteId]
    datSV <- terra::vect(dat, geom = xy, crs = crs)
    posPools <- lapply(sites, function(s) {
      seedRow <- which(dat[, siteId] == s)[1]
      sPool <- findPool2(seedRow, datSV, sites, xy, r, crs)
      n <- length(sPool)
      if (n >= nSite)
        sPool
    })
    names(posPools) <- sites
  } else {
    seedDF <- as.data.frame(matrix("seeds", nrow = nrow(seeding), ncol = ncol(dat)))
    colnames(seedDF) <- colnames(dat)
    seedDF[,which(colnames(dat) %in% xy)] <- seeding[,which(colnames(seeding) %in% xy)]
    seedDF[,which(colnames(dat) %in% siteId)] <- seeding[,which(colnames(seeding) %in% siteId)]
    seedDat <- rbind(dat, seedDF)
    sites <- seedDat[,siteId]
    seeds <- seeding[,siteId]
    datSV <- terra::vect(seedDat, geom = xy, crs = crs)
    posPools <- lapply(seeds, function(s) {
      seedRow <- which(seedDat[, siteId] == s)[1]
      sPool <- findPool2(seedRow, datSV, sites, xy, r, crs)
      sPool <- sPool[-which(sPool %in% seeds)]
      n <- length(sPool)
      if (n >= nSite)
        sPool
    })
    names(posPools) <- seeds
  }
  posPools <- Filter(Negate(is.null), posPools)
  if(length(posPools) > 1) {
    if(oType == "area"){
      if(is.null(seeding)){
        datSVSub <- terra::vect(dat[dat[,siteId] %in% names(posPools),], geom = xy, crs = crs)
      } else {
        datSVSub <- terra::vect(seeding[seeding[,siteId] %in% names(posPools),], geom = xy, crs = crs)
      }
      posPoolsDM <- terra::distance(datSVSub, datSVSub)
      a <- pi*(r^2)
      overlap <- apply(posPoolsDM, c(1,2), function(x) getOverlap(d = x, r = r, a = a))
      diag(overlap) <- 0
      if(oThreshold > 0) {
        OCs <- apply(overlap >= oThreshold, 1, function(x) length(which(x)))
      } else {
        OCs <- apply(overlap > oThreshold, 1, function(x) length(which(x)))
      }
    } else {
      if(oType == "sites"){
        index <- expand.grid(seq(1,length(posPools),1),seq(1,length(posPools),1))
        ShaSi <- matrix(0, length(posPools), length(posPools))
        nSi <- sapply(posPools, function(x) length(x))
        for (i in 1:nrow(index)){
          ShaSi[index[i,1],index[i,2]] <- length(intersect(posPools[[index[i,1]]],posPools[[index[i,2]]]))
        }
        overlap <- t(sapply(1:nrow(ShaSi), function(x) ShaSi[x,]/nSi[x]))
        diag(overlap) <- 0
        if(oThreshold > 0) {
          OCs <- apply(overlap >= oThreshold, 1, function(x) length(which(x)))
        } else {
          OCs <- apply(overlap > oThreshold, 1, function(x) length(which(x)))
        }
      } else {
        stop("Argument oType needs to be 'area' or 'sites'")
      }
    }
    if(sum(OCs) > 0){
      if(oPruningMode == "minOverlap"){
        keepers <- 1:length(OCs)
        while(T){
          max.o <- which(OCs == max(OCs))
          if(length(max.o) > 1){
            ta <- apply(overlap[max.o,], 1, function(x) sum(x))
            max.o <- max.o[sample(which(ta == max(ta)),1)]
          }
          overlap <- overlap[-max.o,-max.o]
          if(oThreshold > 0) {
            OCs <- apply(overlap >= oThreshold, 1, function(x) length(which(x)))
          } else {
            OCs <- apply(overlap > oThreshold, 1, function(x) length(which(x)))
          }
          keepers <- keepers[-max.o]
          if(sum(OCs) == 0){
            break
          }
        }
        finalPools <- posPools[keepers]
      } else {
        if(oPruningMode == "maxOccs"){
          ## initialise droppers vector, get number of occurrences for each potential circular subsample
          droppers <- c()
          posPools.occs <- sapply(1:length(posPools), function(x) nrow(rawData[which(rawData[,uniqID] %in% dat[which(dat[,siteId] %in% posPools[[x]]),uniqID]),]))
          ## Now loop
          while(T){
            ## get each posPool element that overlaps with at least one other element
            potential.droppers <- which(OCs > 0)
            ## get number of occurrences for these overlapping elements
            potential.dropper.occs <- posPools.occs[potential.droppers]
            ## identify overlapping posPool element with least amount of data
            TBDropped.index <- which(potential.dropper.occs == min(potential.dropper.occs))
            ## if tied (i.e., more than 1), pick from selection randomly
            if(length(TBDropped.index) > 1){
              TBDropped.index <- sample(TBDropped.index, 1)
            }
            ## add selected posPool element to droppers vector
            droppers <- c(droppers, potential.droppers[TBDropped.index])
            ## convert rows and columns attributed to selected posPool element in overlap matrix to 0 (i.e, won't contribute in next round to OC value)
            overlap[potential.droppers[TBDropped.index],] <- 0
            overlap[,potential.droppers[TBDropped.index]] <- 0
            ## re-calculate number of times overlap threshold exceeded for each posPool element
            if(oThreshold > 0) {
              OCs <- apply(overlap >= oThreshold, 1, function(x) length(which(x)))
            } else {
              OCs <- apply(overlap > oThreshold, 1, function(x) length(which(x)))
            }
            ## if none of the remaining posPool elements exceed overlap threshold, break
            if(sum(OCs) == 0){
              break
            }
          }
          ## finally, drop posPool elements listed in droppers vector from final list of pools.
          finalPools <- posPools[-droppers]
        } else {
          stop("argument oPruningMode needs to be 'minOverlap' or 'maxOccs'")
        }
      }
    } else {
      finalPools <- posPools
    }
  } else {
    finalPools <- posPools
  }
  return(finalPools)
}
