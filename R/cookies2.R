#' Derive radially constrained spatial subsamples from occurrence data
#'
#' @param dataMat D
#' @param xy D
#' @param r D
#' @param seeding D
#' @param rarefaction D
#' @param reps D
#' @param nSites D
#' @param nOccs D
#' @param oThreshold D
#' @param oType D
#' @param oPruningMode D
#' @param returnSeeds D
#' @param crs D
#' @param output D
#'
#' @return D
#' @export
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
