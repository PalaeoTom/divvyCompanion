#' Find spatially-distinct clusters of radially-constrained subsamples that fulfil all specified criteria
#'
#' @noRd
findClusters <- function(dat, siteId, xy, r, nSite, crs = "EPSG:4326"){
  ## Get grid cell IDs
  sites <- dat[, siteId]
  ## Convert dat into SpatVector
  datSV <- terra::vect(dat, geom = xy, crs = crs)
  ## For each site
  posPools <- lapply(sites, function(s) {
    ## Find site
    seedRow <- which(dat[, siteId] == s)[1]
    ## Find all cells that fall within buffer
    sPool <- findPool2(seedRow, datSV, sites, xy, r, crs)
    ## Get number of cells
    n <- length(sPool)
    ## If number equal to or greater than minimum, return RCR. Otherwise, result is null
    if (n >= nSite)
      sPool
  })
  ## Name list of all potential RCRs
  names(posPools) <- sites
  ## Drop those that include fewer than nSite cells
  posPools <- Filter(Negate(is.null), posPools)
  ## If more than one viable RCR, need to assign to clusters
  if(length(posPools) > 1) {
    ## Index of all comparisons that need to be made between RCRs
    index <- expand.grid(seq(1,length(posPools),1),seq(1,length(posPools),1))
    ## Name columns and add column to track proportion of cells in RCR 1 (seeded by s1) present in RCR 2 (seeded by s2)
    colnames(index) <- c("s1","s2")
    index$p_s1 <- NA
    ## Get vector counting number of grid cells in each RCR
    nSi <- sapply(posPools, function(x) length(x))
    ## For each pair of RCRs to be compared
    for (r in 1:nrow(index)){
      ## Count the number of grid cells that overlap between two cookies and get proportion
      index[r,3] <- length(intersect(posPools[[index[r,1]]],posPools[[index[r,2]]]))/nSi[index[r,1]]
    }
    ## Drop all rows with matching s1 and s2
    index <- index[-which(apply(index, 1, function(x) x[1]==x[2])),]
    ## Intialise clusters
    clusters <- list()
    ## tracker for cell incorporation
    tracker <- rep(NA,length(posPools))
    names(tracker) <- names(posPools)
    ## For each RCR
    for(i in 1:length(posPools)){
    #for(i in 1:57){
      ## Isolate relevant data
      sub <- index[which(index[,1]==i),]
      ## Get overlapping RCRs
      ov <- sub[which(sub[,3] > 0),2]
      #tracker[ov]
      ## Check to see if RCR has already been clustered.
      if(!is.na(tracker[i])){
        ## If it has, check that all constituent grid cells are assigned to the same cluster
        c <- unique(tracker[ov])
        c <- c[!is.na(c)]
        ## If c = 1, skip, as all is well. If not, merge clusters linked by grid cell.
        if(length(c)>1){
          ## Isolate clusters to be merged
          tbm <- clusters[c]
          clusters <- clusters[-c]
          ## Create new cluster
          clusters <- c(clusters, list(c(unlist(tbm,recursive = F), posPools[ov[is.na(tracker[ov])]])))
          ## Create new tracker
          tracker <- rep(NA,length(posPools))
          names(tracker) <- names(posPools)
          ## Populate new tracker
          for(t in 1:length(clusters)){
            tracker[names(tracker) %in% names(clusters[[t]])] <- t
          }
        } else {
          next
        }
      } else {
        ## If at least 1 RCR overlaps, assign overlapping RCRs to cluster. Otherwise, solo RCR cluster.
        if(length(ov) > 0){
          ## If any of overlapping RCRs have already been assigned to a cluster
          if(any(!is.na(tracker[ov]))){
          ## Add all unassigned RCRs present in ov to this cluster, as the already assigned RCR forms a bridge.
            # Get cluster id
            c <- unique(tracker[ov])
            c <- c[!is.na(c)]
            ## If more than one, merge these clusters, then re-do tracker
            if(length(c)>1){
              ### Isolate clusters to be merged
              tbm <- clusters[c]
              clusters <- clusters[-c]
              ## Create new cluster
              clusters <- c(clusters, list(c(unlist(tbm,recursive = F), posPools[c(i,ov[is.na(tracker[ov])])])))
              ## Create new tracker
              tracker <- rep(NA,length(posPools))
              names(tracker) <- names(posPools)
              ## Populate new tracker
              for(t in 1:length(clusters)){
                tracker[names(tracker) %in% names(clusters[[t]])] <- t
              }
            } else {
              ## Otherwise, add i and ov[is.na(ov)] to existing cluster
              clusters[[c]] <- c(clusters[[c]],posPools[c(i,ov[is.na(tracker[ov])])])
              ## Then update tracker
              tracker[c(i,ov)] <- c
            }
          } else {
            ## Otherwise, add as new cluster
            tracker[c(i,ov)] <- length(clusters)+1
            clusters <- c(clusters,list(posPools[c(i,ov)]))
          }
        } else {
          ## Add RCR to clusters as a single-RCR cluster and record in tracker number of latest cluster
          tracker[i] <- length(clusters)+1
          clusters <- c(clusters, list(posPools[i]))
        }
      }
    }
  } else {
    ## Just one viable RCR. Therefore there can only be one cluster.
    clusters <- list(posPools)
  }
  ## Name clusters
  names(clusters) <- paste0("clus",seq(1:length(clusters)))
  return(clusters)
}
