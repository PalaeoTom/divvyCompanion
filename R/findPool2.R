#' Find the other populated cells that fall within a circle centered on a cell.
#'
#' @noRd
findPool2 <- function (seedRow, datSV, sites, xy, r, crs){
  seedpt <- datSV[seedRow, ]
  buf <- terra::buffer(seedpt, width = r)
  if (crs != "EPSG:4326") {
    buf <- terra::project(buf, y = "EPSG:4326")
    datSV <- terra::project(datSV, y = "EPSG:4326")
  }
  poolBool <- terra::extract(buf, datSV)
  pool <- sites[which(!is.na(poolBool[, 2]))]
  return(pool)
}
