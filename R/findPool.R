#' Find the other populated cells that fall within a circle centered on a cell.
#'
#' @noRd
findPool <- function (seedRow, coordsSV, sites, xy, r, crs){
  seedpt <- coordsSV[seedRow, ]
  buf <- terra::buffer(seedpt, width = r)
  if (crs != "EPSG:4326") {
    buf <- terra::project(buf, y = "EPSG:4326")
    coordsSV <- terra::project(coordsSV, y = "EPSG:4326")
  }
  poolBool <- terra::extract(buf, coordsSV)
  pool <- sites[which(!is.na(poolBool[, 2]))]
  return(pool)
}
