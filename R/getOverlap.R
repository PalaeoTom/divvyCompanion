#' Efficient calculation of percentage of areas of two identical circles that overlap, given their radii, their distance between them, and their circumference
#'
#' @param d D
#' @param r D
#' @param a D
#'
#' @noRd
#'
#' @return D
#' @export
#'
#' @examples
#' print("example")
getOverlap <- function(d, r, a){
  if(d >= (2*r)){
    o <- 0
  }
  if(d == 0){
    o <- 1
  }
  if(d > 0 && d < (2*r)){
    o <- (2*(((r^2)*acos(d/(2*r)))-((d/4)*(sqrt((4*(r^2))-(d^2))))))/a
  }
  return(o)
}
