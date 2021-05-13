
#' doubledeepms__calculate_angle
#'
#' Calculate angle between vectors in degrees
#'
#' @param a vector a (required)
#' @param b vector b (required)
#'
#' @return angle in degrees
#' @export
doubledeepms__calculate_angle <- function(
  a,
  b){
	acos((a[1]*b[1] + a[2]*b[2] + a[3]*b[3])/(sqrt(a[1]^2 + a[2]^2 + a[3]^2) * sqrt(b[1]^2 + b[2]^2 + b[3]^2)))*180/pi
}
