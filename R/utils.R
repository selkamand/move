#' Normalize a numeric vector
#'
#' Scales a vector to unit length by dividing each element by its Euclidean norm.
#'
#' @param x Numeric vector.
#' @returns A numeric vector of the same length as `x`, scaled so that its magnitude is 1.
#' @examples
#' normalise(c(3, 4))  # returns c(0.6, 0.8)
#' @export
normalise <- function(x) {
  x / sqrt(sum(x^2))
}

#' Compute the magnitude (Euclidean norm) of a vector
#'
#' Calculates the length of a numeric vector using the Euclidean (L2) norm.
#'
#' @param x Numeric vector.
#' @returns Numeric scalar representing the magnitude of `x`.
#' @examples
#' magnitude(c(3, 4))  # returns 5
#' @export
magnitude <- function(x){
  sqrt(sum(x^2))
}

#' Clamp values to a specified range
#'
#' Restricts numeric values to lie between given minimum and maximum bounds.
#'
#' @param x Numeric vector.
#' @param min,max Numeric scalars specifying the lower and upper limits.
#' @returns A numeric vector where all values lie within `[min, max]`.
#' @examples
#' clamp(c(-2, 0.5, 3), 0, 1)
#' @export
clamp <- function(x, min = -1, max = 1){
  pmax(min, pmin(max, x))
}

#' Check if two 3D vectors are parallel
#'
#' Determines whether two 3D vectors are parallel or anti-parallel
#' based on the magnitude of their cross product.
#'
#' @param v,t Numeric length-3 vectors.
#' @param tol Numeric tolerance for detecting near-zero values.
#' @returns Logical; `TRUE` if `v` and `t` are (anti-)parallel within tolerance.
#' @examples
#' is_parallel(c(1,0,0), c(2,0,0))      # TRUE
#' is_parallel(c(1,0,0), c(-1,0,0))     # TRUE
#' is_parallel(c(1,0,0), c(0,1,0))      # FALSE
#' @export
is_parallel <- function(v, t, tol = 1e-8) {
  if (sqrt(sum(v^2)) < tol || sqrt(sum(t^2)) < tol)
    stop("Cannot check parallelism for zero-length vectors")

  cross <- pracma::cross(v, t)
  norm_cross <- sqrt(sum(cross^2))
  norm_cross < tol
}

#' Check if two 3D vectors are anti-parallel
#'
#' Tests whether two 3D vectors point in opposite directions (180° apart).
#'
#' @param v,t Numeric length-3 vectors.
#' @param tol Numeric tolerance for floating-point comparisons.
#' @returns Logical; `TRUE` if `v` and `t` are anti-parallel within tolerance.
#' @examples
#' is_antiparallel(c(1,0,0), c(-1,0,0))  # TRUE
#' is_antiparallel(c(1,0,0), c(1,0,0))   # FALSE
#' @export
is_antiparallel <- function(v, t, tol = 1e-8) {
  if (length(v) != 3L || length(t) != 3L)
    stop("Both v and t must be numeric vectors of length 3.")

  nv <- sqrt(sum(v^2))
  nt <- sqrt(sum(t^2))

  if (nv < tol || nt < tol)
    stop("Cannot test anti-parallelism for zero-length vectors.")

  v_unit <- v / nv
  t_unit <- t / nt

  cos_theta <- sum(v_unit * t_unit)
  cos_theta <- max(-1, min(1, cos_theta))

  abs(cos_theta + 1) < tol
}

#' Classify directional relationship between two vectors
#'
#' Returns whether two vectors are \emph{parallel}, \emph{anti-parallel}, or \emph{neither}
#' based on their normalized dot product.
#'
#' @param v,t Numeric vectors of equal length.
#' @param tol Numeric tolerance for numerical comparison.
#' @returns Character string: `"parallel"`, `"antiparallel"`, or `"neither"`.
#' @examples
#' vector_alignment(c(1,0,0), c(1,0,0))    # "parallel"
#' vector_alignment(c(1,0,0), c(-1,0,0))   # "antiparallel"
#' vector_alignment(c(1,0,0), c(0,1,0))    # "neither"
#' @export
vector_alignment <- function(v, t, tol = 1e-8) {
  if (sqrt(sum(v^2)) < tol || sqrt(sum(t^2)) < tol)
    stop("Zero-length vector")

  v_unit <- v / sqrt(sum(v^2))
  t_unit <- t / sqrt(sum(t^2))
  cos_theta <- sum(v_unit * t_unit)
  cos_theta <- max(-1, min(1, cos_theta))

  if (abs(cos_theta - 1) < tol) return("parallel")
  if (abs(cos_theta + 1) < tol) return("antiparallel")
  "neither"
}
