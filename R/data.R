#' Benzene molecule atomic coordinates
#'
#' Atomic coordinate data for a benzene molecule (C\eqn{_6}H\eqn{_6}).
#' Contains atom identifiers, element types, coordinates, and associated
#' metadata suitable for use with molecular visualisation or geometry functions.
#'
#' @format A data frame with 12 rows and 10 columns:
#' \describe{
#'   \item{eleno}{Integer atom ID.}
#'   \item{elena}{Atom element symbol.}
#'   \item{x, y, z}{Cartesian coordinates.}
#'   \item{elety}{Atom type (e.g., `"C.ar"` for aromatic carbon).}
#'   \item{resno}{Residue number.}
#'   \item{resid}{Residue identifier.}
#'   \item{charge}{Partial atomic charge.}
#'   \item{statbit}{Optional status field (may be `NA`).}
#' }
#'
#' @examples
#' data(benzene)
#' head(benzene)
#'
#' @source Generated from internal molecular coordinate data.
#' @keywords datasets
"benzene"
