# A series of calculations used in computational chemistry

#' Compute dihedral statistics for A-B-C-D molecular geometry
#'
#' Calculates bond angle, bond length, and torsion angle for a four-atom molecular
#' geometry A-B-C-D. The torsion angle represents the dihedral angle between
#' planes defined by atoms A-B-C and B-C-D.
#'
#' @param a Numeric vector of length 3 giving the xyz coordinates of atom A.
#' @param b Numeric vector of length 3 giving the xyz coordinates of atom B.
#' @param c Numeric vector of length 3 giving the xyz coordinates of atom C.
#' @param d Numeric vector of length 3 giving the xyz coordinates of atom D.
#'
#' @return A list containing:
#' \describe{
#'   \item{a, b, c}{The input atomic coordinates (unchanged).}
#'   \item{bond_angle}{Numeric scalar giving the B-C-D bond angle in degrees.}
#'   \item{bond_length}{Numeric scalar giving the C-D bond length.}
#'   \item{torsion_angle}{Numeric scalar giving the A-B-C-D torsion angle in degrees,
#'     or \code{NA} if the atoms are colinear.}
#' }
#'
#' @details
#' The torsion angle is calculated as the signed angle between the planes
#' defined by atoms A-B-C and B-C-D, using the B-C bond as the reference axis.
#' When atoms A, B, C, D are colinear or nearly colinear, the torsion angle
#' is undefined and the function returns \code{NA} with a warning.
#'
#' @examples
#' # Normal tetrahedral-like geometry
#' a <- c(0, 0, 0)
#' b <- c(1, 0, 0)
#' c <- c(1, 1, 0)
#' d <- c(1, 1, 1)
#' compute_abcd_dihedral_stats(a, b, c, d)
#'
#' # Colinear case (returns NA for torsion angle)
#' a <- c(0, 0, 0)
#' b <- c(1, 0, 0)
#' c <- c(2, 0, 0)
#' d <- c(3, 0, 0)
#' compute_abcd_dihedral_stats(a, b, c, d)
#'
#' @seealso \code{\link{compute_abcd_bond_angle}} for computing atomic positions
#'   from dihedral parameters.
#'
#' @export
compute_abcd_dihedral_stats <- function(a, b, c, d) {
  # Assertions
  if (length(a) != 3) {
    stop("`a` must be a 3d vector representing position of atom A in A-B-C-D(dummy) dihedral. Length should be 3, not: ", length(a))
  }
  if (length(b) != 3) {
    stop("`b` must be a 3d vector representing position of atom B in A-B-C-D(dummy) dihedral. Length should be 3, not: ", length(b))
  }
  if (length(c) != 3) {
    stop("`c` must be a 3d vector representing position of atom C in A-B-C-D(dummy) dihedral. Length should be 3, not: ", length(c))
  }
  if (length(d) != 3) {
    stop("`d` must be a 3d vector representing position of atom D in A-B-C-D(dummy) dihedral. Length should be 3, not: ", length(d))
  }

  # Bond Angle BCD (in degrees)
  cb <- create_vector_from_start_end(start = c, end = b)
  cd <- create_vector_from_start_end(start = c, end = d)
  bond_angle <- measure_angle_between_vectors(cb, cd, degrees = TRUE)


  # Bond Length (CD)
  bond_length <- measure_distance_between_two_points(c, d)

  # Torsion Angle (angle between A and D along axis BC). Sometimes called the dihedral angle.
  # This is actually the angle between planes created by A-B|B-C and BC|CD (in degrees)
  bc <- create_vector_from_start_end(b, c)

  # Check for colinear planes (when BA and BC are parallel, or CB and CD are parallel)
  ba <- create_vector_from_start_end(b, a)
  cb <- create_vector_from_start_end(c, b)
  cd <- create_vector_from_start_end(c, d)

  # Check if vectors are parallel (cross product magnitude near zero)
  cross_abc_magnitude <- magnitude(pracma::cross(ba, bc))
  cross_bcd_magnitude <- magnitude(pracma::cross(cb, cd))

  if (cross_abc_magnitude < 1e-8 || cross_bcd_magnitude < 1e-8) {
    warning("Points A, B, C, D are colinear or nearly colinear. Torsion angle is undefined.")
    torsion_angle <- NA_real_
  } else {
    normal_abc <- compute_plane_normal_from_vectors(ba, bc)
    normal_bcd <- compute_plane_normal_from_vectors(cb, cd)
    torsion_angle <- measure_signed_angle_between_planes(normal_abc, normal_bcd, ref_axis = bc, degrees = TRUE)
  }


  # Return List of Paramaters
  return(
    list(
      a = a,
      b = b,
      c = c,
      bond_angle = bond_angle,
      bond_length = bond_length,
      torsion_angle = torsion_angle
    )
  )
}

#' Compute fourth atom position or bond vector from dihedral parameters
#'
#' Calculates either the bond vector (relative to atom C) or absolute position of atom D
#' in an A-B-C-D molecular geometry given the positions of atoms A, B, C and the
#' desired bond angle, bond length, and torsion angle parameters.
#'
#' @param a Numeric vector of length 3 giving the xyz coordinates of atom A.
#' @param b Numeric vector of length 3 giving the xyz coordinates of atom B.
#' @param c Numeric vector of length 3 giving the xyz coordinates of atom C.
#' @param bond_angle Numeric scalar giving the B-C-D bond angle in degrees.
#' @param bond_length Numeric scalar giving the C-D bond length.
#' @param torsion_angle Numeric scalar giving the A-B-C-D torsion angle in degrees.
#' @param return_absolute_position Logical; if \code{TRUE}, returns the absolute
#'   coordinates of atom D. If \code{FALSE} (default), returns the bond vector
#'   (relative to atom C) that when added to C gives the position of D.
#'
#' @return A numeric vector of length 3 giving either:
#' \describe{
#'   \item{If \code{return_absolute_position = FALSE} (default)}{The bond vector from C to D}
#'   \item{If \code{return_absolute_position = TRUE}}{The absolute coordinates of atom D}
#' }
#'
#' @details
#' This function uses the \code{compas::calCo} function to compute the position
#' of atom D based on the three preceding atoms (A, B, C) and the specified
#' geometric parameters. By default, it returns the bond vector (D - C) which
#' is often more useful for molecular modeling applications.
#'
#' @examples
#' # Define first three atoms
#' a <- c(0, 0, 0)
#' b <- c(1, 0, 0)
#' c <- c(1, 1, 0)
#'
#' # Compute bond vector (default behavior)
#' bond_vector <- compute_fourth_atom_position(a, b, c, bond_angle = 109.5, bond_length = 1.5, torsion_angle = 60)
#' print(bond_vector)
#'
#' # Get absolute position of fourth atom
#' d_position <- compute_fourth_atom_position(a, b, c, bond_angle = 109.5, bond_length = 1.5, torsion_angle = 60, return_absolute_position = TRUE)
#' print(d_position)
#'
#' # Verify: bond_vector + c should equal d_position
#' print(c + bond_vector)
#'
#' @seealso \code{\link{compute_abcd_dihedral_stats}} for computing dihedral
#'   parameters from atomic positions.
#'
#' @export
compute_fourth_atom_position <- function(a, b, c, bond_angle, bond_length, torsion_angle, return_absolute_position = FALSE) {
  rlang::check_installed("compas", reason = "to compute bond angles in `compute_fourth_atom_position()`")

  prev_atoms <- matrix(data = c(a, b, c), byrow = TRUE, ncol = 3)
  d_position <- compas::calCo(prev_atoms = prev_atoms, length = bond_length, bAngle = bond_angle, tAngle = torsion_angle)

  if (return_absolute_position) {
    return(d_position)
  } else {
    # Return bond vector (D - C)
    return(d_position - c)
  }
}
