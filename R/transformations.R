# Rotation ----------------------------------------------------------------
#' Rotate a vector around an axis by a given angle
#'
#' Rotates a 3D vector \code{v} about \code{rotation_axis} by \code{angle} radians
#' using Rodrigues' rotation formula. Positive angles follow the right-hand rule.
#'
#' @param v Numeric length-3 vector to rotate.
#' @param rotation_axis Numeric length-3 vector; rotation axis (will be normalized).
#' @param angle Numeric scalar (radians); rotation angle.
#' @param tol Numeric; tolerance for early return check (uses sum(abs(v)) < tol).
#' @param zap zap small numbers to zero. digits will be based on getOption("digits")
#' @return Rotated vector of same length as \code{v}.
#'
#' @note
#' - Early return uses \code{sum(abs(v)) < tol}, which may skip rotation for some nonzero vectors.
#' - Depends on user-defined \code{normalise()} and \code{pracma::cross()}.
#' - No input validation; assumes 3-element numeric vectors and nonzero rotation axis.
#'
#' @examples
#' # Rotate (1,0,0) by 90° around z-axis -> (0,1,0)
#' v <- c(1, 0, 0)
#' axis <- c(0, 0, 1)
#' rotate_vector_around_axis(v, axis, pi / 2)
#'
#' @export
rotate_vector_around_axis <- function(v, rotation_axis, angle, tol = 1e-8, zap = TRUE) {
  # Unname inputs
  v <- unname(v)
  rotation_axis <- unname(rotation_axis)

  # If zero return untransformed
  if (sum(abs(v)) < tol) {
    return(v)
  }

  # If angle = 0 return untransformed
  if (abs(angle) < tol) {
    return(v)
  }

  # rotation axis must be non-zero
  if (sqrt(sum(rotation_axis^2)) < tol) {
    stop("`rotation_axis` must be non-zero.")
  }

  # Get stats required for rodrigues rotation
  sintheta <- as.vector(sin(angle))
  costheta <- as.vector(cos(angle))

  # normalise the rotation axis
  n <- normalise(rotation_axis)

  # Rotate vector
  v_new <- v * costheta + (pracma::cross(n, v)) * sintheta + n * as.vector(n %*% v) * (1 - costheta)

  # Return vector with very small numbers zapped to zero
  if (zap) v_new <- zapsmall(v_new)

  return(v_new)
}

#' Rotate a vector around an axis passing through an arbitrary point
#'
#' Rotates a 3D vector `v` about the line defined by `rotation_axis`
#' and a point `point_on_axis`. The rotation axis need not pass through the origin.
#'
#' @param v Numeric length-3 vector representing the position to rotate.
#' @param rotation_axis Numeric length-3 vector giving the axis direction (will be normalized).
#' @param point_on_axis Numeric length-3 vector giving a point lying on the rotation axis.
#' @param angle Numeric scalar (radians); rotation angle (positive follows right-hand rule).
#' @param tol Numeric tolerance used for degeneracy checks (passed to `rotate_vector_around_axis`).
#' @param zap Logical; if `TRUE`, small numerical noise is removed via `zapsmall()`.
#'
#' @return Numeric length-3 vector giving the rotated position.
#'
#' @details
#' The algorithm translates `v` so the rotation axis passes through the origin,
#' applies [rotate_vector_around_axis()], then translates back.
#'
#' @examples
#' v <- c(2, 0, 0)
#' axis <- c(0, 0, 1)
#' point <- c(1, 0, 0)
#' rotate_vector_around_axis_through_point(v, axis, point, pi / 2)
#' # Returns c(1, 1, 0)
#'
#' @export
rotate_vector_around_axis_through_point <- function(v,
                                                    rotation_axis,
                                                    point_on_axis,
                                                    angle,
                                                    tol = 1e-8,
                                                    zap = TRUE) {
  if (!is.numeric(v) || length(v) != 3L) stop("`v` must be a numeric vector of length 3.")
  if (!is.numeric(rotation_axis) || length(rotation_axis) != 3L) stop("`rotation_axis` must be a numeric vector of length 3.")
  if (!is.numeric(point_on_axis) || length(point_on_axis) != 3L) stop("`point_on_axis` must be a numeric vector of length 3.")

  v <- unname(v)
  rotation_axis <- unname(rotation_axis)
  point_on_axis <- unname(point_on_axis)

  if (sqrt(sum(rotation_axis^2)) < tol) {
    stop("`rotation_axis` must be non-zero.")
  }

  v_centered <- v - point_on_axis
  rotated_centered <- rotate_vector_around_axis(
    v_centered,
    rotation_axis,
    angle,
    tol = tol,
    zap = FALSE
  )

  v_rotated <- rotated_centered + point_on_axis

  if (zap) {
    v_rotated <- zapsmall(v_rotated)
  }

  return(v_rotated)
}


#' Rotate a 3D vector to align with a target direction (Rodrigues' formula)
#'
#' Rotates a 3D vector \code{v} so that its \emph{direction} aligns with the
#' direction of \code{target}, using Rodrigues' rotation formula.
#' The returned vector preserves the original length of \code{v}.
#'
#' @details
#' Let \eqn{\mathbf{u} = \frac{\mathbf{v}}{\|\mathbf{v}\|}} and
#' \eqn{\mathbf{t} = \frac{\mathbf{target}}{\|\mathbf{target}\|}} be unit vectors.
#' Define \eqn{\theta = \arccos(\mathbf{u}\cdot\mathbf{t})} and rotation axis
#' \eqn{\mathbf{k} = \frac{\mathbf{u}\times\mathbf{t}}{\|\mathbf{u}\times\mathbf{t}\|}} (if not degenerate).
#'
#' Rodrigues' rotation formula rotates any vector \eqn{\mathbf{x}} about unit axis \eqn{\mathbf{k}}
#' by angle \eqn{\theta} via:
#' \deqn{
#'   \mathrm{Rot}_{\mathbf{k},\theta}(\mathbf{x}) =
#'   \mathbf{x}\cos\theta + (\mathbf{k}\times\mathbf{x})\sin\theta + \mathbf{k}\,(\mathbf{k}\cdot\mathbf{x})(1-\cos\theta).
#' }
#' This function applies the formula to \eqn{\mathbf{v}} with \eqn{\mathbf{k}} chosen to rotate
#' \eqn{\mathbf{u}} toward \eqn{\mathbf{t}} (right-hand rule).
#'
#' Edge cases handled:
#' \itemize{
#'   \item \strong{Already aligned} (\eqn{\theta \approx 0}): returns \code{v} unchanged.
#'   \item \strong{Opposite directions} (\eqn{\theta \approx \pi}): chooses a stable axis orthogonal to \code{v}
#'         and applies the \eqn{180^\circ} rotation (closed form).
#'   \item \strong{Numerical safety}: clamps the dot product to \code{[-1, 1]} before \code{acos}.
#' }
#'
#' @param v Numeric length-3 vector. The vector to rotate (need not be unit length).
#' @param target Numeric length-3 vector. Direction to align \code{v} with (need not be unit;
#'        only direction matters).
#' @param tol Numeric scalar. Tolerance for detecting degeneracies (parallel/opposite/zero norms).
#' @param return Character; either \code{"v_new"} (default) to return the rotated vector,
#'   or \code{"axis_plus_angle"} to return a list with rotation axis and angle.
#'
#' @returns
#' Depending on the value of `return`, the function returns:
#'
#' \strong{1. Rotated vector (`return = "v_new"`, default)}
#' A numeric length-3 vector giving the rotated version of `v`.
#' The output satisfies:
#' \itemize{
#'   \item its \emph{direction} is aligned with `target`,
#'   \item its \emph{length} is unchanged (`||v_new|| = ||v||`),
#'   \item it is numerically identical to `v` when no rotation is required.
#' }
#'
#' \strong{2. Rotation parameters (`return = "axis_plus_angle"`)}
#' A list containing:
#' \describe{
#'   \item{axis}{A length-3 \emph{unit} vector describing the rotation axis
#'         (right-hand rule).}
#'   \item{angle}{The rotation angle in radians required to align `v` with `target`.}
#' }
#' These values can be applied to any other 3D vector using
#' `rotate_vector_around_axis()`, making it easy to transform multiple
#' points consistently.
#'
#' In either mode, the function handles aligned, anti-aligned, and
#' numerically degenerate cases robustly (e.g., 180° flips, tiny cross products).
#'
#' @section Potential pitfalls:
#' \itemize{
#'   \item \strong{Zero vectors:} If \code{v} or \code{target} has (near) zero norm,
#'         rotation is undefined;
#'   \item \strong{Axis normalization:} The rotation axis must be a \emph{unit} vector. We automatically scale non-unit target vectors
#'   \item \strong{180° case:} When \code{v} and \code{target} are opposite, \code{u × t = 0} and the axis
#'         is not defined. We choose a stable axis orthogonal to \code{v} automatically.
#' }
#'
#' @references
#' Rodrigues, O. (1840). \emph{Des lois géométriques qui régissent les déplacements d’un système solide dans l’espace}.
#'
#' @seealso
#' \code{rotations::genR()} (CRAN) for generating rotation matrices via Rodrigues’ formula.
#'
#' @examples
#' # 1) Simple 90° rotation around z: (1,0,0) -> (0,1,0) when target = (0,1,0)
#' rotate_vector_to_align_with_target(c(1, 0, 0), c(0, 1, 0))
#'
#' # 2) Align arbitrary v with (1,1,1): direction matches, length preserved
#' v <- c(2, -3, 4)
#' target <- c(1, 1, 1)
#' v_rot <- rotate_vector_to_align_with_target(v, target)
#' sqrt(sum(v^2)) # original length
#' sqrt(sum(v_rot^2)) # same length
#' cor(v_rot, target) # ~ 1 (directions aligned)
#'
#' # 3) Opposite directions (180°): target is -v's direction
#' v <- c(1, 2, 3)
#' target <- -v
#' rotate_vector_to_align_with_target(v, target)
#'
#' # 4) Input validation
#' \dontrun{
#' rotate_vector_to_align_with_target(c(0, 0, 0), c(1, 0, 0)) # error: zero-norm v
#' rotate_vector_to_align_with_target(c(1, 0, 0), c(0, 0, 0)) # error: zero-norm target
#' }
#'
#' @export
rotate_vector_to_align_with_target <- function(v, target, tol = 1e-8, return = c("v_new", "axis_plus_angle")) {
  return <- rlang::arg_match(return)

  if (!is.numeric(v) || length(v) != 3L) {
    stop("`v` must be a numeric vector of length 3.")
  }
  if (!is.numeric(target) || length(target) != 3L) {
    stop("`target` must be a numeric vector of length 3.")
  }

  if (sum(abs(v)) < tol) {
    warning("Cannot rotate a zero vector")
    return(NA)
  }
  if (sum(abs(target)) < tol) {
    warning("Cannot align to a zero target vector")
    return(NA)
  }

  # Unname inputs
  v <- unname(v)
  target <- unname(target)

  # Ensure target is a unit vector
  target_unit <- normalise(target)
  v_unit <- normalise(v)

  # Compute rotation axis and angle components
  crossproduct <- pracma::cross(v_unit, target_unit)
  cross_norm <- sqrt(sum(crossproduct^2))
  n <- if (cross_norm < tol) c(NA_real_, NA_real_, NA_real_) else crossproduct / cross_norm

  # Find the rotation angle
  costheta <- as.vector(v_unit %*% target_unit)

  # Clamp result to -1 to 1. It should already be restricted but floating point errors can cause values very slightly outside this range
  costheta <- clamp(costheta, min = -1, max = 1)

  # Rotation in radians
  theta <- acos(costheta)

  # Parallel case (theta ~ 0)
  if (abs(theta) < tol) {
    if (return == "axis_plus_angle") {
      return(list("axis" = c(1, 0, 0), "angle" = 0))
    } else {
      return(v)
    }
  }
  # Anti-parallel case (theta ~ pi): choose robust axis orthogonal to v
  if (abs(theta - pi) < tol) {
    basis <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
    dots <- as.numeric(basis %*% v_unit)
    idx <- which.min(abs(dots))
    cand <- pracma::cross(v_unit, basis[idx, ])
    if (sqrt(sum(cand^2)) < tol) {
      idx <- if (idx == 1) 2 else 1
      cand <- pracma::cross(v_unit, basis[idx, ])
    }
    n <- normalise(cand)
  }

  # Rotation in degrees (theta * 180/pi) if requested by caller elsewhere

  # Early return of paramaters
  if (return == "axis_plus_angle") {
    return(list("axis" = n, "angle" = theta))
  }

  # Apply rodrigues transformation (rotation around axis)
  v_new <- rotate_vector_around_axis(v = v, angle = theta, rotation_axis = n)

  return(v_new)
}

#' Compute the rotation needed to align a 3D vector with a plane
#'
#' Computes the minimal rotation required to move a 3D vector `v` so that it lies
#' within the plane defined by its normal vector `plane_normal`. The rotation
#' preserves the magnitude of `v`, unlike a projection which shortens it.
#'
#' @param v numeric(3) Vector to rotate.
#' @param plane_normal numeric(3) Plane normal (not necessarily unit length).
#' @param return whether to return the new (transformed) vector or a list of paramaters required to perform the rotation (axis and angle).
#' @return either a rotated vector or a list with elements `axis` (unit vector) and `angle` (radians) depending on value of `return` param.
#' @examples
#' rotate_vector_into_a_plane(c(1, 2, 3), c(0, 1, -1))
#' @export
rotate_vector_into_a_plane <- function(v, plane_normal, return = c("v_new", "axis_plus_angle")) {
  # Argument Checks
  return <- rlang::arg_match(return)

  # Ensure both vectors are numeric, no names attached
  v <- unname(v)
  plane_normal <- unname(plane_normal)

  # Normalise plane normal
  if (sqrt(sum(plane_normal^2)) == 0) {
    stop("`plane_normal` must be non-zero.")
  }
  plane_normal_norm <- normalise(plane_normal)

  # --- Compute in-plane component of v ---
  # Convert dot product to numeric scalar to avoid array recycling warning
  v_dot_n <- as.numeric(v %*% plane_normal_norm)
  v_inplane <- v - v_dot_n * plane_normal_norm

  # --- Compute rotation axis ---
  # Axis perpendicular to both v and its in-plane projection
  u <- pracma::cross(v, v_inplane)
  u_norm <- sqrt(sum(u^2))
  k <- if (u_norm == 0) c(1, 0, 0) else normalise(u) # unit rotation axis
  u_l2 <- sqrt(sum(u^2)) # magnitude of u (for angle)

  # --- Compute rotation angle (radians) ---
  theta <- atan2(u_l2, as.numeric(v %*% v_inplane))

  # Return rotation axis and angle
  if (return == "axis_plus_angle") {
    return(list(axis = k, angle = theta))
  }

  # Apply rodriguez rotation
  v_new <- rotate_vector_around_axis(v, rotation_axis = k, angle = theta)

  return(v_new)
}

# Projection --------------------------------------------------------------
#' Project one vector onto another
#'
#' Computes the vector projection of `a` onto `b`, i.e. the component of `a`
#' that lies along the direction of `b`.
#'
#' @param a,b Numeric vectors of equal length.
#' @returns A numeric vector parallel to `b`.
#' @examples
#' project_vector_into_vector(c(2, 3, 4), c(1, 0, 0)) # -> c(2,0,0)
#'
#' @export
project_vector_into_vector <- function(a, b) {
  if (length(a) != length(b)) stop("Vectors must have the same length.")
  if (sum(b^2) == 0) stop("Cannot project onto a zero-length vector.")

  scalar_coeff <- sum(a * b) / sum(b * b) # (a·b)/(b·b)
  projection <- scalar_coeff * b # multiply scalar by vector
  return(projection)
}

#' Compute the scalar projection (signed length) of one vector onto another
#'
#' Calculates how much of vector `a` lies along vector `b` (positive if same
#' direction, negative if opposite).
#'
#' @param a,b Numeric vectors of equal length.
#' @returns Numeric scalar giving the signed magnitude of the projection.
#' @examples
#' compute_scalar_projection(c(2, 3, 4), c(1, 0, 0)) # -> 2
#'
#' @export
compute_scalar_projection <- function(a, b) {
  if (length(a) != length(b)) stop("Vectors must have the same length.")
  if (sum(b^2) == 0) stop("Cannot project onto a zero-length vector.")

  scalar_proj <- sum(a * b) / sqrt(sum(b^2)) # (a·b)/‖b‖
  return(scalar_proj)
}


#' Project a vector into a plane (remove normal component)
#'
#' Returns the component of vector `v` that lies within the plane whose normal is
#' `plane_normal`. The normal need not be unit length.
#'
#' Mathematically: \code{v_in_plane = v - ((v·n)/(n·n)) n}.
#'
#' @param v Numeric vector: the vector to project.
#' @param plane_normal Numeric vector: the plane's normal (same length as `v`; need not be unit).
#'
#' @returns A numeric vector (same length as `v`) representing the projection of `v` into the plane.
#'
#' @details
#' - If `plane_normal` is unit length, the term \code{(v·n)/(n·n)} simplifies to \code{v·n}.
#' - This is a rigid *projection* (not a rotation); the returned vector generally has a smaller magnitude than `v`.
#'
#' @examples
#' # 3D example: drop the z-component (plane is xy-plane, normal = (0,0,1))
#' project_vector_into_plane(c(1, 2, 3), c(0, 0, 1))
#' # -> c(1, 2, 0)
#'
#' # Arbitrary plane normal (not unit):
#' v <- c(1, 2, 3)
#' n <- c(0, 1, -1) # plane normal
#' v_plane <- project_vector_into_plane(v, n)
#' # v decomposes as v_plane (in-plane) + v_normal (along n)
#'
#' @export
project_vector_into_plane <- function(v, plane_normal) {
  if (!is.numeric(v) || !is.numeric(plane_normal)) {
    stop("`v` and `plane_normal` must be numeric vectors.")
  }
  if (length(v) != length(plane_normal)) {
    stop("`v` and `plane_normal` must have the same length.")
  }
  nn <- sum(plane_normal^2)
  if (nn == 0) {
    stop("`plane_normal` must be non-zero.")
  }
  v - (sum(v %*% plane_normal) / nn) * plane_normal
}


# Translation -------------------------------------------------------------
#' Translate a position by a vector
#'
#' Shifts a position in space by adding a translation vector.
#' Both `position` and `vector` must have the same dimensionality.
#'
#' @param position Numeric vector giving the original coordinates.
#' @param vector Numeric vector giving the translation offset (same length as `position`).
#'
#' @returns A numeric vector representing the translated position.
#'
#' @examples
#' translate_position_by_vector(c(1, 2, 3), c(0.5, 0, -1))
#' # Returns: c(1.5, 2, 2)
#'
#' @export
translate_position_by_vector <- function(position, vector) {
  if (length(position) != length(vector)) {
    stop("To translate a position by a vector, both must have the same number of elements.")
  }
  position + vector
}

#' Translate a position along a direction by a given distance
#'
#' Moves a position along a specified direction vector by a given magnitude.
#' The direction vector is normalized internally, so only its orientation matters.
#'
#' @param position Numeric vector giving the starting position.
#' @param direction Numeric vector giving the direction of translation (need not be unit length).
#' @param magnitude Numeric scalar giving the distance to move along `direction`.
#'
#' @returns A numeric vector representing the new translated position.
#'
#' @examples
#' translate_position_in_direction(c(0, 0, 0), c(1, 1, 0), 10)
#' # Returns approximately c(7.07, 7.07, 0)
#'
#' @export
translate_position_in_direction <- function(position, direction, magnitude) {
  if (sqrt(sum(direction^2)) == 0) stop("`direction` must be non-zero.")
  translation_vector <- normalise(direction) * magnitude
  position + translation_vector
}


#' Compute the translation vector from one position to another
#'
#' @param position numeric vector; starting position.
#' @param target numeric vector; destination position (same length).
#' @return numeric vector: \code{target - position}.
#' @examples
#' compute_translation_vector(c(1, 2, 3), c(4, 6, 3)) # -> c(3,4,0)
#' @export
compute_translation_vector <- function(position, target) {
  if (length(position) != length(target)) {
    stop("Both positions must have the same number of elements.")
  }
  target - position
}


# Planes --------------------------------------------------------------

#' Unit normal of the plane spanned by two vectors
#' @param a,b numeric(3) spanning vectors (non-parallel).
#' @return numeric(3) unit normal perpendicular to both.
#' @examples
#' compute_plane_normal_from_vectors(c(1, 0, 0), c(0, 1, 0)) # c(0,0,1)
#' @export
compute_plane_normal_from_vectors <- function(a, b) {
  cp <- pracma::cross(a, b)
  if (sqrt(sum(cp^2)) == 0) stop("Input vectors must be non-parallel and non-zero.")
  normalise(cp)
}


#' Fit a plane to 3D points
#'
#' Fits a geometric plane to a set of 3D points given in columns `x`, `y`, and `z`.
#' For exactly three non-collinear points, the plane normal is computed analytically
#' from the cross product of two spanning vectors. For more than three points, a
#' least-squares best-fit plane is obtained via singular value decomposition (SVD)
#' of the centered coordinates.
#'
#' The fitted plane is of the form
#' \deqn{n \cdot x + d = 0}
#' where \code{n} is a unit normal vector and \code{d} is the scalar offset.
#'
#' @param points A data frame or matrix with numeric columns \code{x}, \code{y},
#'   and \code{z} giving point coordinates.
#' @param careful Logical; if \code{TRUE} (default), run basic input assertions
#'   (presence of \code{x}, \code{y}, \code{z} columns and at least one row).
#'
#' @return
#' Returns a list with components:
#' \describe{
#'   \item{normal}{Numeric length-3 unit normal vector (named \code{x}, \code{y}, \code{z}).}
#'   \item{offset}{Numeric scalar \code{d} in the plane equation \code{normal · x + d = 0}.}
#'   \item{centroid}{Numeric length-3 centroid of all points}
#' }
#'
#' @details
#' For \code{n > 3} points, the coordinates are first translated so that their
#' geometric center lies at the origin, and SVD is applied to the centered
#' \code{(x,y,z)} matrix. The plane normal is taken as the right singular vector
#' corresponding to the smallest singular value; the offset is then chosen so
#' that the plane passes through the original centroid.
#'
#' The function errors if fewer than three points are supplied.
#'
#' @examples
#' # Simple plane: z = 0 (points in the xy-plane)
#' pts <- data.frame(
#'   x = c(0, 1, 2, 0),
#'   y = c(0, 0, 1, 2),
#'   z = 0
#' )
#' compute_plane_from_points(pts)
#'
#' # Exactly three points defining a plane
#' pts3 <- data.frame(
#'   x = c(0, 1, 0),
#'   y = c(0, 0, 1),
#'   z = c(0, 0, 0)
#' )
#' compute_plane_from_points(pts3)
#'
#' @export
compute_plane_from_points <- function(points, careful = TRUE) {
  n_points <- nrow(points)

  if (careful) {
    if(!all(c("x", "y", "z") %in% colnames(points))) stop("points must be a data.frame/matrix with columns 'x', 'y' and 'z'")
    if(nrow(points) <= 0) stop("Cannot fit a plane to only [", n_points ,"]")
  }

  xyz <- points[, c("x", "y", "z"), drop = FALSE]
  center <- colMeans(xyz)

  if (n_points < 3) {
    stop("Cannot fit a plane to only [", n_points, "] points")
  }

  # If you only have 3 points then directly compute the plane
  if (n_points == 3) {
    a <- xyz[1, , drop = TRUE]
    b <- xyz[2, , drop = TRUE]
    c <- xyz[3, , drop = TRUE]

    ba <- create_vector_from_start_end(start = b, end = a)
    bc <- create_vector_from_start_end(start = b, end = c)

    plane <- compute_plane_normal_from_vectors(ba, bc)
    plane_position <- b
    plane_normal_offset <- convert_plane_point_normal_to_normal_offset(
      normal = plane,
      point = plane_position
    )
    plane_normal_offset$centroid <- center
    return(plane_normal_offset)
  }

  # If you have more than 3 points, use SVD to get the best fitting plane

  # Move geometric center to origin
  xyz[, "x"] <- xyz[,"x"] - center["x"]
  xyz[, "y"] <- xyz[,"y"] - center["y"]
  xyz[, "z"] <- xyz[,"z"] - center["z"]

  # Do singular value decomposition
  sv <- svd(xyz)

  # normal vector is last column of V (smallest singular value direction)
  normal <- sv$v[, ncol(sv$v)]
  normal <- normal / sqrt(sum(normal^2))
  names(normal) <- c("x", "y", "z")

  # plane constant d such that normal·x + d = 0
  d <- -sum(normal * center)

  return(list(normal = normal, offset = d, centroid = center))
}

#' Convert a plane from (point, normal) to (unit normal, offset)
#'
#' Converts a plane defined by a point lying on it and a (possibly non-unit) normal
#' into an equivalent representation defined by a **unit normal** and a **signed offset**
#' along that normal.
#'
#' The plane is expressed in **Hesse normal form**:
#' \deqn{\hat{n} \cdot x = s}
#' where \eqn{\hat{n}} is a unit normal vector and \eqn{s} is the signed distance
#' from the origin to the plane along \eqn{\hat{n}}.
#'
#' @param normal Numeric vector of length 3 — the plane's normal (does not need to be unit length).
#' @param point Numeric vector of length 3 — a known point lying on the plane.
#'
#' @return A list with:
#' \describe{
#'   \item{normal}{Unit normal vector (length 3).}
#'   \item{offset}{Signed distance from the origin along the normal.}
#' }
#'
#' @details
#' The offset is computed as \eqn{s = \hat{n} \cdot p}, where \eqn{p} is the supplied point
#' and \eqn{\hat{n}} = \eqn{normal / ||normal||}.
#' Positive values of \code{offset} mean the plane lies in the direction of the normal
#' from the origin.
#'
#' @examples
#' convert_plane_point_normal_to_normal_offset(c(0, 0, 2), c(2, 3, 5))
#' # $normal = c(0, 0, 1)
#' # $offset = 5
#'
#' # Relationship with rgl::planes3d():
#' x <- convert_plane_point_normal_to_normal_offset(c(1, 0, 0), c(3, 0, 0))
#' # rgl::planes3d(a = x$normal[1], b = x$normal[2], c = x$normal[3], d = -x$offset)
#' @export
convert_plane_point_normal_to_normal_offset <- function(normal, point) {
  stopifnot(length(normal) == 3L, length(point) == 3L)
  if (sqrt(sum(normal^2)) == 0) stop("`normal` must be non-zero.")
  n_hat <- normalise(normal)
  offset <- sum(n_hat * point) # s = n_hat · point
  list(normal = n_hat, offset = as.numeric(offset))
}


#' Convert a plane from (unit normal, offset) to (point, normal)
#'
#' Converts a plane defined by a **unit normal** and a **signed offset**
#' (Hesse normal form) into an equivalent representation defined by
#' a **point on the plane** and the same **unit normal**.
#'
#' The returned point is the one on the plane **closest to the origin**, given by:
#' \deqn{p = s \hat{n}}
#'
#' @param normal Numeric vector of length 3 — plane's unit normal.
#' @param offset Numeric scalar — signed distance from the origin along the normal.
#'
#' @return A list with:
#' \describe{
#'   \item{point}{A point on the plane (closest to the origin).}
#'   \item{normal}{Unit normal vector (length 3).}
#' }
#'
#' @examples
#' convert_plane_normal_offset_to_point(c(0, 0, 1), 5)
#' # $point = c(0, 0, 5)
#' # $normal = c(0, 0, 1)
#' @export
convert_plane_normal_offset_to_point <- function(normal, offset) {
  stopifnot(length(normal) == 3L, length(offset) == 1L)
  if (sqrt(sum(normal^2)) == 0) stop("`normal` must be non-zero.")
  n_hat <- normalise(normal)
  point_on_plane <- as.numeric(offset) * n_hat
  list(point = point_on_plane, normal = n_hat)
}


# Points ------------------------------------------------------------------

#' Create a direction vector from two positions
#'
#' Computes the vector from a starting point to an ending point.
#' Optionally returns a unit (normalized) vector.
#'
#' @param start Numeric vector giving the start position.
#' @param end Numeric vector giving the end position.
#' @param unit Logical; if \code{TRUE}, return a unit-length vector.
#'
#' @return A numeric vector representing the direction from \code{start} to \code{end}.
#' @examples
#' create_vector_from_start_end(c(0, 0, 0), c(1, 2, 2))
#' create_vector_from_start_end(c(0, 0, 0), c(1, 2, 2), unit = TRUE)
#'
#' @export
create_vector_from_start_end <- function(start, end, unit = FALSE) {
  if (length(start) != length(end)) stop("Start and end positions need to be the same length")
  vec <- end - start

  if (unit) {
    if (sqrt(sum(vec^2)) == 0) {
      warning("Cannot normalise zero-length vector (start == end). Returning NA")
      return(c(NaN, NaN, NaN))
    }
    vec <- normalise(vec)
  }

  return(vec)
}

# Locate ------------------------------------------------------------------
#' Locate the geometric center (centroid) of 3D points
#'
#' Computes the arithmetic mean position of a set of 3D coordinates,
#' returning the centroid (center point). This is equivalent to taking
#' the average of all x, y, and z values separately.
#'
#' @param ... Either:
#'   \itemize{
#'     \item Three numeric vectors \code{x}, \code{y}, \code{z}, all of equal length, or
#'     \item A single numeric matrix or data frame with columns named \code{x}, \code{y}, and \code{z}.
#'   }
#'
#' @return Numeric vector of length 3 giving the centroid coordinates
#'   \code{c(x_center, y_center, z_center)}.
#'
#' @examples
#' pts <- data.frame(
#'   x = c(0, 2, 4),
#'   y = c(1, 3, 5),
#'   z = c(2, 4, 6)
#' )
#' locate_center(pts)
#' # Returns c(2, 3, 4)
#'
#' locate_center(c(0, 2, 4), c(1, 3, 5), c(2, 4, 6))
#' # Same result
#' @export
locate_center <- function(...) {
  args <- list(...)

  # Case 1: single data frame or matrix
  if (length(args) == 1 && is.data.frame(args[[1]])) {
    df <- args[[1]]
    if (!all(c("x", "y", "z") %in% names(df))) {
      stop("Data frame must have columns named 'x', 'y', and 'z'.")
    }
    x <- df$x
    y <- df$y
    z <- df$z
  } else if (length(args) == 1 && is.matrix(args[[1]])) {
    mat <- args[[1]]
    if (ncol(mat) < 3) stop("Matrix must have at least 3 columns for x, y, z.")
    x <- mat[, 1]
    y <- mat[, 2]
    z <- mat[, 3]
  } else if (length(args) == 3) {
    x <- args[[1]]
    y <- args[[2]]
    z <- args[[3]]
    if (!(is.numeric(x) && is.numeric(y) && is.numeric(z))) {
      stop("x, y, z must all be numeric vectors.")
    }
    if (!(length(x) == length(y) && length(y) == length(z))) {
      stop("x, y, z must have the same length.")
    }
  } else {
    stop("Provide either a data frame/matrix with x, y, z columns, or three numeric vectors.")
  }

  center <- c(mean(x), mean(y), mean(z))
  names(center) <- c("x", "y", "z")
  return(center)
}



# Measure ---------------------------------------------------------------

#' Measure the angle between two vectors
#'
#' Computes the angle between two numeric vectors using the dot product formula:
#' \deqn{\theta = \arccos\left( \frac{a \cdot b}{\|a\|\|b\|} \right)}
#' The angle can be returned in either radians (default) or degrees.
#'
#' @param a,b Numeric vectors of equal length.
#' @param degrees Logical; if \code{TRUE}, return the angle in degrees
#'   instead of radians.
#'
#' @return Numeric scalar giving the angle between \code{a} and \code{b}.
#'   Returns \code{NA} if either vector has zero magnitude.
#'
#' @examples
#' measure_angle_between_vectors(c(1, 0, 0), c(0, 1, 0)) # pi/2 radians
#' measure_angle_between_vectors(c(1, 0, 0), c(0, 1, 0), degrees = TRUE) # 90 degrees
#' measure_angle_between_vectors(c(1, 0, 0), c(1, 0, 0)) # 0
#'
#' @seealso [radians_to_degrees()], [degrees_to_radians()]
#' @export
measure_angle_between_vectors <- function(a, b, degrees = FALSE) {
  na <- magnitude(a)
  nb <- magnitude(b)
  if (na == 0 || nb == 0) {
    return(NA_real_)
  }
  costheta <- as.numeric((a %*% b) / (na * nb))
  costheta <- clamp(costheta, -1, 1)
  theta <- as.numeric(acos(costheta))
  if (degrees) {
    theta <- radians_to_degrees(theta)
  }
  return(theta)
}

#' Compute signed angle between two planes
#'
#' Calculates the signed angle between two planes given their normal vectors.
#' The sign is determined by the right-hand rule about the specified reference
#' axis, allowing differentiation between clockwise and counter-clockwise
#' rotations. Returns an angle in degrees by default.
#'
#' @param n1,n2 Numeric vectors (length 3) giving the normal vectors of the planes.
#' @param ref_axis Numeric vector (length 3); defines the positive rotation direction
#'   according to the right-hand rule.
#' @param degrees Logical; if \code{TRUE}, return the angle in degrees (default).
#'
#' @return A numeric scalar giving the signed angle between the two planes.
#'
#' @details
#' The function computes an oriented angle in the range \eqn{(-180, 180]} using
#' \code{atan2}, with the rotation direction determined by the \code{ref_axis}.
#' This is useful in geometric or molecular contexts (e.g., dihedral angles)
#' where distinguishing rotation sense is important.
#'
#' @examples
#' n1 <- c(0, 0, 1)
#' n2 <- c(0, 1, 0)
#' measure_signed_angle_between_planes(n1, n2, ref_axis = c(1, 0, 0))
#'
#' @export
measure_signed_angle_between_planes <- function(n1, n2, ref_axis, degrees = TRUE) {
  if (sqrt(sum(n1^2)) == 0 || sqrt(sum(n2^2)) == 0) stop("Plane normals must be non-zero.")
  if (sqrt(sum(ref_axis^2)) == 0) stop("Reference axis must be non-zero.")
  n1u <- n1 / sqrt(sum(n1^2))
  n2u <- n2 / sqrt(sum(n2^2))
  axis_u <- ref_axis / sqrt(sum(ref_axis^2))

  m1 <- pracma::cross(axis_u, n1u)
  x <- sum(n1u * n2u) # dot product
  y <- sum(m1 * n2u) # dot product

  angle <- atan2(y, x)
  if (degrees) angle <- radians_to_degrees(angle)
  angle
}

#' Compute angle between two planes
#'
#' Calculates the (unsigned) smallest angle between two planes given their
#' normal vectors. Returns an angle in degrees by default.
#'
#' @param n1,n2 Numeric vectors (length 3) giving the normal vectors of the planes.
#' @param degrees Logical; if \code{TRUE}, return the angle in degrees (default).
#'
#' @return A numeric scalar giving the unsigned angle between the two planes.
#'
#' @details
#' The returned value is always in the range \eqn{[0, 180]} and does not encode
#' rotation direction. For an oriented (signed) version, see
#' \code{\link{measure_signed_angle_between_planes}}.
#'
#' @examples
#' n1 <- c(0, 0, 1)
#' n2 <- c(0, 1, 0)
#' measure_angle_between_planes(n1, n2)
#'
#' @seealso [measure_signed_angle_between_planes()]
#' @export
measure_angle_between_planes <- function(n1, n2, degrees = TRUE) {
  if (sqrt(sum(n1^2)) == 0 || sqrt(sum(n2^2)) == 0) stop("Plane normals must be non-zero.")
  n1u <- n1 / sqrt(sum(n1 * n1))
  n2u <- n2 / sqrt(sum(n2 * n2))
  dp <- sum(n1u * n2u)
  dp <- max(-1, min(1, dp))
  angle <- acos(dp)
  if (degrees) angle <- radians_to_degrees(angle)
  angle
}


#' Measure the distance between two points
#'
#' Computes the Euclidean distance between two points by taking the
#' magnitude of the difference vector (\code{b - a}).
#'
#' @param a Numeric vector giving the first point.
#' @param b Numeric vector giving the second point.
#'
#' @return A numeric scalar giving the Euclidean distance between
#'   \code{a} and \code{b}.
#'
#' @examples
#' measure_distance_between_two_points(c(0, 0, 0), c(3, 4, 0)) # 5
#'
#' @export
measure_distance_between_two_points <- function(a, b) {
  magnitude(b - a)
}

# Apply transformations to tables -----------------------------------------

# Table is either a data.frame or matrix with column names (x,y, z)
# Function expect at least one input: a vector of length 3 3 dimensions: (x, y, z) and return a 3 column data.frame with cols (x, y, z)
# in that order

#' Apply a 3D transformation function to coordinate tables
#'
#' Applies a user-supplied transformation function to each row of a table
#' containing 3D coordinates (`x`, `y`, `z`), returning an updated table
#' with transformed coordinates.
#'
#' @param table A data frame or matrix containing at least the columns
#'   `x`, `y`, and `z`.
#' @param f A function that takes a numeric vector or named list of length 3
#'   (`x`, `y`, `z`) and returns a structure (e.g., named numeric or list)
#'   with elements `x`, `y`, and `z` in that order.
#'
#' @return A data frame with the same columns as `table`, but with the
#'   `x`, `y`, and `z` coordinates replaced by the transformed values.
#'
#' @examples
#' # Example: rotate points around Z-axis
#' rotate_z <- function(p) {
#'   angle <- pi / 4
#'   c(
#'     x = p["x"] * cos(angle) - p["y"] * sin(angle),
#'     y = p["x"] * sin(angle) + p["y"] * cos(angle),
#'     z = p["z"]
#'   )
#' }
#'
#' pts <- data.frame(x = c(1, 0), y = c(0, 1), z = 0)
#' apply_tranformation_to_table(pts, rotate_z)
#'
#' @export
apply_tranformation_to_table <- function(table, f) {
  coords <- table[c("x", "y", "z")]
  res <- apply(X = coords, MARGIN = 1, FUN = f, simplify = FALSE)

  table$x <- vapply(res, FUN = \(d){
    d[1]
  }, FUN.VALUE = numeric(1))
  table$y <- vapply(res, FUN = \(d){
    d[2]
  }, FUN.VALUE = numeric(1))
  table$z <- vapply(res, FUN = \(d){
    d[3]
  }, FUN.VALUE = numeric(1))
  return(table)
}

#' Rotate all points in a table around an axis
#'
#' Rotates each row of a coordinate table (with columns `x`, `y`, `z`) about a
#' specified axis by a given angle. The axis may pass through an arbitrary point
#' (default is the origin).
#'
#' @param table A data frame or matrix containing columns `x`, `y`, `z`.
#' @param rotation_axis Numeric length-3 vector giving the axis direction
#'   (need not be unit; will be normalized internally).
#' @param angle Numeric scalar (radians); rotation angle (right-hand rule).
#' @param point_on_axis Optional numeric length-3 vector for a point on the axis
#'   (default `c(0,0,0)`).
#' @param tol Numeric tolerance for degeneracy checks.
#' @param zap Logical; if `TRUE`, zaps small numerical noise.
#'
#' @return The input `table` with rotated `x`, `y`, `z` columns.
#'
#' @examples
#' pts <- data.frame(x = c(2, 3), y = c(0, 0), z = c(0, 0))
#' axis <- c(0, 0, 1)
#' pivot <- c(1, 0, 0)
#' rotate_table_around_axis(pts, axis, pi / 2, point_on_axis = pivot)
#'
#' @export
rotate_table_around_axis <- function(table,
                                     rotation_axis,
                                     angle,
                                     point_on_axis = c(0, 0, 0),
                                     tol = 1e-8,
                                     zap = TRUE) {
  if (!all(c("x", "y", "z") %in% colnames(table))) {
    stop("`table` must contain columns 'x', 'y', and 'z'.")
  }
  coords <- as.matrix(table[,c("x", "y", "z"), drop=FALSE])
  rotate_one <- function(v) {
    rotate_vector_around_axis_through_point(v,
      rotation_axis = rotation_axis,
      point_on_axis = point_on_axis,
      angle = angle,
      tol = tol,
      zap = zap
    )
  }
  res <- apply(coords, 1, rotate_one, simplify = FALSE)
  res_mx <- as.matrix(do.call("rbind", res))
  table[,"x"] <- res_mx[, 1, drop=TRUE]
  table[,"y"] <- res_mx[, 2, drop=TRUE]
  table[,"z"] <- res_mx[, 3, drop=TRUE]
  table
}


# translate_table_in_direction <- function(table, direction, magnitude, careful = TRUE){
#   if(careful){
#     if (!all(c("x", "y", "z") %in% colnames(table))) {
#       stop("`table` must contain columns 'x', 'y', and 'z'.")
#     }
#   }
#
#   coords <- as.matrix(table[c("x", "y", "z")])
#
#   apply(X = coords, 1, )
# }
