make_pts <- function(xs, ys, zs, extra = NULL) {
  df <- data.frame(x = xs, y = ys, z = zs)
  if (!is.null(extra)) df$label <- extra
  df
}

expect_table_xyz_equal <- function(tbl, xyz_mat, tol = 1e-8) {
  testthat::expect_equal(as.numeric(tbl$x), as.numeric(xyz_mat[, 1]), tolerance = tol)
  testthat::expect_equal(as.numeric(tbl$y), as.numeric(xyz_mat[, 2]), tolerance = tol)
  testthat::expect_equal(as.numeric(tbl$z), as.numeric(xyz_mat[, 3]), tolerance = tol)
}

test_that("rotates around origin matches per-row rotate_vector_around_axis", {
  pts <- make_pts(c(1, 0), c(0, 1), c(0, 0))
  axis <- c(0, 0, 1)
  angle <- pi / 2
  rotated_tbl <- rotate_table_around_axis(pts, axis, angle)

  # Compute expected row-wise
  coords <- as.matrix(pts[, c("x", "y", "z")])
  exp_rows <- t(apply(coords, 1, function(v) rotate_vector_around_axis(v, axis, angle)))

  expect_table_xyz_equal(rotated_tbl, exp_rows)
})

test_that("rotates around arbitrary point matches translate-rotate-translate", {
  pts <- make_pts(c(2, 3), c(0, 0), c(0, 0))
  axis <- c(0, 0, 1)
  pivot <- c(1, 0, 0)
  angle <- pi / 2
  rotated_tbl <- rotate_table_around_axis(pts, axis, angle, point_on_axis = pivot)

  # Manual expected: shift -> rotate -> shift back
  coords <- as.matrix(pts[, c("x", "y", "z")])
  exp_rows <- t(apply(coords, 1, function(v) {
    v_c <- v - pivot
    v_r <- rotate_vector_around_axis(v_c, axis, angle)
    v_r + pivot
  }))
  expect_table_xyz_equal(rotated_tbl, exp_rows)
})

test_that("preserves non-coordinate columns", {
  pts <- make_pts(c(2, 3), c(0, 0), c(0, 0), extra = c("a", "b"))
  out <- rotate_table_around_axis(pts, c(0, 0, 1), pi / 2)
  testthat::expect_true("label" %in% names(out))
  testthat::expect_equal(out$label, pts$label)
})

test_that("validation: missing columns and zero axis", {
  pts <- data.frame(a = 1, b = 2, c = 3)
  testthat::expect_error(rotate_table_around_axis(pts, c(0, 0, 1), pi / 2))
  pts2 <- make_pts(1, 0, 0)
  testthat::expect_error(rotate_table_around_axis(pts2, c(0, 0, 0), pi / 2))
})

test_that("zap removes small numerical noise", {
  pts <- make_pts(1, 0, 0)
  axis <- c(0, 0, 1)
  # Rotate by 2*pi should return exactly same when zapped
  out <- rotate_table_around_axis(pts, axis, 2 * pi, zap = TRUE)
  testthat::expect_equal(out$x, 1)
  testthat::expect_equal(out$y, 0)
  testthat::expect_equal(out$z, 0)
})
