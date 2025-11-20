
# Helpers -----------------------------------------------------------------

make_pts <- function(xs, ys, zs, extra = NULL) {
  df <- data.frame(x = xs, y = ys, z = zs)
  if (!is.null(extra)) df$label <- extra
  df
}

expect_table_xyz_equal <- function(tbl, xyz_mat, tol = 1e-8) {
  testthat::expect_equal(as.numeric(tbl[, "x", drop=TRUE]), as.numeric(xyz_mat[, 1]), tolerance = tol)
  testthat::expect_equal(as.numeric(tbl[, "y", drop=TRUE]), as.numeric(xyz_mat[, 2]), tolerance = tol)
  testthat::expect_equal(as.numeric(tbl[, "z", drop=TRUE]), as.numeric(xyz_mat[, 3]), tolerance = tol)
}


# apply_tranformation_to_table --------------------------------------------


test_that("apply_tranformation_to_table works with a matrix input", {
  # Base points as data.frame
  pts_df <- make_pts(
    xs = c(1, 2),
    ys = c(3, 4),
    zs = c(5, 6)
  )

  # Same coordinates as a matrix with x,y,z colnames
  pts_mat <- as.matrix(pts_df[, c("x", "y", "z")])
  colnames(pts_mat) <- c("x", "y", "z")

  # Simple affine transform: translate by (1, 2, 3)
  shift_fun <- function(p) {
    c(
      x = p["x"] + 1,
      y = p["y"] + 2,
      z = p["z"] + 3
    )
  }

  out_df  <- apply_tranformation_to_table(pts_df,  shift_fun)
  out_mat <- apply_tranformation_to_table(pts_mat, shift_fun)

  # Both should be data.frames with same xyz numerics
  testthat::expect_s3_class(out_df,  "data.frame")
  testthat::expect_true(is.matrix(out_mat))
  expect_table_xyz_equal(
    out_mat,
    as.matrix(out_df[, c("x", "y", "z")])
  )
})

# rotate_table_around_axis -----------------------------------------------

test_that("rotate_table_around_axis works with matrix input", {
  pts_df <- make_pts(c(1, 0), c(0, 1), c(0, 0))
  pts_mat <- as.matrix(pts_df[, c("x", "y", "z")])
  colnames(pts_mat) <- c("x", "y", "z")

  axis  <- c(0, 0, 1)
  angle <- pi / 2

  out_df  <- rotate_table_around_axis(pts_df,  axis, angle)
  out_mat <- rotate_table_around_axis(pts_mat, axis, angle)

  expect_table_xyz_equal(
    out_mat,
    as.matrix(out_df[, c("x", "y", "z")])
  )
})

test_that("rotate_table_around_axis validation works for matrices", {
  # Matrix missing required columns (no x,y,z names)
  bad_mat <- matrix(1:9, ncol = 3)
  colnames(bad_mat) <- c("a", "b", "c")

  testthat::expect_error(
    rotate_table_around_axis(bad_mat, c(0, 0, 1), pi / 2),
    "must contain columns 'x', 'y', and 'z'"
  )
})

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

# translate_table_in_direction --------------------------------------------

test_that("translate_table_in_direction translates data.frames correctly", {
  pts <- make_pts(
    xs = c(0, 1),
    ys = c(0, 0),
    zs = c(0, 0)
  )

  direction <- c(1, 0, 0)
  magnitude <- 2

  out <- translate_table_in_direction(pts, direction, magnitude)

  # Expected: x shifted by +2 along (1,0,0)
  expected <- make_pts(
    xs = c(2, 3),
    ys = c(0, 0),
    zs = c(0, 0)
  )

  expect_table_xyz_equal(
    out,
    as.matrix(expected[, c("x", "y", "z")])
  )
})

test_that("translate_table_in_direction preserves non-coordinate columns", {
  pts <- make_pts(
    xs    = c(0, 1),
    ys    = c(0, 0),
    zs    = c(0, 0),
    extra = c("a", "b")
  )

  direction <- c(0, 0, 1)
  magnitude <- 5

  out <- translate_table_in_direction(pts, direction, magnitude)

  testthat::expect_true("label" %in% names(out))
  testthat::expect_equal(out$label, pts$label)
})

test_that("translate_table_in_direction works with matrix input", {
  pts_df <- make_pts(
    xs = c(0, 1),
    ys = c(0, 0),
    zs = c(0, 0)
  )
  pts_mat <- as.matrix(pts_df[, c("x", "y", "z")])
  colnames(pts_mat) <- c("x", "y", "z")

  direction <- c(0, 1, 0)
  magnitude <- 3

  out_df  <- translate_table_in_direction(pts_df,  direction, magnitude)
  out_mat <- translate_table_in_direction(pts_mat, direction, magnitude)

  expect_table_xyz_equal(
    out_mat,
    as.matrix(out_df[, c("x", "y", "z")])
  )
})

test_that("translate_table_in_direction validation works for matrices", {
  bad_mat <- matrix(1:9, ncol = 3)
  colnames(bad_mat) <- c("a", "b", "c")

  testthat::expect_error(
    translate_table_in_direction(bad_mat, c(1, 0, 0), magnitude = 1),
    "must contain columns 'x', 'y', and 'z'"
  )
})

# Rowname & Column Preservation -------------------------------------------------

## apply_tranformation_to_table ----------------------------------------------

test_that("apply_tranformation_to_table preserves rownames & extra cols (data.frame)", {
  pts <- make_pts(
    xs    = c(1, 2),
    ys    = c(3, 4),
    zs    = c(5, 6),
    extra = c("a", "b")
  )
  rownames(pts) <- c("row1", "row2")

  shift_fun <- function(p) {
    c(x = p["x"] + 1,
      y = p["y"] + 2,
      z = p["z"] + 3)
  }

  out <- apply_tranformation_to_table(pts, shift_fun)

  testthat::expect_identical(rownames(out), rownames(pts))
  testthat::expect_true("label" %in% names(out))
  testthat::expect_identical(out$label, pts$label)
})

test_that("apply_tranformation_to_table preserves rownames & extra cols (matrix)", {
  pts_mat <- cbind(
    x     = c(1, 2),
    y     = c(3, 4),
    z     = c(5, 6),
    extra = c(10, 20)
  )
  rownames(pts_mat) <- c("row1", "row2")

  shift_fun <- function(p) {
    c(x = p["x"] + 1,
      y = p["y"] + 2,
      z = p["z"] + 3)
  }

  out <- apply_tranformation_to_table(pts_mat, shift_fun)

  testthat::expect_identical(rownames(out), rownames(pts_mat))
  testthat::expect_true("extra" %in% colnames(out))
  testthat::expect_equal(out[, "extra"], pts_mat[, "extra"])
})

## rotate_table_around_axis -------------------------------------------------

test_that("rotate_table_around_axis preserves rownames & extra cols (data.frame)", {
  pts <- make_pts(
    xs    = c(2, 3),
    ys    = c(0, 0),
    zs    = c(0, 0),
    extra = c("a", "b")
  )
  rownames(pts) <- c("row1", "row2")

  out <- rotate_table_around_axis(pts, c(0, 0, 1), angle = pi / 2)

  testthat::expect_identical(rownames(out), rownames(pts))
  testthat::expect_true("label" %in% names(out))
  testthat::expect_identical(out$label, pts$label)
})

test_that("rotate_table_around_axis preserves rownames & extra cols (matrix)", {
  pts_mat <- cbind(
    x     = c(2, 3),
    y     = c(0, 0),
    z     = c(0, 0),
    extra = c(10, 20)
  )
  rownames(pts_mat) <- c("row1", "row2")

  out <- rotate_table_around_axis(pts_mat, c(0, 0, 1), angle = pi / 2)

  testthat::expect_identical(rownames(out), rownames(pts_mat))
  testthat::expect_true("extra" %in% colnames(out))
  testthat::expect_equal(out[, "extra"], pts_mat[, "extra"])
})

## translate_table_in_direction ---------------------------------------------

test_that("translate_table_in_direction preserves rownames & extra cols (data.frame)", {
  pts <- make_pts(
    xs    = c(0, 1),
    ys    = c(0, 0),
    zs    = c(0, 0),
    extra = c("a", "b")
  )
  rownames(pts) <- c("row1", "row2")

  out <- translate_table_in_direction(pts, direction = c(1, 0, 0), magnitude = 2)

  testthat::expect_identical(rownames(out), rownames(pts))
  testthat::expect_true("label" %in% names(out))
  testthat::expect_identical(out$label, pts$label)
})

test_that("translate_table_in_direction preserves rownames & extra cols (matrix)", {
  pts_mat <- cbind(
    x     = c(0, 1),
    y     = c(0, 0),
    z     = c(0, 0),
    extra = c(10, 20)
  )
  rownames(pts_mat) <- c("row1", "row2")

  out <- translate_table_in_direction(pts_mat, direction = c(0, 1, 0), magnitude = 3)

  testthat::expect_identical(rownames(out), rownames(pts_mat))
  testthat::expect_true("extra" %in% colnames(out))
  testthat::expect_equal(out[, "extra"], pts_mat[, "extra"])
})



