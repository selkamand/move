test_that("compute_plane_normal_from_vectors returns unit normal and correct direction up to sign", {
  a <- c(1, 0, 0)
  b <- c(0, 1, 0)
  n <- compute_plane_normal_from_vectors(a, b)
  # unit length
  testthat::expect_equal(sqrt(sum(n^2)), 1, tolerance = 1e-12)
  # parallel to z-axis (either +z or -z)
  testthat::expect_true(all.equal(abs(n), c(0, 0, 1)) == TRUE)
})

test_that("plane conversions are consistent (round-trip)", {
  normal <- c(0, 0, 2)
  point <- c(3, -1, 5)

  x <- convert_plane_point_normal_to_normal_offset(normal, point)
  # x$normal is unit, x$offset is signed distance
  testthat::expect_equal(sqrt(sum(x$normal^2)), 1, tolerance = 1e-12)
  # Convert back
  y <- convert_plane_normal_offset_to_point(x$normal, x$offset)
  # y$point should be closest point to origin on plane (on the normal line)
  testthat::expect_equal(sum(x$normal * y$point), x$offset, tolerance = 1e-12)

  # Check original point lies on plane: n·p == offset
  testthat::expect_equal(sum(x$normal * point), x$offset, tolerance = 1e-12)
})

test_that("locate_center computes centroid from data.frame, matrix, or vectors", {
  df <- data.frame(x = c(0, 2, 4), y = c(1, 3, 5), z = c(2, 4, 6))
  mat <- as.matrix(df)
  cx <- locate_center(df)
  cy <- locate_center(mat)
  cz <- locate_center(df$x, df$y, df$z)
  testthat::expect_equal(cx, c(x = 2, y = 3, z = 4))
  testthat::expect_equal(cy, c(x = 2, y = 3, z = 4))
  testthat::expect_equal(cz, c(x = 2, y = 3, z = 4))
})
