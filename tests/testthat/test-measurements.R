test_that("perpendicular vectors return pi/2 (radians) and 90 (degrees)", {
  a <- c(1, 0, 0)
  b <- c(0, 1, 0)

  theta_rad <- measure_angle_between_vectors(a, b)
  theta_deg <- measure_angle_between_vectors(a, b, degrees = TRUE)

  expect_equal(theta_rad, pi/2, tolerance = 1e-12)
  expect_equal(theta_deg, 90, tolerance = 1e-10)
})

test_that("parallel vectors return 0", {
  a <- c(1, 2, 3)
  b <- 5 * a # same direction

  theta <- measure_angle_between_vectors(a, b)
  expect_equal(theta, 0, tolerance = 1e-12)
})

test_that("anti-parallel vectors return pi (or 180 degrees)", {
  a <- c(1, 2, 3)
  b <- -a

  theta_rad <- measure_angle_between_vectors(a, b)
  theta_deg <- measure_angle_between_vectors(a, b, degrees = TRUE)

  expect_equal(theta_rad, pi, tolerance = 1e-12)
  expect_equal(theta_deg, 180, tolerance = 1e-10)
})

test_that("angle is symmetric and scale-invariant", {
  a <- c(2, -3, 4)
  b <- c(-1, 5, 2)

  # symmetry: angle(a,b) == angle(b,a)
  expect_equal(
    measure_angle_between_vectors(a, b),
    measure_angle_between_vectors(b, a),
    tolerance = 1e-12
  )

  # scale-invariant: scaling either vector doesn't change the angle (unless vector direction gets flipped)
  expect_equal(
    measure_angle_between_vectors(3 * a, b),
    measure_angle_between_vectors(a, 2 * b), # scale & flip b (flip should complement to pi)
    tolerance = 1e-12
  )
})

test_that("degrees flag matches radians_to_degrees of radians result", {
  a <- c(1.2, -0.7, 3.4)
  b <- c(-2.1, 0.5, 0.9)

  rad <- measure_angle_between_vectors(a, b)
  deg <- measure_angle_between_vectors(a, b, degrees = TRUE)

  expect_equal(deg, radians_to_degrees(rad), tolerance = 1e-10)
})

test_that("works in 2D as well as 3D (any equal-length vectors)", {
  a2 <- c(1, 0)
  b2 <- c(0, 1)
  expect_equal(measure_angle_between_vectors(a2, b2), pi/2, tolerance = 1e-12)

  a4 <- c(1, 0, 0, 0)
  b4 <- c(0, 1, 0, 0)
  expect_equal(measure_angle_between_vectors(a4, b4), pi/2, tolerance = 1e-12)
})

test_that("result is within [0, pi] radians (or [0, 180] degrees)", {
  set.seed(123)
  for (i in 1:20) {
    a <- rnorm(3)
    b <- rnorm(3)
    th <- measure_angle_between_vectors(a, b)
    th_deg <- measure_angle_between_vectors(a, b, degrees = TRUE)
    expect_true(is.finite(th))
    expect_true(th >= 0 - 1e-12 && th <= pi + 1e-12)
    expect_true(th_deg >= 0 - 1e-10 && th_deg <= 180 + 1e-10)
  }
})

test_that("zero vector input yields NA/NaN angle (division by zero path)", {
  # Behavior depends on your magnitude() helper; with zero magnitude,
  # the dot/denominator is undefined, so acos(NaN) -> NaN (is.na() is TRUE for NaN).
  a <- c(0, 0, 0)
  b <- c(1, 2, 3)

  th1 <- suppressWarnings(measure_angle_between_vectors(a, b))
  th2 <- suppressWarnings(measure_angle_between_vectors(b, a))

  expect_true(is.na(th1))
  expect_true(is.na(th2))
})

test_that("numerical stability for nearly parallel vectors", {
  # Construct nearly parallel unit vectors
  a <- c(1, 1e-12, 0); a <- a / sqrt(sum(a^2))
  b <- c(1, 0,      0); b <- b / sqrt(sum(b^2))

  th <- measure_angle_between_vectors(a, b)
  expect_gte(th, 0)
  expect_lte(th, 1e-6)  # should be extremely small
})

# tests/testthat/test-measure_distance_between_two_points.R

test_that("computes Euclidean distance in 2D", {
  a <- c(0, 0)
  b <- c(3, 4)
  expect_equal(measure_distance_between_two_points(a, b), 5)
})

test_that("computes Euclidean distance in 3D", {
  a <- c(1, 2, 0)
  b <- c(4, 6, 3)
  # sqrt((3)^2 + (4)^2 + (3)^2) = sqrt(34)
  expect_equal(
    measure_distance_between_two_points(a, b),
    sqrt(34),
    tolerance = 1e-12
  )
})

test_that("is symmetric: d(a,b) == d(b,a)", {
  a <- c(-2, 7, 1)
  b <- c(5, -1, 4)
  dab <- measure_distance_between_two_points(a, b)
  dba <- measure_distance_between_two_points(b, a)
  expect_equal(dab, dba, tolerance = 1e-12)
})

test_that("zero distance for identical points", {
  p <- c(2.5, -3, 10)
  expect_equal(measure_distance_between_two_points(p, p), 0)
})

test_that("handles higher-dimensional vectors", {
  set.seed(1)
  a <- rnorm(10)
  b <- rnorm(10)
  # Reference using explicit formula
  ref <- sqrt(sum((b - a)^2))
  expect_equal(
    measure_distance_between_two_points(a, b),
    ref,
    tolerance = 1e-12
  )
})

test_that("distance is non-negative", {
  a <- c(-1, -2, -3)
  b <- c(4, 0, 5)
  d <- measure_distance_between_two_points(a, b)
  expect_gte(d, 0)
})

