# tests/testthat/test-angle-conversions.R

test_that("radians_to_degrees matches known values", {
  expect_equal(radians_to_degrees(0), 0)
  expect_equal(radians_to_degrees(pi/2), 90, tolerance = 1e-12)
  expect_equal(radians_to_degrees(pi), 180, tolerance = 1e-12)
  expect_equal(radians_to_degrees(2*pi), 360, tolerance = 1e-12)
})

test_that("degrees_to_radians matches known values", {
  expect_equal(degrees_to_radians(0), 0)
  expect_equal(degrees_to_radians(90),  pi/2, tolerance = 1e-12)
  expect_equal(degrees_to_radians(180), pi,   tolerance = 1e-12)
  expect_equal(degrees_to_radians(360), 2*pi, tolerance = 1e-12)
})

test_that("vectorization works for both directions", {
  rads <- c(0, pi/2, pi, 2*pi)
  degs <- c(0, 90, 180, 360)
  expect_equal(radians_to_degrees(rads), degs, tolerance = 1e-12)
  expect_equal(degrees_to_radians(degs), rads, tolerance = 1e-12)
})

test_that("round-trip conversion is identity (within tolerance)", {
  set.seed(123)
  rads <- runif(100, -10*pi, 10*pi)
  degs <- runif(100, -720, 720)

  expect_equal(degrees_to_radians(radians_to_degrees(rads)), rads, tolerance = 1e-12)
  expect_equal(radians_to_degrees(degrees_to_radians(degs)), degs, tolerance = 1e-12)
})

test_that("handles special values (NA, NaN, Inf) by propagation", {
  specials <- c(NA_real_, NaN, Inf, -Inf)
  expect_true(all(is.na(radians_to_degrees(specials)[1])))
  expect_true(is.nan(radians_to_degrees(specials)[2]))
  expect_equal(radians_to_degrees(specials)[3], Inf)
  expect_equal(radians_to_degrees(specials)[4], -Inf)

  expect_true(all(is.na(degrees_to_radians(specials)[1])))
  expect_true(is.nan(degrees_to_radians(specials)[2]))
  expect_equal(degrees_to_radians(specials)[3], Inf)
  expect_equal(degrees_to_radians(specials)[4], -Inf)
})

test_that("zero-length input returns zero-length numeric", {
  expect_equal(length(radians_to_degrees(numeric(0))), 0L)
  expect_equal(length(degrees_to_radians(numeric(0))), 0L)
  expect_true(is.numeric(radians_to_degrees(numeric(0))))
  expect_true(is.numeric(degrees_to_radians(numeric(0))))
})

# tests/testthat/test-wrap-angles.R

test_that("wrap_to_360 handles common cases", {
  expect_equal(wrap_to_360(0), 0)
  expect_equal(wrap_to_360(360), 0)
  expect_equal(wrap_to_360(370), 10)
  expect_equal(wrap_to_360(-30), 330)
  expect_equal(wrap_to_360(-370), 350)

  # Vectorization
  expect_equal(
    wrap_to_360(c(-30, 0, 370, 725)),
    c(330, 0, 10, 5)
  )

  # Range is [0, 360)
  x <- wrap_to_360(seq(-1000, 1000, by = 37))
  expect_true(all(x >= 0 & x < 360))
})

test_that("wrap_to_2pi handles common cases", {
  two_pi <- 2 * pi

  expect_equal(wrap_to_2pi(0), 0)
  expect_equal(wrap_to_2pi(two_pi), 0, tolerance = 1e-12)
  expect_equal(wrap_to_2pi(two_pi + pi/3), pi/3, tolerance = 1e-12)
  expect_equal(wrap_to_2pi(-pi/2), two_pi - pi/2, tolerance = 1e-12)
  expect_equal(wrap_to_2pi(-5*pi), two_pi - pi, tolerance = 1e-12)  # -5π ≡ π

  # Vectorization
  expect_equal(
    wrap_to_2pi(c(-pi/2, 0, 5*pi)),
    c(2*pi - pi/2, 0, pi),
    tolerance = 1e-12
  )

  # Range is [0, 2π)
  x <- wrap_to_2pi(seq(-1000, 1000, length.out = 101) * pi)
  expect_true(all(x >= 0 & x < two_pi + 1e-15))
})

test_that("wrapping is consistent with repeated application (idempotent)", {
  degs <- c(-720, -361, -1, 0, 1, 359, 360, 721)
  rads <- c(-4*pi, -2*pi - 0.1, -1e-9, 0, 1e-9, 2*pi - 1e-9, 2*pi, 7.123)

  expect_equal(wrap_to_360(wrap_to_360(degs)), wrap_to_360(degs))
  expect_equal(wrap_to_2pi(wrap_to_2pi(rads)), wrap_to_2pi(rads), tolerance = 1e-12)
})

test_that("special values propagate as expected", {
  specials <- c(NA_real_, NaN, Inf, -Inf)

  out_deg <- wrap_to_360(specials)
  expect_true(is.na(out_deg[1]))
  expect_true(is.nan(out_deg[2]))
  expect_equal(out_deg[3], NaN)     # Inf %% 360 is NaN in R
  expect_equal(out_deg[4], NaN)     # -Inf %% 360 is NaN in R

  out_rad <- wrap_to_2pi(specials)
  expect_true(is.na(out_rad[1]))
  expect_true(is.nan(out_rad[2]))
  expect_equal(out_rad[3], NaN)     # Inf %% (2π) is NaN
  expect_equal(out_rad[4], NaN)     # -Inf %% (2π) is NaN
})

test_that("zero-length input returns zero-length numeric", {
  expect_equal(length(wrap_to_360(numeric(0))), 0L)
  expect_true(is.numeric(wrap_to_360(numeric(0))))

  expect_equal(length(wrap_to_2pi(numeric(0))), 0L)
  expect_true(is.numeric(wrap_to_2pi(numeric(0))))
})

