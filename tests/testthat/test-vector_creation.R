test_that("computes raw direction vector (3D)", {
  start <- c(1, 2, 3)
  end   <- c(4,  6, 9)
  expect_equal(
    create_vector_from_start_end(start, end),
    c(3, 4, 6)
  )
})

test_that("computes raw direction vector (2D)", {
  start <- c(-2, 5)
  end   <- c( 3, 1)
  expect_equal(
    create_vector_from_start_end(start, end),
    c(5, -4)
  )
})

test_that("errors when start and end lengths differ", {
  expect_error(
    create_vector_from_start_end(c(0, 1, 2), c(1, 2)),
    "same length"
  )
})

test_that("unit = TRUE returns a unit-length vector with correct direction", {
  start <- c(0, 0, 0)
  end   <- c(3, 1, 5)

  v_raw  <- end - start
  v_unit <- create_vector_from_start_end(start, end, unit = TRUE)

  # Length ≈ 1
  expect_equal(sqrt(sum(v_unit^2)), 1, tolerance = 1e-12)

  # Parallel and same orientation as raw vector:
  # v_unit should equal v_raw / ||v_raw||
  expect_equal(
    v_unit,
    v_raw / sqrt(sum(v_raw^2)),
    tolerance = 1e-12
  )
})

test_that("zero displacement returns zero vector when unit = FALSE", {
  start <- c(2, -3, 7)
  end   <- start
  expect_equal(
    create_vector_from_start_end(start, end, unit = FALSE),
    c(0, 0, 0)
  )
})

test_that("zero displacement with unit = TRUE returns missing data", {
  start <- c(2, -3, 7)
  end   <- start
  expect_equal(
    create_vector_from_start_end(start, end, unit = TRUE), expected = c(NaN, NaN, NaN)
  )
})
