test_that("project_vector_into_vector computes correct projection", {
  a <- c(2, 3, 4)
  b <- c(1, 0, 0)
  expect_equal(project_vector_into_vector(a, b), c(2, 0, 0))

  # Non-axis-aligned
  a <- c(1, 2)
  b <- c(1, 1)
  proj <- project_vector_into_vector(a, b)
  # Expected: ((a·b)/(b·b)) * b
  exp <- sum(a * b) / sum(b * b) * b
  expect_equal(proj, exp, tolerance = 1e-12)
})

test_that("compute_scalar_projection returns signed magnitude", {
  a <- c(2, 3, 4)
  b <- c(1, 0, 0)
  expect_equal(compute_scalar_projection(a, b), 2)

  # Opposite direction yields negative
  a2 <- c(-2, 2)
  b2 <- c(2, 0)
  expect_equal(compute_scalar_projection(c(-2, 2), c(2, 0)), expected = -2)
})

test_that("project_vector_into_plane removes normal component", {
  v <- c(1, 2, 3)
  n <- c(0, 0, 1)
  in_plane <- project_vector_into_plane(v, n)
  expect_equal(in_plane, c(1, 2, 0), tolerance = 1e-12)
  # Result is orthogonal to the normal
  expect_equal(sum(in_plane * n), 0, tolerance = 1e-12)

  # Non-unit normal
  n2 <- c(0, 2, -2)
  in_plane2 <- project_vector_into_plane(v, n2)
  expect_equal(sum(in_plane2 * n2), 0, tolerance = 1e-12)
})

test_that("translate_position_by_vector shifts correctly", {
  pos <- c(1, 2, 3)
  vec <- c(-1, 0.5, 2)
  expect_equal(translate_position_by_vector(pos, vec), pos + vec)
})

test_that("compute_translation_vector computes target - position", {
  pos <- c(1, 2, 3)
  tgt <- c(4, 6, 3)
  expect_equal(compute_translation_vector(pos, tgt), c(3, 4, 0))
})

test_that("translate_position_in_direction moves along normalized direction", {
  pos <- c(0, 0, 0)
  dir <- c(3, 4, 0) # norm 5
  out <- translate_position_in_direction(pos, dir, 10)
  # Expected translation is unit(dir) * 10 = (3/5,4/5,0)*10
  expect_equal(out, c(6, 8, 0), tolerance = 1e-12)
})
