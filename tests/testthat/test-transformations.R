
test_that("rotate_vector_to_align_with_target handles antiparallel vectors (180 degrees)", {
  v <- c(1, 2, 3)
  target <- -v
  v_rot <- rotate_vector_to_align_with_target(v, target)
  # Should align direction with target and preserve magnitude
  expect_true(all(is.finite(v_rot)))
  expect_equal(as.numeric(sum(v_rot * target)) / (sqrt(sum(v_rot^2)) * sqrt(sum(target^2))), 1, tolerance = 1e-7)
  expect_equal(sqrt(sum(v_rot^2)), sqrt(sum(v^2)), tolerance = 1e-7)
})

test_that("rotate_vector_to_align_with_target returns axis_plus_angle consistently and robustly", {
  v <- c(0, 0, 1)
  target <- c(0, 0, 1)
  params <- rotate_vector_to_align_with_target(v, target, return = "axis_plus_angle")
  expect_true(is.list(params))
  expect_named(params, c("axis", "angle"))
  expect_true(is.numeric(params$axis))
  expect_length(params$axis, 3)
  expect_true(is.finite(params$angle))
  expect_equal(params$angle, 0, tolerance = 1e-8)
})

test_that("rotate_vector_to_align_with_target errors on zero target", {
  v <- c(1, 0, 0)
  target <- c(0, 0, 0)

  expect_warning(rotate_vector_to_align_with_target(c(1, 0, 0), target = c(0, 0, 0)))
})

test_that("rotate_vector_around_axis errors for zero rotation axis", {
  v <- c(1, 0, 0)
  axis <- c(0, 0, 0)
  expect_error(rotate_vector_around_axis(v, axis, pi / 4))
})

test_that("rotate_vector_around_axis treats very small angles as no-op (tolerance)", {
  v <- c(1, 0, 0)
  axis <- c(0, 0, 1)
  out <- rotate_vector_around_axis(v, axis, 1e-12)
  expect_equal(out, v, tolerance = 1e-10)
})

test_that("rotate_vector_into_a_plane validates normal and returns sane axis_plus_angle for 0-angle", {
  v <- c(1, 0, 0)
  # v is already in plane z=0 (normal (0,0,1))
  params <- rotate_vector_into_a_plane(v, c(0, 0, 1), return = "axis_plus_angle")
  expect_true(is.list(params))
  expect_named(params, c("axis", "angle"))
  expect_true(all(is.finite(params$axis)))
  expect_equal(params$angle, 0, tolerance = 1e-8)
  # zero normal should error
  expect_error(rotate_vector_into_a_plane(v, c(0, 0, 0)))
})

test_that("measure_angle_between_vectors clamps and handles zero vectors", {
  a <- c(1, 0, 0)
  b <- c(1, 0, 0)
  expect_equal(measure_angle_between_vectors(a, b), 0)
  expect_true(is.na(measure_angle_between_vectors(a, c(0, 0, 0))))
  # near floating error: cos slightly > 1 should clamp
  a2 <- c(1, 0, 0)
  b2 <- a2 + 1e-14
  expect_true(is.finite(measure_angle_between_vectors(a2, b2)))
})

test_that("measure_angle_between_planes validates inputs and clamps dot", {
  n1 <- c(0, 0, 1)
  n2 <- c(0, 1, 0)
  expect_equal(measure_angle_between_planes(n1, n2), 90)
  expect_error(measure_angle_between_planes(c(0, 0, 0), n2))
  expect_error(measure_angle_between_planes(n1, c(0, 0, 0)))
})

test_that("measure_signed_angle_between_planes validates reference axis and normals", {
  n1 <- c(0, 0, 1)
  n2 <- c(0, 1, 0)
  ref <- c(1, 0, 0)
  ang <- measure_signed_angle_between_planes(n1, n2, ref_axis = ref)
  expect_true(is.finite(ang))
  expect_error(measure_signed_angle_between_planes(c(0, 0, 0), n2, ref_axis = ref))
  expect_error(measure_signed_angle_between_planes(n1, n2, ref_axis = c(0, 0, 0)))
})

test_that("translate_position_in_direction errors on zero direction", {
  expect_error(translate_position_in_direction(c(0, 0, 0), c(0, 0, 0), 1))
})

test_that("create_vector_from_start_end with unit=TRUE handles zero displacement", {
  start <- c(1, 2, 3)
  end <- c(1, 2, 3)
  expect_warning(create_vector_from_start_end(start, end, unit = TRUE))
  # unit=FALSE should return zeros
  expect_equal(create_vector_from_start_end(start, end, unit = FALSE), c(0, 0, 0))
})

test_that("plane conversion functions and plane normal validate zero/parallel inputs", {
  expect_error(convert_plane_point_normal_to_normal_offset(c(0, 0, 0), c(1, 2, 3)))
  expect_error(convert_plane_normal_offset_to_point(c(0, 0, 0), 5))
  expect_error(compute_plane_normal_from_vectors(c(1, 0, 0), c(2, 0, 0))) # parallel
  expect_error(compute_plane_normal_from_vectors(c(0, 0, 0), c(0, 1, 0))) # zero
})
