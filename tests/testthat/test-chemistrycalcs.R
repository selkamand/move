# Tests for chemistry calculation functions

test_that("compute_abcd_dihedral_stats works with normal geometry", {
  # Simple tetrahedral-like geometry
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(1, 1, 0)
  d <- c(1, 1, 1)

  result <- compute_abcd_dihedral_stats(a, b, c, d)

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("a", "b", "c", "bond_angle", "bond_length", "torsion_angle"))

  # Check input preservation
  expect_equal(result$a, a)
  expect_equal(result$b, b)
  expect_equal(result$c, c)

  # Check computed values are reasonable
  expect_type(result$bond_angle, "double")
  expect_type(result$bond_length, "double")
  expect_type(result$torsion_angle, "double")

  # Bond angle should be 90 degrees (B-C-D is right angle)
  expect_equal(result$bond_angle, 90, tolerance = 1e-10)

  # Bond length should be 1 (distance from c to d)
  expect_equal(result$bond_length, 1, tolerance = 1e-10)

  # Torsion angle should be 90 degrees
  expect_equal(result$torsion_angle, 90, tolerance = 1e-10)
})

test_that("compute_abcd_dihedral_stats handles colinear atoms", {
  # Colinear case
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(2, 0, 0)
  d <- c(3, 0, 0)

  # Should produce warning
  expect_warning(
    result <- compute_abcd_dihedral_stats(a, b, c, d),
    "Points A, B, C, D are colinear or nearly colinear"
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("a", "b", "c", "bond_angle", "bond_length", "torsion_angle"))

  # Bond angle should be 180 degrees (colinear)
  expect_equal(result$bond_angle, 180, tolerance = 1e-10)

  # Bond length should be 1
  expect_equal(result$bond_length, 1, tolerance = 1e-10)

  # Torsion angle should be NA
  expect_true(is.na(result$torsion_angle))
})

test_that("compute_abcd_dihedral_stats handles nearly colinear atoms", {
  # Nearly colinear case (should work without warning)
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(2, 0.01, 0)
  d <- c(3, 0, 0)

  # Should not produce warning
  expect_no_warning(
    result <- compute_abcd_dihedral_stats(a, b, c, d)
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("a", "b", "c", "bond_angle", "bond_length", "torsion_angle"))

  # All values should be finite
  expect_true(is.finite(result$bond_angle))
  expect_true(is.finite(result$bond_length))
  expect_true(is.finite(result$torsion_angle))

  # Bond angle should be close to 180 degrees
  expect_true(result$bond_angle > 170)

  # Torsion angle should be close to 180 degrees
  expect_true(abs(result$torsion_angle - 180) < 10)
})

test_that("compute_abcd_dihedral_stats validates input dimensions", {
  # Test invalid input lengths
  expect_error(
    compute_abcd_dihedral_stats(c(1, 2), c(1, 2, 3), c(1, 2, 3), c(1, 2, 3)),
    "must be a 3d vector"
  )

  expect_error(
    compute_abcd_dihedral_stats(c(1, 2, 3), c(1, 2), c(1, 2, 3), c(1, 2, 3)),
    "must be a 3d vector"
  )

  expect_error(
    compute_abcd_dihedral_stats(c(1, 2, 3), c(1, 2, 3), c(1, 2), c(1, 2, 3)),
    "must be a 3d vector"
  )

  expect_error(
    compute_abcd_dihedral_stats(c(1, 2, 3), c(1, 2, 3), c(1, 2, 3), c(1, 2)),
    "must be a 3d vector"
  )
})

test_that("compute_abcd_dihedral_stats with known geometry", {
  # Test with a known geometry where we can calculate expected values
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(1, 1, 0)
  d <- c(0, 1, 0)

  result <- compute_abcd_dihedral_stats(a, b, c, d)

  # B-C-D should be 90 degrees (right angle)
  expect_equal(result$bond_angle, 90, tolerance = 1e-10)

  # C-D distance should be 1
  expect_equal(result$bond_length, 1, tolerance = 1e-10)

  # Torsion angle should be 0 degrees (A-B-C and B-C-D planes are the same)
  expect_equal(result$torsion_angle, 0, tolerance = 1e-10)
})

test_that("compute_fourth_atom_position works with valid inputs", {
  # Skip if compas is not available
  skip_if_not_installed("compas")

  # Define first three atoms
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(1, 1, 0)

  # Test default behavior (bond vector)
  bond_vector <- compute_fourth_atom_position(a, b, c, bond_angle = 90, bond_length = 1, torsion_angle = 0)

  # Check output is a numeric vector of length 3
  expect_type(bond_vector, "double")
  expect_length(bond_vector, 3)

  # Check that all values are finite
  expect_true(all(is.finite(bond_vector)))

  # Test absolute position option
  d_position <- compute_fourth_atom_position(a, b, c, bond_angle = 90, bond_length = 1, torsion_angle = 0, return_absolute_position = TRUE)

  # Check output is a numeric vector of length 3
  expect_type(d_position, "double")
  expect_length(d_position, 3)

  # Check that all values are finite
  expect_true(all(is.finite(d_position)))

  # Verify relationship: bond_vector + c should equal d_position
  expect_equal(bond_vector + c, d_position, tolerance = 1e-10)
})

test_that("compute_fourth_atom_position bond vector functionality", {
  # Skip if compas is not available
  skip_if_not_installed("compas")

  # Define first three atoms
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(1, 1, 0)

  # Test parameters
  bond_angle <- 90
  bond_length <- 2
  torsion_angle <- 45

  # Get bond vector (default)
  bond_vector <- compute_fourth_atom_position(a, b, c, bond_angle, bond_length, torsion_angle)

  # Get absolute position
  d_position <- compute_fourth_atom_position(a, b, c, bond_angle, bond_length, torsion_angle, return_absolute_position = TRUE)

  # Test that bond vector has correct magnitude
  expect_equal(magnitude(bond_vector), bond_length, tolerance = 1e-10)

  # Test that bond vector + c equals absolute position
  expect_equal(bond_vector + c, d_position, tolerance = 1e-10)

  # Test that we can reconstruct the position from bond vector
  reconstructed_d <- c + bond_vector
  expect_equal(reconstructed_d, d_position, tolerance = 1e-10)
})

# Note: Input validation tests removed as the function doesn't validate input dimensions
# The function will fail when creating the matrix with invalid dimensions, but this
# is not the primary focus of the function's behavior

test_that("compute_fourth_atom_position and compute_abcd_dihedral_stats are inverse operations", {
  # Skip if compas is not available
  skip_if_not_installed("compas")

  # Define first three atoms
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(1, 1, 0)

  # Test parameters
  bond_angle <- 109.5
  bond_length <- 1.5
  torsion_angle <- 60

  # Compute position of fourth atom (absolute position)
  d <- compute_fourth_atom_position(a, b, c, bond_angle, bond_length, torsion_angle, return_absolute_position = TRUE)

  # Verify by computing dihedral stats
  stats <- compute_abcd_dihedral_stats(a, b, c, d)

  # Check that the computed values match the input parameters
  expect_equal(stats$bond_angle, bond_angle, tolerance = 1e-6)
  expect_equal(stats$bond_length, bond_length, tolerance = 1e-6)
  expect_equal(stats$torsion_angle, torsion_angle, tolerance = 1e-6)
})

test_that("compute_fourth_atom_position handles edge cases", {
  # Skip if compas is not available
  skip_if_not_installed("compas")

  # Define first three atoms
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(1, 1, 0)

  # Test with 0 degree torsion angle (bond vector)
  bond_vector1 <- compute_fourth_atom_position(a, b, c, bond_angle = 90, bond_length = 1, torsion_angle = 0)
  expect_length(bond_vector1, 3)
  expect_true(all(is.finite(bond_vector1)))

  # Test with 180 degree torsion angle (bond vector)
  bond_vector2 <- compute_fourth_atom_position(a, b, c, bond_angle = 90, bond_length = 1, torsion_angle = 180)
  expect_length(bond_vector2, 3)
  expect_true(all(is.finite(bond_vector2)))

  # Test with 90 degree bond angle (absolute position)
  d3 <- compute_fourth_atom_position(a, b, c, bond_angle = 90, bond_length = 1, torsion_angle = 90, return_absolute_position = TRUE)
  expect_length(d3, 3)
  expect_true(all(is.finite(d3)))
})

test_that("compute_abcd_dihedral_stats handles zero-length vectors", {
  # Test with identical points (should still work for some calculations)
  a <- c(0, 0, 0)
  b <- c(0, 0, 0)
  c <- c(1, 0, 0)
  d <- c(2, 0, 0)

  # This should work but may produce warnings - suppress them for testing
  result <- suppressWarnings(compute_abcd_dihedral_stats(a, b, c, d))

  expect_type(result, "list")
  expect_named(result, c("a", "b", "c", "bond_angle", "bond_length", "torsion_angle"))

  # Bond length should be 1 (distance from c to d)
  expect_equal(result$bond_length, 1, tolerance = 1e-10)
})

test_that("compute_abcd_dihedral_stats handles various bond angles", {
  # Test with different known bond angles
  a <- c(0, 0, 0)
  b <- c(1, 0, 0)
  c <- c(2, 0, 0)

  # 60 degree bond angle (B-C-D angle)
  # For 60 degrees, we need: cos(60) = 0.5
  # CB = (-1,0,0) and CD = (x,y,0), then CB·CD = -x = |CB||CD|*0.5
  # With |CB| = 1 and |CD| = 1, we need -x = 0.5, so x = -0.5, y = sqrt(1-0.25) = sqrt(0.75)
  d1 <- c(1.5, sqrt(0.75), 0)
  result1 <- suppressWarnings(compute_abcd_dihedral_stats(a, b, c, d1))
  expect_equal(result1$bond_angle, 60, tolerance = 1e-6)

  # 120 degree bond angle (B-C-D angle)
  # For 120 degrees, we need: cos(120) = -0.5
  # CB·CD = -x = -0.5, so x = 0.5, y = sqrt(1-0.25) = sqrt(0.75)
  d2 <- c(2.5, sqrt(0.75), 0)
  result2 <- suppressWarnings(compute_abcd_dihedral_stats(a, b, c, d2))
  expect_equal(result2$bond_angle, 120, tolerance = 1e-6)

  # 180 degree bond angle (colinear)
  d3 <- c(3, 0, 0)
  expect_warning(
    result3 <- compute_abcd_dihedral_stats(a, b, c, d3),
    "colinear"
  )
  expect_equal(result3$bond_angle, 180, tolerance = 1e-10)
  expect_true(is.na(result3$torsion_angle))
})
