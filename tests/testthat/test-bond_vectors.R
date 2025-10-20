test_that("Bond angle calculation works", {
  positions <- list(
    a = c(x = 32.9024, y = 20.7048, z = 31.1175),
    b = c(x = 33.7508, y = 20.7697, z = 32.2441),
    c = c(x = 33.6998, y = 19.7405, z = 33.0971)
  )

  bond_vector <- bond_vector_from_internal(
    positions$a,
    positions$b,
    positions$c,
    bond_length = 1.987,
    bond_angle = 121.05,
    torsion_angle = 176.88,
    degrees = TRUE
  )

  expect_equal(bond_vector, c(x = 1.07089999999999, y = 0.00189999999999912, z = 1.6736))
})
