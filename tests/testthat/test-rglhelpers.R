# tests/testthat/test-convert_tidy_segments_to_interleaved.R

test_that("single 3D segment with metadata is interleaved correctly", {
  df <- data.frame(
    x = 0, y = 0, z = 0,
    xend = 1, yend = 1, zend = 1,
    id = 42,
    stringsAsFactors = FALSE
  )

  out <- convert_tidy_segments_to_interleaved(df)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_equal(names(out), c("point", "x", "y", "z", "id"))

  expect_equal(out$point, c("start", "end"))
  expect_equal(out$x, c(0, 1))
  expect_equal(out$y, c(0, 1))
  expect_equal(out$z, c(0, 1))
  expect_equal(out$id, c(42, 42))
})

test_that("multiple segments are interleaved in start/end order with metadata replicated", {
  df <- data.frame(
    x = c(0, 10), y = c(0, 10), z = c(0, 10),
    xend = c(1, 11), yend = c(2, 12), zend = c(3, 13),
    grp = factor(c("a", "b")),
    label = c("first", "second"),
    stringsAsFactors = FALSE
  )

  out <- convert_tidy_segments_to_interleaved(df)

  expect_equal(nrow(out), 4)
  expect_equal(out$point, c("start", "end", "start", "end"))

  # metadata replicated in pairs
  expect_identical(as.character(out$grp), c("a", "a", "b", "b"))
  expect_equal(out$label, c("first", "first", "second", "second"))

  # coordinates interleave per row
  expect_equal(out$x, c(0, 1, 10, 11))
  expect_equal(out$y, c(0, 2, 10, 12))
  expect_equal(out$z, c(0, 3, 10, 13))
})

test_that("custom 2D coordinates and custom end suffix work", {
  df <- data.frame(
    x = c(0, 1),
    y = c(0, 1),
    x_end = c(2, 3),
    y_end = c(2, 3),
    id = c("s1", "s2"),
    stringsAsFactors = FALSE
  )

  out <- convert_tidy_segments_to_interleaved(
    df,
    coord = c("x", "y"),
    end_suffix = "_end"
  )

  expect_equal(names(out), c("point", "x", "y", "id"))
  expect_equal(nrow(out), 4)
  expect_equal(out$point, c("start", "end", "start", "end"))
  expect_equal(out$x, c(0, 2, 1, 3))
  expect_equal(out$y, c(0, 2, 1, 3))
  expect_equal(out$id, c("s1", "s1", "s2", "s2"))
})

test_that("matrix input with proper colnames is supported", {
  m <- cbind(
    x = c(0, 5),
    y = c(0, 6),
    z = c(0, 7),
    xend = c(1, 8),
    yend = c(2, 9),
    zend = c(3, 10)
  )
  # matrix -> data.frame inside function via subsetting; here we pass matrix directly
  out <- convert_tidy_segments_to_interleaved(m, coord = c("x", "y", "z"))

  expect_equal(nrow(out), 4)
  expect_equal(names(out), c("point", "x", "y", "z"))
  expect_equal(out$point, c("start", "end", "start", "end"))
  expect_equal(out$x, c(0, 1, 5, 8))
  expect_equal(out$y, c(0, 2, 6, 9))
  expect_equal(out$z, c(0, 3, 7, 10))
})

test_that("error is thrown when required columns are missing", {
  df_missing <- data.frame(
    x = 0, y = 0, z = 0,
    xend = 1, zend = 1,  # yend missing
    stringsAsFactors = FALSE
  )
  expect_error(
    convert_tidy_segments_to_interleaved(df_missing),
    "Missing required columns"
  )
})

test_that("zero-row input returns zero rows with correct columns", {
  df0 <- data.frame(
    x = numeric(0), y = numeric(0), z = numeric(0),
    xend = numeric(0), yend = numeric(0), zend = numeric(0),
    stringsAsFactors = FALSE
  )

  out <- convert_tidy_segments_to_interleaved(df0)

  expect_equal(nrow(out), 0)
  expect_equal(names(out), c("point", "x", "y", "z"))
})

test_that("output column order is point, coords, then meta", {
  df <- data.frame(
    x = 0, y = 0, z = 0,
    xend = 1, yend = 1, zend = 1,
    meta1 = "A",
    meta2 = 99L,
    stringsAsFactors = FALSE
  )

  out <- convert_tidy_segments_to_interleaved(df)

  expect_equal(names(out), c("point", "x", "y", "z", "meta1", "meta2"))
})

test_that("metadata types are preserved", {
  df <- data.frame(
    x = 0, y = 0, z = 0,
    xend = 1, yend = 1, zend = 1,
    f = factor("lvl1"),
    i = 7L,
    d = as.Date("2020-01-01"),
    stringsAsFactors = FALSE
  )

  out <- convert_tidy_segments_to_interleaved(df)

  expect_true(is.factor(out$f))
  expect_true(is.integer(out$i))
  expect_s3_class(out$d, "Date")
  # replicated values
  expect_equal(out$f, factor(c("lvl1", "lvl1")))
  expect_equal(out$i, c(7L, 7L))
  expect_equal(out$d, as.Date(c("2020-01-01", "2020-01-01")))
})
