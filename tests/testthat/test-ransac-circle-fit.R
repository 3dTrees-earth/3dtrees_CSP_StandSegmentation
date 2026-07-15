test_that("RANSAC circle fitting avoids per-iteration output redirection", {
  body_text <- paste(deparse(body(ransac_circle_fit)), collapse = "\n")
  expect_false(grepl("suppress_cat", body_text, fixed = TRUE))
})

test_that("RANSAC circle fitting preserves partial-arc results without loop overhead", {
  theta <- seq(0, pi, length.out = 50L)
  points <- cbind(
    10 + 0.3 * cos(theta) + sin(seq_along(theta)) * 0.001,
    20 + 0.3 * sin(theta) + cos(seq_along(theta)) * 0.001
  )

  set.seed(42L)
  result <- ransac_circle_fit(
    points,
    n_iterations = 100L,
    distance_threshold = 0.01,
    min_inliers = 3L
  )

  expect_equal(as.numeric(result$circle), c(10.00015, 20.00363, 0.2983673), tolerance = 1e-5)
  expect_equal(result$inliers, 50L)
  expect_equal(result$angle_segs, 20L)
})
