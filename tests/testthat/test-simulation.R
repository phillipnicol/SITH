
test_that("C++ code rejects improper parameter values", {
  expect_error(simulateTumor(N = 0))
  expect_error(simulateTumor(s = -0.1))
  expect_error(simulateTumor(b = -0.1))
  expect_error(simulateTumor(d = -0.1))
  expect_error(simulateTumor(u = -0.1))
  expect_error(simulateTumor(du = -0.1))
  expect_error(simulateTumor(du = 1.1))
})

test_that("Simulation behaves as expected", {
  set.seed(116776544, kind = "Mersenne-Twister", normal.kind = "Inversion")
  out <- simulateTumor(N = 200, verbose = FALSE, du = 0.5)
  
  expect_equal(nrow(out$cell_ids), 200)
  expect_equal(mean(out$cell_ids$nmuts), 1.1)
  expect_equal(length(out$drivers), 6)
  expect_equal(nrow(out$muts), 11)
  expect_equal(out$muts$MAF[5], 0.045)
  
  out <- simulateTumor(N = 200, verbose = F, du = 1.0)
  expect_equal(nrow(out$muts) - 1, length(out$drivers))
})

test_that("The alleles data frame is properly constructed", {
  out <- simulateTumor(N = 50, verbose = F)
  
  expect_equal(50, sum(out$alleles$count))
})