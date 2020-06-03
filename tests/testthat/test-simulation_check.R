
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
  set.seed(116776544)
  out <- simulateTumor(N = 200, verbose = F, du = 0.1)
  
  expect_equal(nrow(out$cell_ids), 200)
  expect_equal(mean(out$cell_ids$nmuts), 1.2)
  expect_equal(length(out$drivers), 1)
  expect_equal(nrow(out$muts), 11)
  expect_equal(out$muts$MAF[4], 0.12)
})

test_that("The alleles data frame is properly constructed", {
  out <- simulateTumor(N = 50, verbose = F)
  
  expect_equal(50, sum(out$alleles$count))
})