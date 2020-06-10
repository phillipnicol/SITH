
test_that("C++ code rejects improper parameter values", {
  expect_error(simulateTumor(N = 0))
  expect_error(simulateTumor(s = -0.1))
  expect_error(simulateTumor(b = -0.1))
  expect_error(simulateTumor(d = -0.1))
  expect_error(simulateTumor(u = -0.1))
  expect_error(simulateTumor(du = -0.1))
  expect_error(simulateTumor(du = 1.1))
})

test_that("Simulation is reproducible", {
  set.seed(116776544)
  out <- simulateTumor(N = 200, verbose = FALSE, du = 0.5)
  
  expect_equal(nrow(out$cell_ids), 200)
  expect_equal(mean(out$cell_ids$nmuts), 1.675)
  expect_equal(length(out$drivers), 12)
  expect_equal(nrow(out$muts), 19)
  expect_equal(out$muts$MAF[5], 0)
  expect_equal(out$time, 43.879, tolerance = 0.1)
  
  out <- simulateTumor(N = 200, verbose = F, du = 1.0)
  expect_equal(nrow(out$muts) - 1, length(out$drivers))
})

test_that("The alleles data frame is properly constructed", {
  out <- simulateTumor(N = 50, verbose = F)
  
  expect_equal(50, sum(out$alleles$count))
})

test_that("The wild type allele is grey", {
  out <- simulateTumor(N = 50, verbose = F)
  
  expect_equal(out$color_scheme[1], "#808080FF")
})