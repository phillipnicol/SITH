
test_that("C++ code rejects improper parameter values", {
  expect_error(simulateTumor(max_pop = 0))
  expect_error(simulateTumor(selective_adv = -0.1))
  expect_error(simulateTumor(div_rate = -0.1))
  expect_error(simulateTumor(death_rate = -0.1))
  expect_error(simulateTumor(mut_rate = -0.1))
  expect_error(simulateTumor(driver_prob = -0.1))
  expect_error(simulateTumor(driver_prob = 1.1))
})

test_that("Simulation is reproducible", {
  set.seed(116776544)
  out <- simulateTumor(max_pop = 200, verbose = FALSE, driver_prob = 0.5)
  
  expect_equal(nrow(out$cell_ids), 200)
  expect_equal(mean(out$cell_ids$nmuts), 1.11)
  expect_equal(length(out$drivers), 10)
  expect_equal(nrow(out$muts), 18)
  expect_equal(out$muts$MAF[5], 0)
  expect_equal(out$time, 105.7691, tolerance = 0.1)
  
  out <- simulateTumor(max_pop = 200, verbose = F, driver_prob = 1.0)
  expect_equal(nrow(out$muts) - 1, length(out$drivers))
})

test_that("The genotypes data frame is properly constructed", {
  out <- simulateTumor(max_pop = 50, verbose = F)
  
  expect_equal(50, sum(out$genotypes$count))
})

test_that("The wild type allele is grey", {
  out <- simulateTumor(max_pop = 50, verbose = F)
  
  expect_equal(out$color_scheme[1], "#808080FF")
})