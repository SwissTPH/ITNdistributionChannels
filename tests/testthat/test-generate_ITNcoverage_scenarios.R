test_that("calculate_cov returns correct structure and values", {
  # Test typical usage
  result <- calculate_cov(L = 2, kappa = 2)
  expect_type(result, "list")
  expect_named(result, c("ITNdecay", "average_cov"))
  expect_s3_class(result$ITNdecay, "data.frame")
  expect_equal(nrow(result$ITNdecay), 20 * 365)
  expect_equal(names(result$ITNdecay), c("time", "decay"))
  expect_type(result$average_cov, "double")
  expect_gt(result$average_cov, 0)
  expect_lt(result$average_cov, 1)

  # Test edge case: very short half-life
  result_short <- calculate_cov(L = 0.5, kappa = 2)
  expect_lt(result_short$average_cov, result$average_cov)

  # Test edge case: very long half-life
  result_long <- calculate_cov(L = 10, kappa = 2)
  expect_gt(result_long$average_cov, result$average_cov)

  # Test edge case: kappa = 1 (exponential decay)
  result_exp <- calculate_cov(L = 2, kappa = 1)
  expect_gt(result_exp$average_cov, 0)
})


test_that("generate_ITN_distribution works as expected", {
  # Test typical usage
  result <- generate_ITN_distribution(halflife_weibull = 2, shape_weibull = 2, coverage = 0.8, frequency = 3)
  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("time", "cov"))
  expect_equal(nrow(result), 20 * 365)

  # Test edge case: zero coverage
  result_zero <- generate_ITN_distribution(halflife_weibull = 2, shape_weibull = 2, coverage = 0, frequency = 3)
  expect_equal(sum(result_zero$cov), 0)

  # Test edge case: very frequent distribution
  result_freq <- generate_ITN_distribution(halflife_weibull = 2, shape_weibull = 2, coverage = 0.8, frequency = 1)
  expect_gt(mean(result_freq$cov), mean(result$cov))
})


test_that("generate_massCampaign_UniformAllAges returns correct output", {
  # Test typical usage
  result <- generate_massCampaign_UniformAllAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, frequency = 3)
  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("time", "cov","use_total", "reach", "max_usage", "halflife"))
  expect_equal(nrow(result), 20 * 365)
  expect_lt(max(result$use_total), 0.81)

  # Test edge case: ageGroupProp < 1
  result_prop <- generate_massCampaign_UniformAllAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, frequency = 3, ageGroupProp = 0.5)
  expect_lt(max(result_prop$use_total), max(result$use_total))

  # Test edge case: lower reach
  result_07 <- generate_massCampaign_UniformAllAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.7, max_usage = 0.8, frequency = 3, ageGroupProp = 0.5)
  expect_lt(max(result_07$use_total), max(result$use_total))
})


test_that("generate_massCampaign_SpecificAges works for target and other groups", {
  # Test typical usage
  result <- generate_massCampaign_SpecificAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, frequency = 3, ageGroup = "U5", otherageGroup = "others", ageGroupProp = 0.16)
  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("time", "use_U5", "use_others", "reach", "max_usage", "halflife"))
  expect_equal(nrow(result), 20 * 365)
  expect_gt(mean(result$use_U5), mean(result$use_others))

  # Test full_output = TRUE
  result_full <- generate_massCampaign_SpecificAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, frequency = 3, ageGroup = "U5", otherageGroup = "others", ageGroupProp = 0.16, full_output = TRUE)
  expect_equal(names(result_full), c("time", "totcov_targetgroup","use_targetpop","totcov_othergroup","spillover_cov_othergroup" ,"use_otherpop"))
})



test_that("generate_continuousDistr_uniformAllAges returns correct output", {
  # Test typical usage
  result <- generate_continuousDistr_uniformAllAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8)
  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("time","cov", "use_total", "reach", "max_usage", "halflife"))
  expect_equal(nrow(result), 20 * 365)
  expect_lt(max(result$use_total), 0.8)
})


test_that("generate_continuousDistr_SpecificAges works for target and other groups", {
  # Test typical usage
  result <- generate_continuousDistr_SpecificAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, ageGroup = "U1", otherageGroup = "others", ageGroupProp = 0.035)
  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("time", "use_U1", "use_others", "reach", "max_usage", "halflife"))
  expect_equal(nrow(result), 20 * 365)
  expect_gt(mean(result$use_U1), mean(result$use_others))

  # Test full_output = TRUE
  result_full <- generate_continuousDistr_SpecificAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, ageGroup = "U1", otherageGroup = "others", ageGroupProp = 0.035, full_output = TRUE)
  expect_equal(names(result_full), c("time", "totcov_targetgroup","use_targetpop","totcov_othergroup","spillover_cov_othergroup" ,"use_otherpop"))
})


test_that("get_average_coverage_continuous returns correct average", {
  scenario <- data.frame(
    time = 1:730,
    halflife = rep(2, 730),
    reach = rep(0.8, 730),
    max_usage = rep(0.8, 730),
    use_total = runif(730, 0, 0.8)
  )
  result <- get_average_coverage_continuous(scenario, year = 1)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(names(result), c("halflife", "reach", "max_usage", "use_total"))
})

test_that("get_start_coverage_mass filters correctly", {
  scenario <- data.frame(
    time = 1:730,
    use_total = runif(730, 0, 0.8)
  )
  result <- get_start_coverage_mass(scenario, year = 1)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$time, 366)
})

test_that("get_average_coverage_mass returns correct average", {
  scenario <- data.frame(
    time = 1:730,
    halflife = rep(2, 730),
    reach = rep(0.8, 730),
    max_usage = rep(0.8, 730),
    use_total = runif(730, 0, 0.8)
  )
  result <- get_average_coverage_mass(scenario, year_start = 0, year_end = 1)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(names(result), c("halflife", "reach", "max_usage", "use_total"))
})



test_that("generate_continuousDistr_uniformAllAges and generate_continuousDistr_SpecificAges align for newborns", {
  age_prop=0.035
  # Test typical usage
  result_unif <- generate_continuousDistr_uniformAllAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8,ageGroupProp=age_prop)
  result_targ <- generate_continuousDistr_SpecificAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, ageGroup = "U1", otherageGroup = "others", ageGroupProp = age_prop)

  result_targ$use_total=result_targ$use_U1*age_prop + result_targ$use_others*(1-age_prop)
  expect_equal(result_unif$use_total, result_targ$use_total)
})


test_that("generate_continuousDistr_uniformAllAges and generate_continuousDistr_SpecificAges align for school", {
  age_prop=0.27*0.5
  # Test typical usage
  result_unif <- generate_continuousDistr_uniformAllAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8,ageGroupProp=age_prop)
  result_targ <- generate_continuousDistr_SpecificAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, ageGroup = "U1", otherageGroup = "others", ageGroupProp = age_prop)

  result_targ$use_total=result_targ$use_U1*age_prop + result_targ$use_others*(1-age_prop)
  expect_equal(result_unif$use_total, result_targ$use_total)
})




test_that("generate_continuousDistr_uniformAllAges and generate_continuousDistr_SpecificAges align for newborns", {
  age_prop=0.035
  # Test typical usage
  result_unif <- generate_continuousDistr_uniformAllAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8,ageGroupProp=age_prop)
  result_targ <- generate_continuousDistr_SpecificAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, ageGroup = "U1", otherageGroup = "others", ageGroupProp = age_prop)

  result_targ$use_total=result_targ$use_U1*age_prop + result_targ$use_others*(1-age_prop)
  expect_equal(result_unif$use_total, result_targ$use_total)
})


test_that("generate_massCampaign_UniformAllAges and generate_massCampaign_SpecificAges align for u5", {
  age_prop=0.16
  # Test typical usage
  result_unif <- generate_massCampaign_UniformAllAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, frequency = 3, ageGroupProp = age_prop)
  result_targ <- generate_massCampaign_SpecificAges(halflife_weibull = 2, shape_weibull = 2, reach = 0.8, max_usage = 0.8, frequency = 3, ageGroup = "U5", otherageGroup = "others", ageGroupProp = age_prop)

  result_targ$use_total=result_targ$use_U5*age_prop + result_targ$use_others*(1-age_prop)
  expect_equal(result_unif$use_total, result_targ$use_total)
})
