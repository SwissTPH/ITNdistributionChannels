#####################################
# check functions individually
#####################################
test_that("generate_MassAndContinuousDistr_uniformAllAges: typical usage", {
  result <- generate_MassAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.035
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("time", "use_total", "reach_mass", "reach_cont", "max_usage", "halflife") %in% colnames(result)))
  expect_true(all(result$use_total <= 0.8)) # Should not exceed max_usage
})

test_that("generate_MassAndContinuousDistr_uniformAllAges: zero reach", {
  result <- generate_MassAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_mass = 0,
    reach_cont = 0,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.035
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(result$use_total == 0)) # No reach, no usage
})

test_that("generate_MassAndContinuousDistr_SpecificAges: typical usage", {
  result <- generate_MassAndContinuousDistr_SpecificAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.035,
    ageGroup_cont = "U1",
    otherageGroup_cont = "others"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("time", "use_U1", "use_others", "reach_mass", "reach_cont", "max_usage", "halflife") %in% colnames(result)))
  expect_true(all(result$use_U1 <= 0.8))
  expect_true(all(result$use_others <= 0.8))
})

test_that("generate_ContAndContinuousDistr_uniformAllAges: typical usage", {
  result <- generate_ContAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_cont1 = 0.8,
    reach_cont2 = 0.5,
    max_usage = 0.8,
    ageGroupProp_cont1 = 0.035,
    ageGroupProp_cont2 = 0.27
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(result$use_total <= 0.8))
})

test_that("generate_ContAndContinuousDistr_SpecificAges: typical usage", {
  result <- generate_ContAndContinuousDistr_SpecificAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_cont1 = 0.8,
    reach_cont2 = 0.5,
    max_usage = 0.8,
    ageGroupProp_cont1 = 0.035,
    ageGroupProp_cont2 = 0.27,
    ageGroupCont1 = "U1",
    otherageGroupCont1 = "others",
    ageGroupCont2 = "school",
    otherageGroupCont2 = "others"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("time", "use_U1", "use_school", "use_others") %in% colnames(result)))
  expect_true(all(result$use_U1 <= 0.8))
  expect_true(all(result$use_school <= 0.8))
  expect_true(all(result$use_others <= 0.8))
})

test_that("generate_AlternateMassCont_uniformAllAges: typical usage", {
  result <- generate_AlternateMassCont_uniformAllAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(result$use_total <= 0.8))
})

test_that("generate_AlternateMassCont_SpecificAges: typical usage", {
  result <- generate_AlternateMassCont_SpecificAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.035,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("time", "use_school", "use_others", "reach_mass", "reach_cont", "max_usage", "halflife") %in% colnames(result)))
  expect_true(all(result$use_school <= 0.8))
  expect_true(all(result$use_others <= 0.8))
})

test_that("generate_AlternateMassContRoutine_uniformAllAges: typical usage", {
  result <- generate_AlternateMassContRoutine_uniformAllAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(result$use_total <= 0.8))
})

test_that("generate_AlternateMassContRoutine_SpecificAges: typical usage", {
  result <- generate_AlternateMassContRoutine_SpecificAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("time", "use_school", "use_U1", "use_others") %in% colnames(result)))
  expect_true(all(result$use_school <= 0.8))
  expect_true(all(result$use_U1 <= 0.8))
  expect_true(all(result$use_others <= 0.8))
})

test_that("generate_AlternateContRoutine_uniformAllAges: typical usage", {
  result <- generate_AlternateContRoutine_uniformAllAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(result$use_total <= 0.8))
})

test_that("generate_AlternateContRoutine_SpecificAges: typical usage", {
  result <- generate_AlternateContRoutine_SpecificAges(
    halflife_weibull = 2.1,
    shape_weibull = 2,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("time", "use_school", "use_U1", "use_others") %in% colnames(result)))
  expect_true(all(result$use_school <= 0.8))
  expect_true(all(result$use_U1 <= 0.8))
  expect_true(all(result$use_others <= 0.8))
})


#####################################################
# Consistency check: uniform vs specific age groups
#####################################################

# Consistency check: Mass + Continuous
test_that("generate_MassAndContinuousDistr_uniformAllAges and generate_MassAndContinuousDistr_SpecificAges align for school", {
  age_prop = 0.27
  result_unif <- generate_MassAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop
  )%>% dplyr::arrange(time)

  result_targ <- generate_MassAndContinuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  )
  result_targ$use_total = result_targ$use_school * age_prop + result_targ$use_others * (1 - age_prop)
  expect_equal(result_unif$use_total, result_targ$use_total, tolerance = 1e-6)
})


# Consistency check: Continuous + Continuous
test_that("generate_ContAndContinuousDistr_uniformAllAges and generate_ContAndContinuousDistr_SpecificAges align for school and U1", {
  age_prop1 = 0.035
  age_prop2 = 0.27
  result_unif <- generate_ContAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont1 = 0.8,
    reach_cont2 = 0.5,
    max_usage = 0.8,
    ageGroupProp_cont1 = age_prop1,
    ageGroupProp_cont2 = age_prop2
  ) %>% dplyr::arrange(time)

  result_targ <- generate_ContAndContinuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont1 = 0.8,
    reach_cont2 = 0.5,
    max_usage = 0.8,
    ageGroupProp_cont1 = age_prop1,
    ageGroupProp_cont2 = age_prop2,
    ageGroupCont1 = "U1",
    otherageGroupCont1 = "others",
    ageGroupCont2 = "school",
    otherageGroupCont2 = "others"
  )
  result_targ$use_total = result_targ$use_U1 * age_prop1 + result_targ$use_school * age_prop2 + result_targ$use_others * (1 - age_prop1 - age_prop2)
  expect_equal(result_unif$use_total, result_targ$use_total, tolerance = 1e-3)
})

# Consistency check: Alternate Mass + Continuous
test_that("generate_AlternateMassCont_uniformAllAges and generate_AlternateMassCont_SpecificAges align for school", {
  age_prop = 0.27
  result_unif <- generate_AlternateMassCont_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop
  ) %>% dplyr::arrange(time)
  result_targ <- generate_AlternateMassCont_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  )
  result_targ$use_total = result_targ$use_school * age_prop + result_targ$use_others * (1 - age_prop)
  expect_equal(result_unif$use_total, result_targ$use_total, tolerance = 1e-6)
})

# Consistency check: Alternate Mass + Continuous + Routine
test_that("generate_AlternateMassContRoutine_uniformAllAges and generate_AlternateMassContRoutine_SpecificAges align for school and U1", {
  age_prop_cont = 0.27
  age_prop_routine = 0.035
  result_unif <- generate_AlternateMassContRoutine_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop_cont,
    ageGroupProp_routine = age_prop_routine
  ) %>% dplyr::arrange(time)
  result_targ <- generate_AlternateMassContRoutine_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop_cont,
    ageGroupProp_routine = age_prop_routine,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  )
  result_targ$use_total = result_targ$use_school * age_prop_cont + result_targ$use_U1 * age_prop_routine + result_targ$use_others * (1 - age_prop_cont - age_prop_routine)
  expect_equal(result_unif$use_total, result_targ$use_total, tolerance = 1e-4)
})

# Consistency check: Alternate Continuous + Routine
test_that("generate_AlternateContRoutine_uniformAllAges and generate_AlternateContRoutine_SpecificAges align for school and U1", {
  age_prop_cont = 0.27
  age_prop_routine = 0.035
  result_unif <- generate_AlternateContRoutine_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop_cont,
    ageGroupProp_routine = age_prop_routine
  ) %>% dplyr::arrange(time)
  result_targ <- generate_AlternateContRoutine_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop_cont,
    ageGroupProp_routine = age_prop_routine,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  )
  result_targ$use_total = result_targ$use_school * age_prop_cont + result_targ$use_U1 * age_prop_routine + result_targ$use_others * (1 - age_prop_cont - age_prop_routine)
  expect_equal(result_unif$use_total, result_targ$use_total, tolerance = 1e-3)
})



# Consistency check: Alternate Continuous
test_that("generate_AlternateCont_uniformAllAges and generate_AlternateCont_SpecificAges align for school and U1", {
  age_prop_cont = 0.27
  result_unif <- generate_AlternateCont_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop_cont
  ) %>% dplyr::arrange(time)
  result_targ <- generate_AlternateCont_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = age_prop_cont,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  )
  result_targ$use_total = result_targ$use_school * age_prop_cont  + result_targ$use_others * (1 - age_prop_cont)
  expect_equal(result_unif$use_total, result_targ$use_total, tolerance = 1e-6)
})
##############################


# Consistency check: Mass + Continuous (uniform) with reach_mass=0 should equal Continuous only
test_that("generate_MassAndContinuousDistr_uniformAllAges with reach_mass=0 equals generate_continuousDistr_uniformAllAges", {
  result_mass_cont <- generate_MassAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27
  ) %>% dplyr::arrange(time)
  result_cont <- generate_continuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.5,
    max_usage = 0.8,
    ageGroupProp = 0.27
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont$use_total, result_cont$use_total, tolerance = 1e-6)
})

# Consistency check: Mass + Continuous (specific) with reach_mass=0 should equal Continuous only
test_that("generate_MassAndContinuousDistr_SpecificAges with reach_mass=0 equals generate_continuousDistr_SpecificAges", {
  result_mass_cont <- generate_MassAndContinuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  ) %>% dplyr::arrange(time)
  result_cont <- generate_continuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.5,
    max_usage = 0.8,
    ageGroup = "school",
    otherageGroup = "others",
    ageGroupProp = 0.27
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont$use_school, result_cont$use_school, tolerance = 1e-6)
  expect_equal(result_mass_cont$use_others, result_cont$use_others, tolerance = 1e-6)
})

# Consistency check: Alternate Mass + Continuous (uniform) with reach_mass=0 should equal Alternate Continuous only
test_that("generate_AlternateMassCont_uniformAllAges with reach_mass=0 equals generate_AlternateCont_uniformAllAges", {
  result_mass_cont <- generate_AlternateMassCont_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27
  ) %>% dplyr::arrange(time)
  result_cont <- generate_AlternateCont_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    rename = TRUE
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont$use_total, result_cont$use_total, tolerance = 1e-6)
})

# Consistency check: Alternate Mass + Continuous (specific) with reach_mass=0 should equal Alternate Continuous only
test_that("generate_AlternateMassCont_SpecificAges with reach_mass=0 equals generate_AlternateCont_SpecificAges", {
  result_mass_cont <- generate_AlternateMassCont_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  ) %>% dplyr::arrange(time)

  result_cont <- generate_AlternateCont_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    rename = TRUE
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont$use_school, result_cont$use_school, tolerance = 1e-6)
  expect_equal(result_mass_cont$use_others, result_cont$use_others, tolerance = 1e-6)
})

# Consistency check: Alternate Mass + Continuous + Routine (uniform) with reach_mass=0 should equal Alternate Continuous + Routine only
test_that("generate_AlternateMassContRoutine_uniformAllAges with reach_mass=0 equals generate_AlternateContRoutine_uniformAllAges", {
  result_mass_cont_routine <- generate_AlternateMassContRoutine_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035
  ) %>% dplyr::arrange(time)
  result_cont_routine <- generate_AlternateContRoutine_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont_routine$use_total, result_cont_routine$use_total, tolerance = 1e-6)
})

# Consistency check: Alternate Mass + Continuous + Routine (specific) with reach_mass=0 should equal Alternate Continuous + Routine only
test_that("generate_AlternateMassContRoutine_SpecificAges with reach_mass=0 equals generate_AlternateContRoutine_SpecificAges", {
  result_mass_cont_routine <- generate_AlternateMassContRoutine_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  ) %>% dplyr::arrange(time)
  result_cont_routine <- generate_AlternateContRoutine_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0.5,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont_routine$use_school, result_cont_routine$use_school, tolerance = 1e-6)
  expect_equal(result_mass_cont_routine$use_U1, result_cont_routine$use_U1, tolerance = 1e-6)
  expect_equal(result_mass_cont_routine$use_others, result_cont_routine$use_others, tolerance = 1e-6)
})


# Consistency check: Mass + Continuous (uniform) with reach_cont=0 should equal Mass only
test_that("generate_MassAndContinuousDistr_uniformAllAges with reach_cont=0 equals generate_massCampaign_UniformAllAges", {
  result_mass_cont <- generate_MassAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27
  ) %>% dplyr::arrange(time)
  result_mass <- generate_massCampaign_UniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.8,
    max_usage = 0.8,
    frequency = 3,
    ageGroupProp = 1
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont$use_total, result_mass$use_total, tolerance = 1e-6)
})

# Consistency check: Mass + Continuous (specific) with reach_cont=0 should equal Mass only
test_that("generate_MassAndContinuousDistr_SpecificAges with reach_cont=0 equals generate_massCampaign_UniformAllAges", {
  result_mass_cont <- generate_MassAndContinuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  ) %>% dplyr::arrange(time)
  result_mass <- generate_massCampaign_UniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.8,
    max_usage = 0.8,
    frequency = 3,
    ageGroupProp = 1
  ) %>% dplyr::arrange(time)
  # For specific ages, use_total is not returned, so we calculate it from the mass campaign
  result_mass_cont$use_total = result_mass_cont$use_school * 0.27 + result_mass_cont$use_others * 0.73
  expect_equal(result_mass_cont$use_total, result_mass$use_total, tolerance = 1e-6)
})

# Consistency check: Continuous + Continuous (uniform) with reach_cont2=0 should equal Continuous only (cont1)
test_that("generate_ContAndContinuousDistr_uniformAllAges with reach_cont2=0 equals generate_continuousDistr_uniformAllAges (cont1)", {
  result_cont_cont <- generate_ContAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont1 = 0.8,
    reach_cont2 = 0,
    max_usage = 0.8,
    ageGroupProp_cont1 = 0.035,
    ageGroupProp_cont2 = 0.27
  ) %>% dplyr::arrange(time)
  result_cont <- generate_continuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.8,
    max_usage = 0.8,
    ageGroupProp = 0.035
  ) %>% dplyr::arrange(time)
  expect_equal(result_cont_cont$use_total, result_cont$use_total, tolerance = 1e-6)
})

# Consistency check: Continuous + Continuous (specific) with reach_cont2=0 should equal Continuous only (cont1)
test_that("generate_ContAndContinuousDistr_SpecificAges with reach_cont2=0 equals generate_continuousDistr_SpecificAges (cont1)", {
  result_cont_cont <- generate_ContAndContinuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont1 = 0.8,
    reach_cont2 = 0,
    max_usage = 0.8,
    ageGroupProp_cont1 = 0.035,
    ageGroupProp_cont2 = 0.27,
    ageGroupCont1 = "U1",
    otherageGroupCont1 = "others",
    ageGroupCont2 = "school",
    otherageGroupCont2 = "others"
  ) %>% dplyr::arrange(time)
  result_cont <- generate_continuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.8,
    max_usage = 0.8,
    ageGroup = "U1",
    otherageGroup = "others",
    ageGroupProp = 0.035
  ) %>% dplyr::arrange(time)
  expect_equal(result_cont_cont$use_U1, result_cont$use_U1, tolerance = 1e-6)
  expect_equal(result_cont_cont$use_others, result_cont$use_others, tolerance = 1e-6)
})




# Consistency check: Continuous + Continuous (uniform) with reach_cont1=0 should equal Continuous only (cont1)
test_that("generate_ContAndContinuousDistr_uniformAllAges with reach_cont1=0 equals generate_continuousDistr_uniformAllAges (cont1)", {
  result_cont_cont <- generate_ContAndContinuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont1 = 0,
    reach_cont2 = 0.8,
    max_usage = 0.8,
    ageGroupProp_cont1 = 0.035,
    ageGroupProp_cont2 = 0.27
  ) %>% dplyr::arrange(time)
  result_cont <- generate_continuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.8,
    max_usage = 0.8,
    ageGroupProp = 0.27
  ) %>% dplyr::arrange(time)
  expect_equal(result_cont_cont$use_total, result_cont$use_total, tolerance = 1e-6)
})

# Consistency check: Continuous + Continuous (specific) with reach_cont1=0 should equal Continuous only (cont1)
test_that("generate_ContAndContinuousDistr_SpecificAges with reach_cont1=0 equals generate_continuousDistr_SpecificAges (cont1)", {
  result_cont_cont <- generate_ContAndContinuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont1 = 0,
    reach_cont2 = 0.8,
    max_usage = 0.8,
    ageGroupProp_cont1 = 0.035,
    ageGroupProp_cont2 = 0.27,
    ageGroupCont1 = "U1",
    otherageGroupCont1 = "others",
    ageGroupCont2 = "school",
    otherageGroupCont2 = "others"
  ) %>% dplyr::arrange(time)
  result_cont <- generate_continuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.8,
    max_usage = 0.8,
    ageGroup = "school",
    otherageGroup = "others",
    ageGroupProp = 0.27
  ) %>% dplyr::arrange(time)
  expect_equal(result_cont_cont$use_school, result_cont$use_school, tolerance = 1e-6)
  expect_equal(result_cont_cont$use_others, result_cont$use_others, tolerance = 1e-6)
})


# Consistency check: Alternate Mass + Continuous (uniform) with reach_cont=0 should equal Mass only
test_that("generate_AlternateMassCont_uniformAllAges with reach_cont=0 equals generate_massCampaign_UniformAllAges", {
  result_mass_cont <- generate_AlternateMassCont_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27
  ) %>% dplyr::arrange(time)
  result_mass <- generate_massCampaign_UniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.8,
    max_usage = 0.8,
    frequency = 3,
    ageGroupProp = 1
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont$use_total, result_mass$use_total, tolerance = 1e-6)
})

# Consistency check: Alternate Mass + Continuous (specific) with reach_cont=0 should equal Mass only
test_that("generate_AlternateMassCont_SpecificAges with reach_cont=0 equals generate_massCampaign_UniformAllAges", {
  result_mass_cont <- generate_AlternateMassCont_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  ) %>% dplyr::arrange(time)
  result_mass <- generate_massCampaign_UniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.8,
    max_usage = 0.8,
    frequency = 3,
    ageGroupProp = 1
  ) %>% dplyr::arrange(time)
  result_mass_cont$use_total = result_mass_cont$use_school * 0.27 + result_mass_cont$use_others * 0.73
  expect_equal(result_mass_cont$use_total, result_mass$use_total, tolerance = 1e-6)
})



# Consistency check: Alternate Continuous + Routine (uniform) with reach_cont=0 should equal Routine only
test_that("generate_AlternateContRoutine_uniformAllAges with reach_cont=0 equals generate_continuousDistr_uniformAllAges (routine)", {
  result_cont_routine <- generate_AlternateContRoutine_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035
  ) %>% dplyr::arrange(time)
  result_routine <- generate_continuousDistr_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.6,
    max_usage = 0.8,
    ageGroupProp = 0.035
  ) %>% dplyr::arrange(time)
  expect_equal(result_cont_routine$use_total, result_routine$use_total, tolerance = 1e-6)
})

# Consistency check: Alternate Continuous + Routine (specific) with reach_cont=0 should equal Routine only
test_that("generate_AlternateContRoutine_SpecificAges with reach_cont=0 equals generate_continuousDistr_SpecificAges (routine)", {
  result_cont_routine <- generate_AlternateContRoutine_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_cont = 0,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  ) %>% dplyr::arrange(time)
  result_routine <- generate_continuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach = 0.6,
    max_usage = 0.8,
    ageGroup = "U1",
    otherageGroup = "others",
    ageGroupProp = 0.035
  ) %>% dplyr::arrange(time)
  expect_equal(result_cont_routine$use_U1, result_routine$use_U1, tolerance = 1e-6)
  expect_equal(result_cont_routine$use_others, result_routine$use_others, tolerance = 1e-6)
})


# Consistency check: Mass + Cont + Routine (uniform) with reach_routine=0 should equal Mass + Cont only
test_that("generate_AlternateMassContRoutine_uniformAllAges with reach_routine=0 equals generate_AlternateMassCont_uniformAllAges", {
  result_mass_cont_routine <- generate_AlternateMassContRoutine_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    reach_routine = 0,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035
  ) %>% dplyr::arrange(time)
  result_mass_cont <- generate_AlternateMassCont_uniformAllAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont_routine$use_total, result_mass_cont$use_total, tolerance = 1e-6)
})

# Consistency check: Mass + Cont + Routine (specific) with reach_routine=0 should equal Mass + Cont only
test_that("generate_AlternateMassContRoutine_SpecificAges with reach_routine=0 equals generate_AlternateMassCont_SpecificAges", {
  result_mass_cont_routine <- generate_AlternateMassContRoutine_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    reach_routine = 0,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  ) %>% dplyr::arrange(time)
  result_mass_cont <- generate_AlternateMassCont_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.5,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroup_cont = "school",
    otherageGroup_cont = "others"
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont_routine$use_school, result_mass_cont$use_school, tolerance = 1e-6)
  expect_equal(result_mass_cont_routine$use_others, result_mass_cont$use_others, tolerance = 1e-6)
})

# Consistency check: Mass + Cont + Routine (specific) with reach_cont=0 should equal Mass + Routine only
test_that("generate_AlternateMassContRoutine_SpecificAges with reach_cont=0 equals Mass + Routine only", {
  result_mass_cont_routine <- generate_AlternateMassContRoutine_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0,
    reach_routine = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.27,
    ageGroupProp_routine = 0.035,
    ageGroup_cont = "school",
    otherageGroup_cont = "others",
    ageGroup_routine = "U1",
    otherageGroup_routine = "others"
  ) %>% dplyr::arrange(time)
  # Use the mass + cont function with cont=0 and routine proportion as cont proportion
  result_mass_routine <- generate_MassAndContinuousDistr_SpecificAges(
    halflife_weibull = 2,
    shape_weibull = 2,
    reach_mass = 0.8,
    reach_cont = 0.6,
    max_usage = 0.8,
    frequency_mass = 3,
    ageGroupProp_cont = 0.035,
    ageGroup_cont = "U1",
    otherageGroup_cont = "others"
  ) %>% dplyr::arrange(time)
  expect_equal(result_mass_cont_routine$use_U1, result_mass_routine$use_U1, tolerance = 1e-6)
  expect_equal(result_mass_cont_routine$use_others, result_mass_routine$use_others, tolerance = 1e-6)
})

