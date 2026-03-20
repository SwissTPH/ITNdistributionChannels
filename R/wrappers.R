#' Coverage from Mass Campaign (All Ages, Decay Fit)
#'
#' Simulates a mass campaign with uniform coverage across all ages, fits a decay
#' model to the usage over years 12–15, and returns updated durability and coverage parameters.
#'
#' @param halflife Numeric. Initial half-life of the intervention.
#' @param reach Numeric. Proportion of population reached.
#' @param max_usage Numeric. Maximum achievable usage.
#' @param shape_weibull Numeric. Shape parameter of Weibull decay (default = 2).
#'
#' @return Named numeric vector with updated half-life, kappa, and usage.
#' @export
#'
getCov_MassCampaign=function(halflife, reach, max_usage, shape_weibull=2){
  scenario=generate_massCampaign_UniformAllAges(halflife_weibull = halflife,
                                                shape_weibull=shape_weibull,
                                                reach=reach, max_usage = max_usage,
                                                frequency=3)
  scenario2=scenario |>
    dplyr::filter(time>=12*365+1, time<15*365+1)

  df_for_decay=data.frame(value= scenario2$use_total, time=(1:length(scenario2$use_total))/365)
  my_updated_param=calculate_durability_param(df_for_decay)[1:3]
  plot_decay(as.numeric(unlist(my_updated_param)), df_for_decay)

  names(my_updated_param)=c("updated_halflife", "updated_kappa", "use_total")

  return(my_updated_param)
}

#' Coverage from Mass Campaign (Uniform, Under-5 Proportion)
#'
#' Same as \code{getCov_MassCampaign} but allows specifying the proportion of
#' under-5 population in a uniform-age simulation.
#'
#' @inheritParams getCov_MassCampaign
#' @param popU5 Numeric. Proportion of population under 5 (default = 0.16).
#'
#' @return Named numeric vector with updated half-life, kappa, and usage.
#' @export
getCov_MassCampaignU5unif=function(halflife, reach, max_usage,popU5=0.16, shape_weibull=2){

  scenario=generate_massCampaign_UniformAllAges(halflife_weibull = halflife,
                                                shape_weibull=shape_weibull,
                                                reach=reach, max_usage = max_usage,
                                                ageGroupProp=popU5,
                                                frequency=3)

  scenario2=scenario |>
    dplyr::filter(time>=12*365+1, time<15*365+1)

  df_for_decay=data.frame(value= scenario2$use_total, time=(1:length(scenario2$use_total))/365)
  my_updated_param=calculate_durability_param(df_for_decay)[1:3]
  plot_decay(as.numeric(unlist(my_updated_param)), df_for_decay)

  names(my_updated_param)=c("updated_halflife", "updated_kappa", "use_total")


  return(my_updated_param)
}


#' Coverage from Mass Campaign (Age-Specific: Under-5 vs Others)
#'
#' Simulates a mass campaign with separate coverage for under-5 and other age groups,
#' fits decay models, and returns updated parameters for both groups.
#'
#' @inheritParams getCov_MassCampaign
#' @param popU5 Numeric. Proportion of population under 5 (default = 0.16).
#'
#' @return Named numeric vector containing parameters for both groups.
#' @export
getCov_MassCampaignU5=function(halflife, reach, max_usage, popU5=0.16, shape_weibull=2){

  scenario=generate_massCampaign_SpecificAges(halflife_weibull = halflife,
                                              shape_weibull=shape_weibull,
                                              ageGroup="U5",otherageGroup="others",
                                              ageGroupProp=popU5, max_usage = max_usage,
                                              reach = reach, frequency=3)

  scenario2=scenario |>
    dplyr::filter(time>=12*365+1, time<15*365+1)

  df_for_decay_U5=data.frame(value= scenario2$use_U5, time=(1:length(scenario2$use_U5))/365)
  df_for_decay_others=data.frame(value= scenario2$use_others, time=(1:length(scenario2$use_others))/365)

  my_updated_param_U5=calculate_durability_param(df_for_decay_U5)[1:3]
  plot_decay(as.numeric(unlist(my_updated_param_U5)), df_for_decay_U5)
  names(my_updated_param_U5)=c("updated_halflife_U5", "updated_kappa_U5", "use_U5")

  my_updated_param_others=calculate_durability_param(df_for_decay_others)[1:3]
  plot_decay(as.numeric(unlist(my_updated_param_others)), df_for_decay_others)
  names(my_updated_param_others)=c("updated_halflife_others", "updated_kappa_others", "use_others")

  return(c(my_updated_param_U5,my_updated_param_others))
}


#' Continuous Distribution Coverage (Uniform, Under-1)
#'
#' Computes average coverage at year 12 for a continuous distribution scenario
#' assuming uniform age structure.
#'
#' @inheritParams getCov_MassCampaign
#' @param popU1 Numeric. Proportion under 1 year (default = 0.035).
#'
#' @return Numeric. Average usage coverage.
#' @export
getCov_ContinuousU1unif=function(halflife, reach, max_usage, popU1=0.035, shape_weibull=2){

  scenario=generate_continuousDistr_uniformAllAges(halflife_weibull = halflife,
                                                   shape_weibull=shape_weibull,
                                                   reach = reach,
                                                   ageGroupProp=popU1, max_usage = max_usage)
  summary_use=get_average_coverage_continuous(scenario, year=12)

  return(as.numeric(summary_use$use_total))
}


#' Continuous Distribution Coverage (Uniform, School-Age)
#'
#' Computes average coverage for school-age children with adjustment for
#' proportion attending school.
#'
#' @inheritParams getCov_MassCampaign
#' @param pop6to14 Numeric. Proportion aged 6–14 (default = 0.27).
#' @param prop_school_classes Numeric. Proportion attending school (default = 0.5).
#'
#' @return Numeric. Average usage coverage.
#' @export
getCov_ContinuousSchoolunif=function(halflife, reach, max_usage, pop6to14=0.27, prop_school_classes=0.5, shape_weibull=2){

  scenario=generate_continuousDistr_uniformAllAges(halflife_weibull = halflife,
                                                   shape_weibull=shape_weibull,
                                                   reach = reach*prop_school_classes,
                                                   ageGroupProp=pop6to14, max_usage = max_usage)
  summary_use=get_average_coverage_continuous(scenario, year=12)

  return(as.numeric(summary_use$use_total))
}



#' Continuous Distribution Coverage (Age-Specific: Under-1)
#'
#' Computes coverage separately for under-1 and other age groups.
#'
#' @inheritParams getCov_ContinuousU1unif
#'
#' @return List with coverage for under-1 and others.
#' @export
getCov_ContinuousU1=function(halflife, reach, max_usage, popU1=0.035, shape_weibull=2){

  scenario=generate_continuousDistr_SpecificAges(halflife_weibull = halflife,
                                                 shape_weibull=shape_weibull,
                                                 reach = reach,ageGroup="U1",otherageGroup="others",
                                                 ageGroupProp=popU1, max_usage = max_usage)
  summary_use_U1=get_average_coverage_continuous(scenario, year=12,col = "use_U1")
  summary_use_others=get_average_coverage_continuous(scenario, year=12,col = "use_others")

  return(list("use_U1"=as.numeric(summary_use_U1$use_U1),
                     "use_others"=as.numeric(summary_use_others$use_others)))
}



#' Continuous Distribution Coverage (Age-Specific: Schoolchildren)
#'
#' Computes coverage for schoolchildren and others in a continuous distribution.
#'
#' @inheritParams getCov_ContinuousSchoolunif
#'
#' @return List with coverage for schoolchildren and others.
#' @export
getCov_ContinuousSchool=function(halflife, reach, max_usage, pop6to14=0.27, prop_school_classes=0.5, shape_weibull=2){

  scenario=generate_continuousDistr_SpecificAges(halflife_weibull = halflife,
                                                 shape_weibull=shape_weibull,
                                                 reach = reach*prop_school_classes,ageGroup="schoolchildren",otherageGroup="others",
                                                 ageGroupProp=pop6to14, max_usage = max_usage)
  summary_use_schoolchildren=get_average_coverage_continuous(scenario, year=12,col = "use_schoolchildren")
  summary_use_others=get_average_coverage_continuous(scenario, year=12,col = "use_others")

  return(list("use_schoolchildren"=as.numeric(summary_use_schoolchildren$use_schoolchildren),
                     "use_others"=as.numeric(summary_use_others$use_others)))
}



#' Average Coverage from Mass Campaign
#'
#' Computes average coverage between years 12 and 15 for a mass campaign.
#'
#' @inheritParams getCov_MassCampaign
#'
#' @return Numeric. Average usage coverage.
#' @export
getAverageCov_MassCampaign=function(halflife, reach, max_usage, shape_weibull=2){
  scenario=generate_massCampaign_UniformAllAges(halflife_weibull = halflife,
                                                shape_weibull=shape_weibull,
                                                reach=reach, max_usage = max_usage,
                                                frequency=3)
  summary_use=get_average_coverage_mass(scenario, year_start=12, year_end=15)

  return(as.numeric(summary_use$use_total))
}

#' Average Coverage from Mass Campaign (Uniform, Under-5 Proportion)
#'
#' Simulates a mass campaign with uniform coverage across all ages while specifying
#' the proportion of the population under 5, and computes the average coverage
#' between years 12 and 15.
#'
#' @param halflife Numeric. Initial half-life of the intervention.
#' @param reach Numeric. Proportion of the population reached.
#' @param max_usage Numeric. Maximum achievable usage.
#' @param popU5 Numeric. Proportion of population under 5 (default = 0.16).
#' @param shape_weibull Numeric. Shape parameter of the Weibull decay (default = 2).
#'
#' @return Numeric. Average usage coverage over years 12–15.
#' @export
getAverageCov_MassCampaignU5unif=function(halflife, reach, max_usage,popU5=0.16, shape_weibull=2){

  scenario=generate_massCampaign_UniformAllAges(halflife_weibull = halflife,
                                                shape_weibull=shape_weibull,
                                                reach=reach, max_usage = max_usage,
                                                ageGroupProp=popU5,
                                                frequency=3)
  summary_use=get_average_coverage_mass(scenario, year_start=12, year_end=15)

  return(as.numeric(summary_use$use_total))
}

#' Average Coverage from Mass Campaign (Age-Specific: Under-5 vs Others)
#'
#' Simulates a mass campaign with separate coverage for under-5 and other age groups,
#' and computes average coverage between years 12 and 15 for each group.
#'
#' @inheritParams getAverageCov_MassCampaignU5unif
#'
#' @return List with:
#' \describe{
#'   \item{use_U5}{Numeric. Average coverage for under-5 population.}
#'   \item{use_others}{Numeric. Average coverage for the remaining population.}
#' }
#' @export
getAverageCov_MassCampaignU5=function(halflife, reach, max_usage, popU5=0.16, shape_weibull=2){

  scenario=generate_massCampaign_SpecificAges(halflife_weibull = halflife,
                                              shape_weibull=shape_weibull,
                                              ageGroup="U5",otherageGroup="others",
                                              ageGroupProp=popU5, max_usage = max_usage,
                                              reach = reach, frequency=3)

  summary_use_U5=get_average_coverage_mass(scenario, year_start=12, year_end=15, col="use_U5")
  summary_use_others=get_average_coverage_mass(scenario, year_start=12, year_end=15, col="use_others")

  return(list("use_U5"=summary_use_U5$use_U5,
              "use_others"=summary_use_others$use_others))
}

