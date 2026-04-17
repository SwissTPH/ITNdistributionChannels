#' Continuous Distribution Coverage, school in non campaign years (Age-Specific: Schoolchildren)
#'
#' Computes coverage for schoolchildren and others in a continuous distribution.
#'
#' @inheritParams getCov_ContinuousSchoolunif
#' @frequency_mass Frequency of the mass distribution
#'
#' @return List with coverage for schoolchildren and others.
#' @export
getCov_AlternateSchool=function(halflife, reach, max_usage, pop6to14=0.27, prop_school_classes=0.5, shape_weibull=2,
                                frequency_mass=3, year_start=12){

  scenario=generate_AlternateCont_SpecificAges(halflife_weibull = halflife,
                                                 shape_weibull=shape_weibull,frequency_mass = frequency_mass,
                                                 reach_cont = reach*prop_school_classes,ageGroup_cont ="schoolchildren",otherageGroup_cont = "others",
                                                 ageGroupProp_cont =pop6to14, max_usage = max_usage)%>%
    mutate(reach=reach, max_usage=max_usage, halflife=halflife)


  summary_use_schoolchildren=get_average_coverage_continuous(scenario, year=year_start,col = "use_schoolchildren")
  summary_use_others=get_average_coverage_continuous(scenario, year=year_start,col = "use_others")


  return(list("use_schoolchildren"=as.numeric(summary_use_schoolchildren$use_schoolchildren),
              "use_others"=as.numeric(summary_use_others$use_others)))
}




#' Coverage from Mass Campaign every 5 years (All Ages, Decay Fit)
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
getCov_MassCampaign5years=function(halflife, reach, max_usage, shape_weibull=2){
  scenario=generate_massCampaign_UniformAllAges(halflife_weibull = halflife,
                                                shape_weibull=shape_weibull,
                                                reach=reach, max_usage = max_usage,
                                                frequency=5)
  scenario2=scenario |>
    dplyr::filter(time>=10*365+1, time<15*365+1)

  df_for_decay=data.frame(value= scenario2$use_total, time=(1:length(scenario2$use_total))/365)
  my_updated_param=calculate_durability_param(df_for_decay)[1:3]
  plot_decay(as.numeric(unlist(my_updated_param)), df_for_decay)

  names(my_updated_param)=c("updated_halflife", "updated_kappa", "use_total")

  return(my_updated_param)
}



#' Continuous Distribution Coverage (Age-Specific: Schoolchildren + Routine)
#'
#' Computes coverage for schoolchildren and others in a continuous distribution.
#'
#' @inheritParams getCov_ContinuousSchoolunif
#'
#' @return List with coverage for schoolchildren and others.
#' @export
getCov_ContinuousSchoolRoutine=function(halflife, reach_school, reach_routine, max_usage, pop6to14=0.27,prop_school_classes=0.5, pop_routine=0.035, shape_weibull=2){

  scenario=generate_ContAndContinuousDistr_SpecificAges(halflife_weibull = halflife,shape_weibull = shape_weibull,
                                               reach_cont1 = reach_routine, reach_cont2 = reach_school*prop_school_classes,
                                               ageGroupProp_cont1 = pop_routine,ageGroupProp_cont2 = pop6to14,
                                               ageGroupCont1 ="U1",ageGroupCont2 = "schoolchildren",
                                               max_usage = max_usage)
  summary_use_schoolchildren=get_average_coverage_continuous(scenario %>% mutate(reach=reach_school), year=12,col = "use_schoolchildren")
  summary_use_others=get_average_coverage_continuous(scenario %>% mutate(reach=reach_school), year=12,col = "use_others")
  summary_use_U1=get_average_coverage_continuous(scenario %>% mutate(reach=reach_routine), year=12,col = "use_U1")

  return(list("use_schoolchildren"=as.numeric(summary_use_schoolchildren$use_schoolchildren),
              "use_U1"=as.numeric(summary_use_U1$use_U1),
              "use_others"=as.numeric(summary_use_others$use_others)))
}

