getCov_MassCampaign=function(halflife, reach, max_usage){
  scenario=generate_massCampaign_UniformAllAges(halflife_weibull = halflife,
                                                reach=reach, max_usage = max_usage,
                                                frequency=3)
  summary_use=get_start_coverage_mass(scenario, year=12)

  return(as.numeric(summary_use$use_total))
}


getCov_MassCampaignU5unif=function(halflife, reach, max_usage,popU5=0.16){

  scenario=generate_massCampaign_UniformAllAges(halflife_weibull = halflife,
                                                reach=reach, max_usage = max_usage,
                                                ageGroupProp=popU5,
                                                frequency=3)
  summary_use=get_start_coverage_mass(scenario, year=12)

  return(as.numeric(summary_use$use_total))
}

getCov_MassCampaignU5=function(halflife, reach, max_usage, popU5=0.16){

  scenario=generate_massCampaign_SpecificAges(halflife_weibull = halflife,
                                              ageGroup="U5",otherageGroup="others",
                                              ageGroupProp=popU5, max_usage = max_usage,
                                              reach = reach, frequency=3)

  summary_use=get_start_coverage_mass(scenario, year=12)

  return(list("use_U5"=as.numeric(summary_use$use_U5),
              "use_others"=as.numeric(summary_use$use_others)))
}


getCov_ContinuousU1unif=function(halflife, reach, max_usage, popU1=0.035){

  scenario=generate_continuousDistr_uniformAllAges(halflife_weibull = halflife,
                                                   reach = reach,
                                                   ageGroupProp=popU1, max_usage = max_usage)
  summary_use=get_average_coverage_continuous(scenario, year=12)

  return(as.numeric(summary_use$use_total))
}

getCov_ContinuousSchoolunif=function(halflife, reach, max_usage, pop6to14=0.27, prop_school_classes=0.5){

  scenario=generate_continuousDistr_uniformAllAges(halflife_weibull = halflife,
                                                   reach = reach*prop_school_classes,
                                                   ageGroupProp=pop6to14, max_usage = max_usage)
  summary_use=get_average_coverage_continuous(scenario, year=12)

  return(as.numeric(summary_use$use_total))
}



getCov_ContinuousU1=function(halflife, reach, max_usage, popU1=0.035){

  scenario=generate_continuousDistr_SpecificAges(halflife_weibull = halflife,
                                                   reach = reach,ageGroup="U1",otherageGroup="others",
                                                 ageGroupProp=popU1, max_usage = max_usage)
  summary_use_U1=get_average_coverage_continuous(scenario, year=12,col = "use_U1")
  summary_use_others=get_average_coverage_continuous(scenario, year=12,col = "use_others")

  return(list("use_U1"=as.numeric(summary_use_U1$use_U1),
                     "use_others"=as.numeric(summary_use_others$use_others)))
}




getCov_ContinuousSchool=function(halflife, reach, max_usage, pop6to14=0.27, prop_school_classes=0.5){

  scenario=generate_continuousDistr_SpecificAges(halflife_weibull = halflife,
                                                 reach = reach*prop_school_classes,ageGroup="schoolchildren",otherageGroup="others",
                                                 ageGroupProp=pop6to14, max_usage = max_usage)
  summary_use_schoolchildren=get_average_coverage_continuous(scenario, year=12,col = "use_schoolchildren")
  summary_use_others=get_average_coverage_continuous(scenario, year=12,col = "use_others")

  return(list("use_schoolchildren"=as.numeric(summary_use_schoolchildren$use_schoolchildren),
                     "use_others"=as.numeric(summary_use_others$use_others)))
}
