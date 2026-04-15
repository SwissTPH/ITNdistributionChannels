#' Mass Campaign + Continuous ITN Distribution (Uniform Across All Ages)
#'
#' Simulates periodic mass ITN campaigns to the whole population combined with a continuous distribution. Coverage is uniformly applied
#' across the entire population.
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_mass Numeric. Proportion of the population reached by the mass campaign.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export
generate_MassAndContinuousDistr_uniformAllAges=function(halflife_weibull=2.1,
                                                        shape_weibull=2,max_usage=0.8,
                                                        reach_mass=0.8, frequency_mass=3,
                                                        reach_cont=0.8, ageGroupProp_cont=0.035){
  mass=generate_massCampaign_UniformAllAges(halflife_weibull=halflife_weibull,
                                            shape_weibull=shape_weibull,
                                            reach=reach_mass, max_usage=max_usage,
                                            frequency=frequency_mass, ageGroupProp=1)
  cont=generate_continuousDistr_uniformAllAges(halflife_weibull=halflife_weibull,
                                              shape_weibull=shape_weibull,
                                              reach=reach_cont, max_usage=max_usage,
                                              ageGroupProp=ageGroupProp_cont)

  combined=merge(mass |> dplyr::select(-use_total) |>  dplyr::rename(cov_mass=cov, reach_mass=reach),
                 cont |> dplyr::select(-use_total) |>  dplyr::rename(cov_cont=cov, reach_cont=reach)) |>
    dplyr::mutate(cov=cov_mass+cov_cont,
                  use_total=pmin(cov*max_usage, max_usage))|>
    dplyr::select(-cov_cont, -cov_mass)
  return(combined)
}

#' Mass Campaign + Continuous ITN Distribution (Specific Age Groups)
#'
#' Simulates periodic mass ITN campaigns to the whole population combined with a continuous distribution. Coverage is provided by specific age group
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_mass Numeric. Proportion of the population reached by the mass campaign.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#' @param ageGroup_cont Character. Name of the targeted age group for the continuous distribution.
#' @param otherageGroup_cont Character. Name of the non-targeted age group for the continuous distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export
generate_MassAndContinuousDistr_SpecificAges=function(halflife_weibull=2.1,
                                                        shape_weibull=2,max_usage=0.8,
                                                        reach_mass=0.8, frequency_mass=3,
                                                        reach_cont=0.8, ageGroupProp_cont=0.035,
                                                        ageGroup_cont="U1",otherageGroup_cont="others"){
  mass=generate_massCampaign_UniformAllAges(halflife_weibull=halflife_weibull,
                                          shape_weibull=shape_weibull,
                                          reach=reach_mass, max_usage=max_usage,
                                          frequency=frequency_mass, ageGroupProp=1)

  cont=generate_continuousDistr_SpecificAges(halflife_weibull=halflife_weibull,
                                             shape_weibull=shape_weibull,
                                             reach=reach_cont, max_usage=max_usage,
                                             ageGroupProp=ageGroupProp_cont,
                                           ageGroup=ageGroup_cont,otherageGroup=otherageGroup_cont, full_output = T)

  combined=merge(mass |> dplyr::select(-use_total) |>  dplyr::rename(cov_mass=cov, reach_mass=reach),
               cont |> dplyr::select(-use_targetpop, -use_otherpop)) |>
  dplyr::mutate(cov_targetgroup=cov_mass+totcov_targetgroup,
                spillover2_cov_othergroup=pmax(0, cov_targetgroup-1)*ageGroupProp_cont/(1-ageGroupProp_cont),
                cov_othergroup=cov_mass+totcov_othergroup+spillover2_cov_othergroup,
                use_targetpop=pmin(cov_targetgroup*max_usage, max_usage),
                use_otherpop=pmin(cov_othergroup*max_usage, max_usage)) |>
    dplyr::select(time, use_targetpop, use_otherpop)|>
    dplyr::mutate( reach_mass=reach_mass, reach_cont=reach_cont, max_usage=max_usage , halflife=halflife_weibull)

  names(combined)[c(2,3)]=c(paste0("use_", ageGroup_cont),paste0("use_", otherageGroup_cont))
  return(combined)
}



#' Continuous + Continuous ITN Distribution (Uniform Across All Ages)
#'
#' Simulates combined continuous distribution to different target age groups. Coverage is uniformly applied
#' across the entire population.
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_cont1 Numeric. Proportion of the population reached by the first continuous distribution.
#' @param reach_cont2 Numeric. Proportion of the population reached by the second continuous distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param ageGroupProp_cont1 Numeric. Proportion of population in the targeted age group for the first continuous distribution.
#' @param ageGroupProp_cont2 Numeric. Proportion of population in the targeted age group for the second continuous distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export
generate_ContAndContinuousDistr_uniformAllAges=function(halflife_weibull=2.1,
                                                        shape_weibull=2,max_usage=0.8,
                                                        reach_cont1=0.8, ageGroupProp_cont1=0.035,
                                                        reach_cont2=0.8*0.5, ageGroupProp_cont2=0.27){

  cont1=generate_continuousDistr_uniformAllAges(halflife_weibull=halflife_weibull,
                                               shape_weibull=shape_weibull,
                                               reach=reach_cont1, max_usage=max_usage,
                                               ageGroupProp=ageGroupProp_cont1)
  cont2=generate_continuousDistr_uniformAllAges(halflife_weibull=halflife_weibull,
                                               shape_weibull=shape_weibull,
                                               reach=reach_cont2, max_usage=max_usage,
                                               ageGroupProp=ageGroupProp_cont2)

  combined=merge(cont1 |> dplyr::select(-use_total) |>  dplyr::rename(cov_cont1=cov, reach_cont1=reach),
                 cont2 |> dplyr::select(-use_total) |>  dplyr::rename(cov_cont2=cov, reach_cont2=reach)) |>
    dplyr::mutate(cov=cov_cont1+cov_cont2,
                  use_total=pmin(cov*max_usage, max_usage))|>
    dplyr::select(-cov_cont1, -cov_cont2)
  return(combined)
}

#' Continuous + Continuous ITN Distribution (Specific Age Groups)
#'
#' Simulates combined continuous distribution to different target age groups. Coverage is provided by specific age group
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_cont1 Numeric. Proportion of the population reached by the first continuous distribution.
#' @param reach_cont2 Numeric. Proportion of the population reached by the second continuous distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont1 Numeric. Proportion of population in the targeted age group for the first continuous distribution.
#' @param ageGroupProp_cont2 Numeric. Proportion of population in the targeted age group for the second continuous distribution.
#' @param ageGroup_cont1 Character. Name of the targeted age group for the first continuous distribution.
#' @param otherageGroup_cont1 Character. Name of the non-targeted age group for the first continuous distribution.
#' @param ageGroup_cont2 Character. Name of the targeted age group for the second continuous distribution.
#' @param otherageGroup_cont2 Character. Name of the non-targeted age group for the second continuous distribution.
#' @param overlap Boolean. FALSE (default) indicates that the age groups for the two distribution channels don't overlap. TRUE is not yet supported.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export
generate_ContAndContinuousDistr_SpecificAges=function(halflife_weibull=2.1,
                                                      shape_weibull=2,max_usage=0.8,
                                                      reach_cont1=0.8, ageGroupProp_cont1=0.035,
                                                      reach_cont2=0.8*0.5, ageGroupProp_cont2=0.27,
                                                      ageGroupCont1="U1",otherageGroupCont1="others",
                                                      ageGroupCont2="school",otherageGroupCont2="others", overlap=FALSE){
  if(overlap){
    stop("overlap =TRUE not implemented yet")
  } else {
    cont1=generate_continuousDistr_SpecificAges(halflife_weibull=halflife_weibull,
                                                shape_weibull=shape_weibull,
                                                reach=reach_cont1, max_usage=max_usage,
                                                ageGroupProp=ageGroupProp_cont1,
                                                ageGroup=ageGroupCont1,otherageGroup=paste0(otherageGroupCont1,1), full_output = T)|>
      dplyr::select(-use_targetpop, -use_otherpop, -spillover_cov_othergroup)

    names(cont1)=gsub("group", "group1",names(cont1))

    cont2=generate_continuousDistr_SpecificAges(halflife_weibull=halflife_weibull,
                                                shape_weibull=shape_weibull,
                                                reach=reach_cont2, max_usage=max_usage,
                                                ageGroupProp=ageGroupProp_cont2,
                                                ageGroup=ageGroupCont2,otherageGroup=paste0(otherageGroupCont2,2), full_output = T)|>
      dplyr::select(-use_targetpop, -use_otherpop, -spillover_cov_othergroup)


    names(cont2)=gsub("group", "group2",names(cont2))


    combined=merge(cont1,cont2) |>
      dplyr::mutate(cov_targetgroup1_0=totcov_targetgroup1+totcov_othergroup2,
                    cov_targetgroup2_0=totcov_targetgroup2+totcov_othergroup1,
                    group1maxed=(cov_targetgroup1_0>=1),
                    group2maxed=(cov_targetgroup2_0>=1),
                    denominator_spillover_othergroup1=ifelse(group2maxed,(1-ageGroupProp_cont1-ageGroupProp_cont2),1-ageGroupProp_cont1),
                    denominator_spillover_othergroup2=ifelse(group1maxed,(1-ageGroupProp_cont1-ageGroupProp_cont2),1-ageGroupProp_cont2),
                    spillover2_cov_othergroup1=pmax(0, cov_targetgroup1_0-1)*ageGroupProp_cont1/denominator_spillover_othergroup1,
                    spillover2_cov_othergroup2=pmax(0, cov_targetgroup2_0-1)*ageGroupProp_cont2/denominator_spillover_othergroup2,
                    cov_othergroup=totcov_othergroup1+totcov_othergroup2+spillover2_cov_othergroup1+spillover2_cov_othergroup2,
                    cov_targetgroup1=cov_targetgroup1_0+ifelse(group1maxed,0,spillover2_cov_othergroup2),
                    cov_targetgroup2=cov_targetgroup2_0+ifelse(group2maxed,0,spillover2_cov_othergroup1),
                    use_targetpop1=pmin(cov_targetgroup1*max_usage, max_usage),
                    use_targetpop2=pmin(cov_targetgroup2*max_usage, max_usage),
                    use_otherpop=pmin(cov_othergroup*max_usage, max_usage)) |>
      dplyr::select(time, use_targetpop1,use_targetpop2, use_otherpop)|>
      dplyr::mutate( reach_cont1=reach_cont1, reach_cont2=reach_cont2, max_usage=max_usage , halflife=halflife_weibull)

    names(combined)[c(2,3,4)]=c(paste0("use_", ageGroupCont1),paste0("use_", ageGroupCont2),paste0("use_others"))
  }

  return(combined)
}

#' Continuous ITN Distribution in alternate years (Uniform Across All Ages)
#'
#' Simulates a continuous distribution during the years without mass campaign. Coverage is uniformly applied
#' across the entire population.
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#' @param rename Boolean. If TRUE, final dataset includes ITN use calculations
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export

generate_AlternateCont_uniformAllAges=function(halflife_weibull=2.1,
                                               shape_weibull=2,max_usage=0.8,
                                               frequency_mass=3,
                                               reach_cont=0.8*0.5, ageGroupProp_cont=0.27, rename=TRUE){

  cont=generate_massCampaign_UniformAllAges(halflife_weibull=halflife_weibull,
                                            shape_weibull=shape_weibull,frequency = frequency_mass,
                                            reach=reach_cont, max_usage=max_usage,
                                            ageGroupProp=ageGroupProp_cont) |>
    dplyr::select(-use_total) |>  dplyr::rename(reach_cont=reach)

  cont_tot = cont |> dplyr::rename(cov_cont=cov) |> dplyr::mutate(time=time+365)
  for(i in 2:(frequency_mass-1)){
    cont_tot=merge(cont_tot,
                   cont |> dplyr::rename(cov_contY2=cov) |> dplyr::mutate(time=time+i*365),
                   all=T) |>
      dplyr::mutate(cov_cont=ifelse(is.na(cov_cont),0,cov_cont)+ifelse(is.na(cov_contY2),0,cov_contY2)) |>
      dplyr::select(time, halflife, cov_cont)
  }

  cont_tot_final=cont_tot|>  dplyr::right_join(cont |>  dplyr::select(time, halflife))|> dplyr::mutate(cov_cont=ifelse(is.na(cov_cont),0, cov_cont))|> dplyr::arrange(time)

  if(rename){
    cont_tot_final=cont_tot_final|>
      dplyr::mutate( use_total=pmin(cov_cont*max_usage, max_usage))|> dplyr::rename(cov=cov_cont)
    }

  return(cont_tot_final)
}


#'  Continuous ITN Distribution in alternate years (Specific Age Groups)
#'
#' Simulates a continuous distribution during the years without mass campaign. Coverage is provided by specific age group
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#' @param ageGroup_cont Character. Name of the targeted age group for the continuous distribution.
#' @param otherageGroup_cont Character. Name of the non-targeted age group for the continuous distribution.
#' @param rename Boolean. If TRUE, final dataset is renamed and includes ITN use calculations
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export

generate_AlternateCont_SpecificAges=function(halflife_weibull=2.1,
                                             shape_weibull=2,max_usage=0.8,
                                             frequency_mass=3,
                                             reach_cont=0.8*0.5, ageGroupProp_cont=0.27,
                                             ageGroup_cont="school",otherageGroup_cont="others", rename=TRUE){

  cont=generate_massCampaign_SpecificAges(halflife_weibull=halflife_weibull,
                                          shape_weibull=shape_weibull,
                                          reach=reach_cont, max_usage=max_usage,
                                          frequency=frequency_mass,
                                          ageGroupProp=ageGroupProp_cont,
                                          ageGroup=ageGroup_cont,otherageGroup=otherageGroup_cont, full_output = T)|>
    dplyr::select(time, totcov_targetgroup,totcov_othergroup)

  cont_tot = cont  |> dplyr::mutate(time=time+365)
  for(i in 2:(frequency_mass-1)){
    cont_tot=merge(cont_tot,
                   cont |> dplyr::rename(cov_target_contY2=totcov_targetgroup, cov_other_contY2=totcov_othergroup) |>
                     dplyr:: mutate(time=time+i*365),
                   all=T) |>
      dplyr::mutate(totcov_targetgroup=ifelse(is.na(totcov_targetgroup),0,totcov_targetgroup)+ifelse(is.na(cov_target_contY2),0,cov_target_contY2),
             totcov_othergroup=ifelse(is.na(totcov_othergroup),0,totcov_othergroup)+ifelse(is.na(cov_other_contY2),0,cov_other_contY2)) |>
      dplyr::select(time, totcov_targetgroup, totcov_othergroup)
  }
  cont_tot_final=cont_tot|>  dplyr::right_join(cont |>  dplyr::select(time))|>
    dplyr::mutate(totcov_targetgroup=ifelse(is.na(totcov_targetgroup),0, totcov_targetgroup),
           totcov_othergroup=ifelse(is.na(totcov_othergroup),0, totcov_othergroup))|>
    dplyr::arrange(time)

  if(rename){
    cont_tot_final=cont_tot_final|>
      dplyr::mutate(spillover_cov_othergroup=pmax(0, totcov_targetgroup-1)*ageGroupProp_cont/(1-ageGroupProp_cont),
                    use_targetpop=pmin(totcov_targetgroup*max_usage, max_usage),
                    use_otherpop=pmin((totcov_othergroup+spillover_cov_othergroup)*max_usage, max_usage)) |>
    dplyr::select(time, use_targetpop, use_otherpop)
    names(cont_tot_final)[c(2,3)]=c(paste0("use_", ageGroup_cont),paste0("use_", otherageGroup_cont))}
  return(cont_tot_final)
}





#' Mass Campaign + Continuous ITN Distribution in alternate years (Uniform Across All Ages)
#'
#' Simulates periodic mass ITN campaigns to the whole population combined with a continuous distribution during the years without mass campaign. Coverage is uniformly applied
#' across the entire population.
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_mass Numeric. Proportion of the population reached by the mass campaign.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export
generate_AlternateMassCont_uniformAllAges=function(halflife_weibull=2.1,
                                                    shape_weibull=2,max_usage=0.8,
                                                    reach_mass=0.8, frequency_mass=3,
                                                    reach_cont=0.8*0.5, ageGroupProp_cont=0.27){

  mass=generate_massCampaign_UniformAllAges(halflife_weibull=halflife_weibull,
                                            shape_weibull=shape_weibull,
                                            reach=reach_mass, max_usage=max_usage,
                                            frequency=frequency_mass, ageGroupProp=1)

  cont=generate_AlternateCont_uniformAllAges(halflife_weibull=halflife_weibull,
                                             shape_weibull=shape_weibull,frequency = frequency_mass,
                                             reach=reach_cont, max_usage=max_usage,
                                             ageGroupProp=ageGroupProp_cont, rename=FALSE)

  combined=merge(mass |> dplyr::select(-use_total) |>  dplyr::rename(cov_mass=cov, reach_mass=reach),
                 cont, all.x=T) |>
    dplyr::mutate(cov=cov_mass+ifelse(is.na(cov_cont),0,cov_cont),
                  use_total=pmin(cov*max_usage, max_usage))|>
    dplyr::select(-cov_cont, -cov_mass)|>
    dplyr::mutate(reach_cont=reach_cont)|> dplyr::arrange(time)

  return(combined)
}


#' Mass Campaign + Continuous ITN Distribution in alternate years (Specific Age Groups)
#'
#' Simulates periodic mass ITN campaigns to the whole population combined with a continuous distribution during the years without mass campaign. Coverage is provided by specific age group
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_mass Numeric. Proportion of the population reached by the mass campaign.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#' @param ageGroup_cont Character. Name of the targeted age group for the continuous distribution.
#' @param otherageGroup_cont Character. Name of the non-targeted age group for the continuous distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export
generate_AlternateMassCont_SpecificAges=function(halflife_weibull=2.1,
                                                  shape_weibull=2,max_usage=0.8,
                                                  reach_mass=0.8, frequency_mass=3,
                                                  reach_cont=0.8, ageGroupProp_cont=0.035,
                                                  ageGroup_cont="school",otherageGroup_cont="others"){

  mass=generate_massCampaign_UniformAllAges(halflife_weibull=halflife_weibull,
                                            shape_weibull=shape_weibull,
                                            reach=reach_mass, max_usage=max_usage,
                                            frequency=frequency_mass, ageGroupProp=1)

  cont=generate_AlternateCont_SpecificAges(halflife_weibull=halflife_weibull,
                                           shape_weibull=shape_weibull,
                                           reach=reach_cont, max_usage=max_usage,
                                           frequency=frequency_mass,
                                           ageGroupProp_cont =ageGroupProp_cont,
                                           ageGroup_cont=ageGroup_cont,otherageGroup_cont=otherageGroup_cont, rename = FALSE)

  combined=merge(mass |> dplyr::select(-use_total) |>  dplyr::rename(cov_mass=cov, reach_mass=reach),
                 cont, all.x = T) |>
    dplyr::mutate(cov_targetgroup=cov_mass+ifelse(is.na(totcov_targetgroup),0,totcov_targetgroup),
                  spillover2_cov_othergroup=pmax(0, cov_targetgroup-1)*ageGroupProp_cont/(1-ageGroupProp_cont),
                  cov_othergroup=cov_mass+ifelse(is.na(totcov_othergroup),0,totcov_othergroup)+spillover2_cov_othergroup,
                  use_targetpop=pmin(cov_targetgroup*max_usage, max_usage),
                  use_otherpop=pmin(cov_othergroup*max_usage, max_usage)) |>
    dplyr::select(time, use_targetpop, use_otherpop)|>
    dplyr::mutate( reach_mass=reach_mass, reach_cont=reach_cont, max_usage=max_usage , halflife=halflife_weibull)

  names(combined)[c(2,3)]=c(paste0("use_", ageGroup_cont),paste0("use_", otherageGroup_cont))
  return(combined)
}



#' Mass Campaign + Continuous ITN Distribution in alternate years + routine distribution (Uniform Across All Ages)
#'
#' Simulates periodic mass ITN campaigns to the whole population combined with a continuous distribution during the years without mass campaign and a routine distribution every year. Coverage is uniformly applied
#' across the entire population.
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_mass Numeric. Proportion of the population reached by the mass campaign.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param reach_routine Numeric. Proportion of the population reached by the routine distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#' @param ageGroupProp_routine Numeric. Proportion of population in the targeted age group for the routine distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export


generate_AlternateMassContRoutine_uniformAllAges=function(halflife_weibull=2.1,
                                                           shape_weibull=2,max_usage=0.8,
                                                           reach_mass=0.8, frequency_mass=3,
                                                           reach_cont=0.8*0.5, ageGroupProp_cont=0.27,
                                                           reach_routine=0.8, ageGroupProp_routine=0.035){

  mass=generate_massCampaign_UniformAllAges(halflife_weibull=halflife_weibull,
                                            shape_weibull=shape_weibull,
                                            reach=reach_mass, max_usage=max_usage,
                                            frequency=frequency_mass, ageGroupProp=1)

  cont=generate_AlternateCont_uniformAllAges(halflife_weibull=halflife_weibull,
                                             shape_weibull=shape_weibull,frequency = frequency_mass,
                                             reach=reach_cont, max_usage=max_usage,
                                             ageGroupProp=ageGroupProp_cont, rename=FALSE)

  routine=generate_continuousDistr_uniformAllAges(halflife_weibull=halflife_weibull,
                                                  shape_weibull=shape_weibull,
                                                  reach=reach_routine, max_usage=max_usage,
                                                  ageGroupProp=ageGroupProp_routine)

  combined=merge(mass |> dplyr::select(-use_total) |>  dplyr::rename(cov_mass=cov, reach_mass=reach),
                 merge(cont,
                       routine |> dplyr::select(-use_total) |>  dplyr::rename(cov_routine=cov, reach_routine=reach), all.y=T)) |>
    dplyr::mutate(cov=cov_mass+ifelse(is.na(cov_cont),0,cov_cont)+cov_routine,
                  use_total=pmin(cov*max_usage, max_usage),
                  reach_cont=reach_cont)|>
    dplyr::select(-cov_cont, -cov_mass, -cov_routine)|> dplyr::arrange(time)

  return(combined)
}



#' Mass Campaign + Continuous ITN Distribution in alternate years + routine distribution (Specific Age Groups)
#'
#' Simulates periodic mass ITN campaigns to the whole population combined with a continuous distribution during the years without mass campaign and a routine distribution every year. Coverage is provided by specific age group
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_mass Numeric. Proportion of the population reached by the mass campaign.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param reach_routine Numeric. Proportion of the population reached by the routine distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#' @param ageGroupProp_routine Numeric. Proportion of population in the targeted age group for the routine distribution.
#' @param ageGroup_cont Character. Name of the targeted age group for the continuous distribution.
#' @param otherageGroup_cont Character. Name of the non-targeted age group for the continuous distribution.
#' @param ageGroup_routine Character. Name of the targeted age group for the routine distribution.
#' @param otherageGroup_routine Character. Name of the non-targeted age group for the routine distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export

generate_AlternateMassContRoutine_SpecificAges=function(halflife_weibull=2.1,
                                                         shape_weibull=2,max_usage=0.8,
                                                         reach_mass=0.8, frequency_mass=3,
                                                         reach_cont=0.8*0.5, ageGroupProp_cont=0.27,
                                                         reach_routine=0.8, ageGroupProp_routine=0.035,
                                                         ageGroup_cont="school",otherageGroup_cont="others",
                                                         ageGroup_routine="U1",otherageGroup_routine="others"){

  mass=generate_massCampaign_UniformAllAges(halflife_weibull=halflife_weibull,
                                            shape_weibull=shape_weibull,
                                            reach=reach_mass, max_usage=max_usage,
                                            frequency=frequency_mass, ageGroupProp=1)

  cont=generate_AlternateCont_SpecificAges(halflife_weibull=halflife_weibull,
                                           shape_weibull=shape_weibull,
                                           reach=reach_cont, max_usage=max_usage,
                                           frequency=frequency_mass,
                                           ageGroupProp_cont =ageGroupProp_cont,
                                           ageGroup_cont=ageGroup_cont,otherageGroup_cont=otherageGroup_cont, rename = FALSE)

  routine=generate_continuousDistr_SpecificAges(halflife_weibull=halflife_weibull,
                                                shape_weibull=shape_weibull,
                                                reach=reach_routine, max_usage=max_usage,
                                                ageGroupProp=ageGroupProp_routine,
                                                ageGroup=ageGroup_routine,otherageGroup=otherageGroup_routine, full_output = T)|>
    dplyr::select(-use_targetpop, -use_otherpop, -spillover_cov_othergroup)


  names(routine)=gsub("group", "group_routine",names(routine))


  combined=merge(merge(mass |> dplyr::select(-use_total) |>  dplyr::rename(cov_mass=cov, reach_mass=reach),
                       cont |> dplyr::mutate(reach_cont=reach_cont), all.x = T),
                 routine |> dplyr::mutate(reach_routine=reach_routine))|>
    dplyr::mutate(totcov_targetgroup=ifelse(is.na(totcov_targetgroup),0,totcov_targetgroup),
                  totcov_othergroup=ifelse(is.na(totcov_othergroup),0,totcov_othergroup),
                  cov_targetgroup=cov_mass+totcov_targetgroup+totcov_othergroup_routine,
                  cov_targetgroup_routine=cov_mass+totcov_targetgroup_routine+totcov_othergroup,
                  group1maxed=(cov_targetgroup>=1),
                  group2maxed=(cov_targetgroup_routine>=1),
                  denominator_spillover_othergroup1=ifelse(group2maxed,(1-ageGroupProp_cont-ageGroupProp_routine),1-ageGroupProp_cont),
                  denominator_spillover_othergroup2=ifelse(group1maxed,(1-ageGroupProp_cont-ageGroupProp_routine),1-ageGroupProp_routine),
                  spillover2_cov_othergroup=pmax(0, cov_targetgroup-1)*ageGroupProp_cont/denominator_spillover_othergroup1,
                  spillover2_cov_othergroup_routine=pmax(0, cov_targetgroup_routine-1)*ageGroupProp_routine/denominator_spillover_othergroup2,
                  cov_othergroup=cov_mass+totcov_othergroup+totcov_othergroup_routine+spillover2_cov_othergroup+spillover2_cov_othergroup_routine,
                  cov_targetgroup=cov_targetgroup+ifelse(group1maxed,0,spillover2_cov_othergroup_routine),
                  cov_targetgroup_routine=cov_targetgroup_routine+ifelse(group2maxed,0,spillover2_cov_othergroup),
                  use_targetpop=pmin(cov_targetgroup*max_usage, max_usage),
                  use_targetpop_routine=pmin(cov_targetgroup_routine*max_usage, max_usage),
                  use_otherpop=pmin(cov_othergroup*max_usage, max_usage)) |>
    dplyr::select(time, use_targetpop,use_targetpop_routine, use_otherpop)|>
    dplyr::mutate( reach_mass=reach_mass, reach_cont=reach_cont, reach_routine=reach_routine, max_usage=max_usage , halflife=halflife_weibull)

  names(combined)[c(2,3,4)]=c(paste0("use_", ageGroup_cont),paste0("use_", ageGroup_routine), "use_others")
  return(combined)
}





#' Continuous ITN Distribution in alternate years + routine distribution (Uniform Across All Ages)
#'
#' Simulates a continuous distribution during the years without mass campaign and a routine distribution every year. Coverage is uniformly applied
#' across the entire population.
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param reach_routine Numeric. Proportion of the population reached by the routine distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#' @param ageGroupProp_routine Numeric. Proportion of population in the targeted age group for the routine distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export


generate_AlternateContRoutine_uniformAllAges=function(halflife_weibull=2.1,
                                                          shape_weibull=2,max_usage=0.8,
                                                          frequency_mass=3,
                                                          reach_cont=0.8*0.5, ageGroupProp_cont=0.27,
                                                          reach_routine=0.8, ageGroupProp_routine=0.035){

  cont=generate_AlternateCont_uniformAllAges(halflife_weibull=halflife_weibull,
                                             shape_weibull=shape_weibull,frequency = frequency_mass,
                                             reach=reach_cont, max_usage=max_usage,
                                             ageGroupProp=ageGroupProp_cont, rename=FALSE)

  routine=generate_continuousDistr_uniformAllAges(halflife_weibull=halflife_weibull,
                                                  shape_weibull=shape_weibull,
                                                  reach=reach_routine, max_usage=max_usage,
                                                  ageGroupProp=ageGroupProp_routine)

  combined=merge(cont,
                 routine |> dplyr::select(-use_total) |>  dplyr::rename(cov_routine=cov, reach_routine=reach), all.y=T) |>
    dplyr::mutate(cov=ifelse(is.na(cov_cont),0,cov_cont)+cov_routine,
                  use_total=pmin(cov*max_usage, max_usage),
                  reach_cont=reach_cont)|>
    dplyr::select(-cov_cont, -cov_routine)|> dplyr::arrange(time)

  return(combined)
}



#' Continuous ITN Distribution in alternate years + routine distribution (Specific Age Groups)
#'
#' Simulates a continuous distribution during the years without mass campaign and a routine distribution every year. Coverage is provided by specific age group
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach_cont Numeric. Proportion of the population reached by the continuous distribution.
#' @param reach_routine Numeric. Proportion of the population reached by the routine distribution.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency_mass Numeric. Distribution frequency for the mass campaign(years).
#' @param ageGroupProp_cont Numeric. Proportion of population in the targeted age group for the continuous distribution.
#' @param ageGroupProp_routine Numeric. Proportion of population in the targeted age group for the routine distribution.
#' @param ageGroup_cont Character. Name of the targeted age group for the continuous distribution.
#' @param otherageGroup_cont Character. Name of the non-targeted age group for the continuous distribution.
#' @param ageGroup_routine Character. Name of the targeted age group for the routine distribution.
#' @param otherageGroup_routine Character. Name of the non-targeted age group for the routine distribution.
#'
#' @return A data.frame containing time series of:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{use_total}{Effective ITN usage over time.}
#'   \item{reach}{Campaign reach parameter.}
#'   \item{max_usage}{Maximum usage cap.}
#'   \item{halflife}{ITN half-life used in simulation.}
#' }
#'
#' @export

generate_AlternateContRoutine_SpecificAges=function(halflife_weibull=2.1,
                                                        shape_weibull=2,max_usage=0.8,
                                                        frequency_mass=3,
                                                        reach_cont=0.8*0.5, ageGroupProp_cont=0.27,
                                                        reach_routine=0.8, ageGroupProp_routine=0.035,
                                                        ageGroup_cont="school",otherageGroup_cont="others",
                                                        ageGroup_routine="U1",otherageGroup_routine="others"){

  cont=generate_AlternateCont_SpecificAges(halflife_weibull=halflife_weibull,
                                           shape_weibull=shape_weibull,
                                           reach=reach_cont, max_usage=max_usage,
                                           frequency=frequency_mass,
                                           ageGroupProp_cont =ageGroupProp_cont,
                                           ageGroup_cont=ageGroup_cont,otherageGroup_cont=otherageGroup_cont, rename = FALSE)

  routine=generate_continuousDistr_SpecificAges(halflife_weibull=halflife_weibull,
                                                shape_weibull=shape_weibull,
                                                reach=reach_routine, max_usage=max_usage,
                                                ageGroupProp=ageGroupProp_routine,
                                                ageGroup=ageGroup_routine,otherageGroup=otherageGroup_routine, full_output = T)|>
    dplyr::select(-use_targetpop, -use_otherpop, -spillover_cov_othergroup)


  names(routine)=gsub("group", "group_routine",names(routine))


  combined=merge(cont |> dplyr::mutate(reach_cont=reach_cont),
                 routine |> dplyr::mutate(reach_routine=reach_routine))|>
    dplyr::mutate(totcov_targetgroup=ifelse(is.na(totcov_targetgroup),0,totcov_targetgroup),
                  totcov_othergroup=ifelse(is.na(totcov_othergroup),0,totcov_othergroup),
                  cov_targetgroup=totcov_targetgroup+totcov_othergroup_routine,
                  cov_targetgroup_routine=totcov_targetgroup_routine+totcov_othergroup,
                  group1maxed=(cov_targetgroup>=1),
                  group2maxed=(cov_targetgroup_routine>=1),
                  denominator_spillover_othergroup1=ifelse(group2maxed,(1-ageGroupProp_cont-ageGroupProp_routine),1-ageGroupProp_cont),
                  denominator_spillover_othergroup2=ifelse(group1maxed,(1-ageGroupProp_cont-ageGroupProp_routine),1-ageGroupProp_routine),
                  spillover2_cov_othergroup=pmax(0, cov_targetgroup-1)*ageGroupProp_cont/denominator_spillover_othergroup1,
                  spillover2_cov_othergroup_routine=pmax(0, cov_targetgroup_routine-1)*ageGroupProp_routine/denominator_spillover_othergroup2,
                  cov_othergroup=totcov_othergroup+totcov_othergroup_routine+spillover2_cov_othergroup+spillover2_cov_othergroup_routine,
                  cov_targetgroup=cov_targetgroup+ifelse(group1maxed,0,spillover2_cov_othergroup_routine),
                  cov_targetgroup_routine=cov_targetgroup_routine+ifelse(group2maxed,0,spillover2_cov_othergroup),
                  use_targetpop=pmin(cov_targetgroup*max_usage, max_usage),
                  use_targetpop_routine=pmin(cov_targetgroup_routine*max_usage, max_usage),
                  use_otherpop=pmin(cov_othergroup*max_usage, max_usage)) |>
    dplyr::select(time, use_targetpop,use_targetpop_routine, use_otherpop)|>
    dplyr::mutate(reach_cont=reach_cont, reach_routine=reach_routine, max_usage=max_usage , halflife=halflife_weibull)

  names(combined)[c(2,3,4)]=c(paste0("use_", ageGroup_cont),paste0("use_", ageGroup_routine), "use_others")
  return(combined)
}
