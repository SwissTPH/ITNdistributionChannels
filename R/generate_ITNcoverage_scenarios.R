#' Calculate ITN Coverage Decay Using a Weibull Function
#'
#' Computes the decay in insecticide-treated net (ITN) coverage over a 20-year
#' time horizon assuming a Weibull survival function parameterization.
#'
#' @param L Numeric. Half-life of the ITN in years.
#' @param kappa Numeric. Shape parameter of the Weibull decay function.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{ITNdecay}{A data.frame with daily time (1 to 20*365) and corresponding decay values.}
#'   \item{average_cov}{Numeric. Average decay (mean coverage) over the 20-year period.}
#' }
#'
#' @details
#' The decay function is defined as:
#' \deqn{\exp\left(-\left(\frac{t}{L \times 365}\right)^{\kappa} \log(2)\right)}
#' where \eqn{L} is the half-life and \eqn{\kappa} is the Weibull shape parameter.
#'
calculate_cov=function(L, kappa){
  ITNdecay= data.frame(
    time=1:(20*365))

  ITNdecay$decay=exp( -(ITNdecay$time/(L*365))^kappa * log(2) )

  average_cov= dplyr::summarise(ITNdecay, average_cov=mean(decay))


  return(list("ITNdecay"=ITNdecay, "average_cov"= as.numeric(average_cov))
              )
}


#' Generate Repeated ITN Distribution Coverage
#'
#' Simulates repeated ITN mass distributions over a 20-year horizon and
#' aggregates overlapping decay curves.
#'
#' @param halflife_weibull Numeric. Half-life of ITNs in years.
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param coverage Numeric. Initial coverage level per distribution round.
#' @param frequency Numeric. Distribution frequency in years.
#'
#' @return A data.frame containing:
#' \describe{
#'   \item{time}{Time in days.}
#'   \item{cov}{Total effective coverage at each time point.}
#' }
#'
#' @details
#' The function stacks multiple decay curves according to the specified
#' distribution frequency and sums overlapping coverage.
#'
#' @export
generate_ITN_distribution=function(halflife_weibull, shape_weibull, coverage, frequency){

  ITNdecay=calculate_cov(L=halflife_weibull, kappa=shape_weibull)$ITNdecay

  ITNdecay$cov=ITNdecay$decay*coverage

  ITNdecay_stacked=ITNdecay
  for(i in 1:ceiling(20/frequency)){
    ITNdecay_stacked=rbind(ITNdecay_stacked,
                           dplyr::mutate(ITNdecay, time=time+frequency*i*365))

  }
  ITNdecay_stacked_collapse=ITNdecay_stacked |>
    dplyr::group_by(time) |>
    dplyr::summarise(cov=sum(decay)*coverage) |>
    dplyr::filter(time <=20*365)

  return(ITNdecay_stacked_collapse)
}


#' Mass Campaign ITN Distribution (Uniform Across All Ages)
#'
#' Simulates periodic mass ITN campaigns where coverage is uniformly applied
#' across the entire population.
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach Numeric. Proportion of the population reached.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency Numeric. Distribution frequency (years).
#' @param ageGroupProp Numeric. Proportion of population in the targeted age group.
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
generate_massCampaign_UniformAllAges=function(halflife_weibull=2.1,
                                       shape_weibull=2,
                                       reach=0.8, max_usage=0.8,
                                       frequency=3, ageGroupProp=1){

  multiplier=1.8
  effective_multiplier=ifelse(ageGroupProp==1, 1,multiplier )
  ITNscenario=generate_ITN_distribution(halflife_weibull=halflife_weibull,
                                        shape_weibull=shape_weibull,
                                        coverage=reach*ageGroupProp*effective_multiplier,
                                        frequency=frequency)

  ITNscenario=dplyr::mutate(ITNscenario,
                              use_total=pmin(cov*max_usage, max_usage), reach=reach, max_usage=max_usage , halflife=halflife_weibull)


  return(ITNscenario)
}

#' Mass Campaign ITN Distribution (Specific Age Groups)
#'
#' Simulates periodic mass ITN campaigns targeting a specific age group,
#' including spillover effects to other age groups.
#'
#' @param halflife_weibull Numeric. ITN half-life (years).
#' @param shape_weibull Numeric. Weibull shape parameter.
#' @param reach Numeric. Proportion of the target population reached.
#' @param max_usage Numeric. Maximum achievable ITN usage.
#' @param frequency Numeric. Distribution frequency (years).
#' @param ageGroup Character. Name of the targeted age group.
#' @param otherageGroup Character. Name of the non-targeted age group.
#' @param ageGroupProp Numeric. Proportion of total population in target group.
#' @param full_output Boolean. Whether all intermediate coverage variables are returned (default =FALSE)
#'
#' @return A data.frame with time series of ITN usage for target and
#' non-target age groups.
#'
#' @export
generate_massCampaign_SpecificAges=function(halflife_weibull=2.1,
                                              shape_weibull=2,
                                              reach=0.8, max_usage=0.8,
                                              frequency=3, ageGroup="U5",otherageGroup="others",
                                            ageGroupProp=0.16, full_output=FALSE
                                            ){
  multiplier=1.8
  cov_targetpop=reach
  cov_otherpop=ageGroupProp*reach*max(0,multiplier-1) /(1-ageGroupProp)

  ITNscenario_targetpop=generate_ITN_distribution(halflife_weibull=halflife_weibull,
                                        shape_weibull=shape_weibull,
                                        coverage=cov_targetpop,
                                        frequency=frequency)

  ITNscenario_targetpop=dplyr::mutate(ITNscenario_targetpop,
                                      use_targetpop=pmin(cov*max_usage, max_usage))

  ITNscenario_otherpop=generate_ITN_distribution(halflife_weibull=halflife_weibull,
                                                  shape_weibull=shape_weibull,
                                                  coverage=cov_otherpop,
                                                  frequency=frequency)

  ITNscenario=merge(ITNscenario_targetpop|> dplyr::rename(totcov_targetgroup=cov),
                    ITNscenario_otherpop|> dplyr::rename(totcov_othergroup=cov)) |>
    dplyr::mutate(spillover_cov_othergroup=pmax(0, totcov_targetgroup-1)*ageGroupProp/(1-ageGroupProp),
           use_otherpop=pmin(max_usage*(totcov_othergroup+spillover_cov_othergroup), max_usage))

  if(full_output==FALSE){
    ITNscenario=ITNscenario |>
      dplyr::select(time, use_targetpop,use_otherpop)|>
      dplyr::mutate( reach=reach, max_usage=max_usage , halflife=halflife_weibull)
    names(ITNscenario)[c(2,3)]=c(paste0("use_", ageGroup),paste0("use_", otherageGroup))

  }

  return(ITNscenario)
}


#' Continuous ITN Distribution (Uniform Across All Ages)
#'
#' Simulates continuous annual ITN distribution applied uniformly across
#' the population.
#'
#' @inheritParams generate_massCampaign_UniformAllAges
#'
#' @return A data.frame containing time series of total ITN usage.
#'
#' @export
generate_continuousDistr_uniformAllAges=function(halflife_weibull=2.1,
                                               shape_weibull=2,
                                               reach=0.8, max_usage=0.8,
                                               ageGroupProp=0.035
){
  multiplier=1.8

  ITNscenario=generate_ITN_distribution(halflife_weibull=halflife_weibull,
                                        shape_weibull=shape_weibull,
                                        coverage=reach*ageGroupProp*multiplier,
                                        frequency=1)|>
    dplyr::mutate(use_total=pmin(cov*max_usage, max_usage), reach=reach, max_usage=max_usage , halflife=halflife_weibull)

  return(ITNscenario)
}


#' Continuous ITN Distribution (Specific Age Groups)
#'
#' Simulates continuous annual ITN distribution targeting a specific age group,
#' including spillover to non-targeted groups.
#'
#' @inheritParams generate_massCampaign_SpecificAges
#'
#' @return A data.frame containing time series of ITN usage by age group.
#'
#' @export
generate_continuousDistr_SpecificAges=function(halflife_weibull=2.1,
                                            shape_weibull=2,
                                            reach=0.8, max_usage=0.8,
                                            ageGroup="U1",otherageGroup="others",
                                            ageGroupProp=0.035, full_output=FALSE
){
  multiplier=1.8
  cov_targetpop=reach
  cov_otherpop=ageGroupProp*reach*max(0,multiplier-1) /(1-ageGroupProp)

  ITNscenario_targetpop=generate_ITN_distribution(halflife_weibull=halflife_weibull,
                                                  shape_weibull=shape_weibull,
                                                  coverage=cov_targetpop,
                                                  frequency=1)|>
    dplyr::mutate(use_targetpop=pmin(cov*max_usage, max_usage))

  ITNscenario_otherpop=generate_ITN_distribution(halflife_weibull=halflife_weibull,
                                                 shape_weibull=shape_weibull,
                                                 coverage=cov_otherpop,
                                                 frequency=1)

  ITNscenario=merge(ITNscenario_targetpop|> dplyr::rename(totcov_targetgroup=cov),
                    ITNscenario_otherpop|> dplyr::rename(totcov_othergroup=cov))|>
    dplyr::mutate(spillover_cov_othergroup=pmax(0, totcov_targetgroup-1)*ageGroupProp/(1-ageGroupProp),
           use_otherpop=pmin(max_usage*(totcov_othergroup+spillover_cov_othergroup), max_usage))

  if(full_output==FALSE){
    ITNscenario=ITNscenario |>
      dplyr::select(time, use_targetpop,use_otherpop)|>
      dplyr::mutate( reach=reach, max_usage=max_usage , halflife=halflife_weibull)
    names(ITNscenario)[c(2,3)]=c(paste0("use_", ageGroup),paste0("use_", otherageGroup))

  }

  return(ITNscenario)
}



#' Get average annual coverage for continuous intervention
#'
#' Calculates the mean coverage for a given year from a continuous-use
#' simulation scenario. The function filters the data to the specified year,
#' averages the selected coverage column across time, and returns summary
#' values grouped by \code{halflife}, \code{reach}, and \code{max_usage}.
#'
#' @param scenario A data frame containing simulation output with at least
#'   the columns \code{time}, \code{halflife}, \code{reach}, \code{max_usage},
#'   and the specified coverage column.
#' @param year Integer. The calendar year for which average coverage should
#'   be calculated. Years are derived as \code{floor(time / 365)}.
#' @param col Character. Name of the coverage column to average.
#'   Defaults to \code{"use_total"}.
#'
#' @return A data frame with one row per combination of
#'   \code{halflife}, \code{reach}, and \code{max_usage}, containing the
#'   average coverage for the specified year in the selected column.
#'
#' @examples
#' # get_average_coverage_continuous(sim_data, year = 5)
#'
#' @export
get_average_coverage_continuous=function(scenario, year, col="use_total"){

  names(scenario)[names(scenario) == col]=paste0("use")

  scenario=scenario |>
    dplyr::filter(floor(time/365)==year)|>
    dplyr::group_by(halflife, reach, max_usage)|>
    dplyr::summarise(use=mean(use, na.rm=T))

  names(scenario)[names(scenario) == "use"]=col

  return(scenario)
}


#' Get starting coverage for mass intervention
#'
#' Extracts coverage values at the start of a specified year from a
#' mass-distribution simulation scenario. The function filters the data
#' to the first day of the given year (\code{year * 365 + 1}).
#'
#' @param scenario A data frame containing simulation output with at least
#'   a \code{time} column and coverage columns.
#' @param year Integer. The calendar year for which starting coverage
#'   should be extracted.
#' @param col Character. Name of the coverage column (not directly used in
#'   the calculation but included for interface consistency).
#'   Defaults to \code{"use_total"}.
#'
#' @return A filtered data frame containing only rows corresponding to
#'   the first day of the specified year.
#'
#' @examples
#' # get_start_coverage_mass(sim_data, year = 3)
#'
#' @export
get_start_coverage_mass=function(scenario, year){

  scenario=scenario |>
    dplyr::filter(time==year*365+1)

  return(scenario)
}


#' Get average coverage over a year range for mass intervention
#'
#' Calculates mean coverage between two calendar years for a
#' mass-distribution simulation scenario. The function filters the data
#' between \code{year_start * 365} (inclusive) and
#' \code{year_end * 365} (exclusive), then averages the selected
#' coverage column grouped by \code{halflife}, \code{reach}, and
#' \code{max_usage}.
#'
#' @param scenario A data frame containing simulation output with at least
#'   the columns \code{time}, \code{halflife}, \code{reach}, \code{max_usage},
#'   and the specified coverage column.
#' @param year_start Integer. First year (inclusive) of the averaging window.
#' @param year_end Integer. Final year (exclusive) of the averaging window.
#' @param col Character. Name of the coverage column to average.
#'   Defaults to \code{"use_total"}.
#'
#' @return A data frame with one row per combination of
#'   \code{halflife}, \code{reach}, and \code{max_usage}, containing the
#'   average coverage over the specified time range.
#'
#' @examples
#' # get_average_coverage_mass(sim_data, year_start = 2, year_end = 5)
#'
#' @export
get_average_coverage_mass=function(scenario, year_start, year_end, col="use_total"){

  names(scenario)[names(scenario) == col]=paste0("use")

  scenario=scenario |>
    dplyr::filter(time>=year_start*365, time<year_end*365)|>
    dplyr::group_by(halflife, reach, max_usage)|>
    dplyr::summarise(use=mean(use, na.rm=T))

  names(scenario)[names(scenario) == "use"]=col

  return(scenario)
}

