#' Objective Function for Weibull Decay Fitting
#'
#' Computes the sum of squared differences between observed decay values and
#' a Weibull decay model with scaling parameter.
#'
#' @param param Numeric vector of length 3:
#' \describe{
#'   \item{param[1]}{Half-life (L)}
#'   \item{param[2]}{Shape parameter (k)}
#'   \item{param[3]}{Initial ITN use (a)}
#' }
#' @param df Data frame with columns:
#' \describe{
#'   \item{time}{Time values}
#'   \item{value}{Observed decay values}
#' }
#'
#' @return Data frame with total squared error.
#' @export
optimise_decay_param=function(param, df){
  #if(param[3]<1){
  out=data.frame(
    time=df |> dplyr::pull(time),
    decay_obs=df |> dplyr::pull(value)
  ) |>
    dplyr::mutate(decay=param[3]*exp( -(time/(param[1]))^param[2] * log(2) ),
           diff=(decay-decay_obs)^2) |>
    dplyr::summarise(diff=sum(diff))
  return(out)
}

#' Estimate Weibull Durability Parameters
#'
#' Fits a Weibull decay model to observed data using numerical optimisation.
#'
#' @param df Data frame with columns:
#' \describe{
#'   \item{time}{Time values}
#'   \item{value}{Observed decay values}
#' }
#'
#' @return Named numeric vector:
#' \describe{
#'   \item{L}{Estimated half-life}
#'   \item{k}{Estimated shape parameter}
#'   \item{a}{Initial ITN use}
#' }
#' @details Uses \code{optim()} with method \code{"L-BFGS-B"} and bounded parameters.
#' @export
calculate_durability_param=function( df){
  opti=optim(c(0.448, 1.11, 0.8), optimise_decay_param, df=df , method = "L-BFGS-B", lower=c(0.01, 0, 0), upper = c(Inf, Inf, 1) )$par
  names(opti)=c("L", "k", "a")
  return(opti)
}

#' Plot Observed vs Fitted Decay Curve
#'
#' Generates a line plot comparing observed decay values with a fitted Weibull model.
#'
#' @param param Numeric vector of length 3:
#' \describe{
#'   \item{param[1]}{Half-life (L)}
#'   \item{param[2]}{Shape parameter (k)}
#'   \item{param[3]}{Initial ITN use (a)}
#' }
#' @param df Data frame with columns:
#' \describe{
#'   \item{time}{Time values}
#'   \item{value}{Observed decay values}
#' }
#'
#' @return A \code{ggplot2} object.
#' @export
plot_decay=function(param, df){
  out=data.frame(
    time=df |> dplyr::pull(time),
    decay_obs=df |> dplyr::pull(value)
  )|>
    dplyr::mutate(decay=param[3]*exp( -(time/(param[1]))^param[2] * log(2) ))

  p=ggplot2::ggplot(out)+
    ggplot2::geom_line(ggplot2::aes(x=time, y=decay, col="model"), lwd=1)+
    ggplot2::geom_line(ggplot2::aes(x=time, y=decay_obs, col="observations"), lwd=1)
  return(p)
}



