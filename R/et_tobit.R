#' @title \insertCite{Eaton1995;textual}{gravity} threshold Tobit model (ET Tobit)
#'
#' @description \code{et_tobit} estimates gravity models in their additive form
#' by conducting a left-censored regression.
#'
#' It follows the \insertCite{Eaton1995;textual}{gravity} Tobit model,
#' also called threshold Tobit model, where,
#' instead of adding number \code{1} to the dependent variable as done
#' in \code{\link[gravity]{tobit}}, the constant added to the
#' data is estimated and interpreted as a threshold.
#'
#' For estimating this threshold, we follow \insertCite{Carson2007;textual}{gravity}.
#'
#' @details \code{et_tobit} represents the \insertCite{Eaton1995;textual}{gravity} Tobit model
#' which is often used when several gravity models are compared.
#'
#' When taking the log of the gravity equation flows equal to zero constitute a problem as their
#' log is not defined. Therefore, a constant is added to the flows.
#'
#' This constant, opposed to \code{\link[gravity]{tobit}}, is estimated.
#' Compared to the usual ET-Tobit approaches, in this package, the estimation
#' of the threshold is done before the other parameters are estimated.
#'
#' We follow \insertCite{Carson2007;textual}{gravity}, who show that taking the minimum
#' positive flow value as an estimate of the threshold is super-consistent and that
#' using this threshold estimate ensures that the parameter MLE are asymptotically normal with
#' the asymptotic variance identical to the variance achieved when the threshold is known. Hence, first
#' the threshold is estimated as the minimum positive flow. This threshold is added to the flow variable,
#' it is logged afterwards and taken as the dependent variable.
#'
#' The Tobit estimation is then conducted using the
#' \code{\link[censReg]{censReg}} function and setting the lower bound
#' equal to the log of the minimum positive flow value which was added to all
#' observations.
#'
#' A Tobit regression represents a combination of a binary and a
#' linear regression. This procedure has to be taken into consideration when
#' interpreting the estimated coefficients.
#'
#' The marginal effects of an explanatory variable on the expected value of
#' the dependent variable equals the product of both the probability of the
#' latent variable exceeding the threshold and the marginal effect of the
#' explanatory variable of the expected value of the latent variable.
#'
#' For a more elaborate Tobit function, see \code{\link[gravity]{ek_tobit}}
#' for the Eaton and Kortum (2001) Tobit model where each zero trade volume
#' is assigned a country specific interval with the upper
#' bound equal to the minimum positive trade level of the respective
#' importing country.
#'
#' The function is designed for cross-sectional data, but can be extended to panel data using the
#' \code{\link[censReg]{censReg}} function.
#'
#' A robust estimations is not implemented to the present
#' as the \code{\link[censReg]{censReg}} function is not
#' compatible with the \code{\link[sandwich]{vcovHC}} function.
#'
#' @param dependent_variable name (type: character) of the dependent variable in the dataset
#' \code{data}, e.g. trade flows.
#'
#' Following Carson and Sun (2007), the smallest positive flow value is
#' used as an estimate of the threshold, this value is is added to the \code{dependent_variable},
#' the result is logged and taken as the dependent variable in the Tobit estimation with
#' lower bound equal to the log of the smallest possible flow value.
#'
#' @param regressors name (type: character) of the regressors to include in the model.
#'
#' Include the distance variable in the dataset \code{data} containing a measure of
#' distance between all pairs of bilateral partners and bilateral variables that should
#' be taken as the independent variables in the estimation.
#'
#' Unilateral metric variables such as GDPs can be added but those variables have to be logged first.
#'
#' Interaction terms can be added.
#'
#' Write this argument as \code{c(distance, contiguity, common curreny, ...)}.
#'
#' @param data name of the dataset to be used (type: character).
#'
#' To estimate gravity equations you need a square dataset including bilateral
#' flows defined by the argument \code{dependent_variable}, ISO codes or similar of type character
#' (e.g. \code{iso_o} for the country of origin and \code{iso_d} for the
#' destination country), a distance measure defined by the argument \code{distance}
#' and other potential influences (e.g. contiguity and common currency) given as a vector in
#' \code{regressors} are required.
#'
#' All dummy variables should be of type numeric (0/1).
#'
#' Make sure the ISO codes are of type "character".
#'
#' If an independent variable is defined as a ratio, it should be logged.
#'
#' The user should perform some data cleaning beforehand to remove observations that contain entries that
#' can distort estimates.
#'
#' The function allows zero flows but will remove zero distances.
#'
#' @param ... additional arguments to be passed to \code{et_tobit}.
#'
#' @references
#' For more information on gravity models, theoretical foundations and
#' estimation methods in general see
#'
#' \insertRef{Anderson1979}{gravity}
#'
#' \insertRef{Anderson2001}{gravity}
#'
#' \insertRef{Anderson2010}{gravity}
#'
#' \insertRef{Baier2009}{gravity}
#'
#' \insertRef{Baier2010}{gravity}
#'
#' \insertRef{Head2010}{gravity}
#'
#' \insertRef{Santos2006}{gravity}
#'
#' and the citations therein.
#'
#' Especially for Tobit models see
#'
#' \insertRef{Tobin1958}{gravity}
#'
#' \insertRef{Eaton1995}{gravity}
#'
#' \insertRef{Eaton2001}{gravity}
#'
#' \insertRef{Carson2007}{gravity}
#'
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#'
#' @examples
#' \dontrun{
#' # Example for data with zero trade flows
#' data(gravity_zeros)
#'
#' gravity_zeros <- gravity_zeros %>%
#'     mutate(
#'         lgdp_o = log(gdp_o),
#'         lgdp_d = log(gdp_d)
#'     )
#'
#' et_tobit(dependent_variable = "flow", regressors = c("distw", "rta","lgdp_o","lgdp_d"),
#' data = gravity_zeros)
#' }
#'
#' \dontshow{
#' # examples for CRAN checks:
#' # executable in < 5 sec together with the examples above
#' # not shown to users
#'
#' data(gravity_zeros)
#' gravity_zeros$lgdp_o <- log(gravity_zeros$gdp_o)
#' gravity_zeros$lgdp_d <- log(gravity_zeros$gdp_d)
#'
#' # choose exemplarily 10 biggest countries for check data
#' countries_chosen_zeros <- names(sort(table(gravity_zeros$iso_o), decreasing = TRUE)[1:10])
#' grav_small_zeros <- gravity_zeros[gravity_zeros$iso_o %in% countries_chosen_zeros,]
#' et_tobit(dependent_variable = "flow", regressors = c("distw", "rta","lgdp_o","lgdp_d"),
#' data = grav_small_zeros)
#' }
#'
#' @return
#' The function returns the summary of the estimated gravity model as a
#' \code{\link[censReg]{censReg}}-object.
#'
#' @seealso \code{\link[censReg]{censReg}}
#'
#' @export

et_tobit <- function(dependent_variable, regressors, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)

  # Split input vectors -----------------------------------------------------
  distance <- regressors[1]
  additional_regressors <- regressors[-1]

  # Discarding unusable observations ----------------------------------------
  d <- data %>%
    filter_at(vars(!!sym(distance)), any_vars(!!sym(distance) > 0)) %>%
    filter_at(vars(!!sym(distance)), any_vars(is.finite(!!sym(distance))))

  # Transforming data, logging distances ---------------------------------------
  d <- d %>%
    mutate(
      dist_log = log(!!sym(distance))
    )

  # Transforming data, logging flows -------------------------------------------
  flow_min_log <- filter_at(d, vars(!!sym(dependent_variable)), any_vars(!!sym(dependent_variable) > 0))

  d <- d %>%
    mutate(
      y_log_et = ifelse(!!sym(dependent_variable) > 0, log(!!sym(dependent_variable)), NA)
    )

  # Transforming data, logging flows, distances --------------------------------
  d <- d %>%
    mutate(
      y2 = ifelse(!!sym(dependent_variable) > 0, !!sym(dependent_variable), NA),
      y2_log = log(!!sym("y2"))
    )

  y2min <- min(select(d, !!sym("y2")), na.rm = TRUE)
  y2min_log <- log(y2min)

  d <- d %>%
    rowwise() %>%
    mutate(y_cens_log_et = log(sum(!!sym(dependent_variable), y2min, na.rm = TRUE))) %>%
    ungroup()

  # Model -------------------------------------------------------------------
  vars <- paste(c("dist_log", additional_regressors), collapse = " + ")
  form <- stats::as.formula(paste("y_cens_log_et", "~", vars))
  model_et_tobit <- censReg::censReg(
    formula = form,
    left = y2min_log,
    right = Inf,
    data = d,
    start = rep(0, 2 + length(regressors)),
    method = "BHHH"
  )

  # Return ------------------------------------------------------------------
  return_object <- summary(model_et_tobit)
  return_object$call <- form
  return(return_object)
}
