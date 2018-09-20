#' @title Fixed Effects
#'
#' @description \code{fixed_effects} estimates gravity models via
#' OLS and fixed effects for the countries of origin and destination.
#' These effects catch country specific effects.
#'
#' @details To account for MR terms, Feenstra (2002) and Feenstra (2004) propose to use
#' importer and exporter fixed effects. Due to the use of these effects, all
#' unilateral influences such as GDPs can no longer be estimated.
#'
#' A disadvantage of the use of \code{fixed_effects} is that, when applied to
#' panel data, the number of country-year or country-pair fixed effects can be
#' too high for estimation. In addition, no comparative statistics are
#' possible with \code{fixed_effects} as the Multilateral Resistance terms are not estimated
#' explicitly. Nevertheless, \insertCite{Head2014;textual}{gravity} highlight the importance of
#' the use of fixed effects.
#'
#' By including country specific fixed effects, all monadic effects
#' are captured, including Multilateral Resistance terms.
#' Therefore, no other unilateral variables such as GDP can be
#' included as independent variables in the estimation.
#'
#' \code{fixed_effects} estimation can be used for both, cross-sectional as well as
#' panel data.
#'
#' Nonetheless, the function is designed to be consistent with the
#' Stata code for cross-sectional data provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' The function \code{fixed_effects} was therefore tested for
#' cross-sectional data. Its up to the user to ensure that the functions can be applied
#' to panel data.
#'
#' Depending on the panel dataset and the variables -
#' specifically the type of fixed effects -
#' included in the model, it may easily occur that the model is not computable.
#'
#' Also, note that by including bilateral fixed effects such as country-pair
#' effects, the coefficients of time-invariant observables such as distance
#' can no longer be estimated.
#'
#' Depending on the specific model, the code of the
#' respective function may has to be changed in order to exclude the distance
#' variable from the estimation.
#'
#' At the very least, the user should take special
#' care with respect to the meaning of the estimated coefficients and variances
#' as well as the decision about which effects to include in the estimation.
#'
#' When using panel data, the parameter and variance estimation of the models
#' may have to be changed accordingly.
#'
#' For a comprehensive overview of gravity models for panel data
#' see \insertCite{Egger2003;textual}{gravity}, \insertCite{Gomez-Herrera2013;textual}{gravity} and
#' \insertCite{Head2010;textual}{gravity} as well as the references therein.
#'
#' @param dependent_variable name (type: character) of the dependent variable in the dataset
#' \code{data} (e.g. trade flows).
#'
#' This variable is logged and then used as the dependent variable in the estimation.
#'
#' @param regressors name (type: character) of the regressors to include in the model.
#'
#' Include the distance variable in the dataset \code{data} containing a measure of
#' distance between all pairs of bilateral partners and bilateral variables that should
#' be taken as the independent variables in the estimation.
#'
#' The distance is logged automatically when the function is executed.
#'
#' Fixed effects catch all unilateral effects. Therefore, no other unilateral variables such as
#' GDP can be included as independent variables in the estimation.
#'
#' Write this argument as \code{c(distance, contiguity, common curreny, ...)}.
#'
#' @param codes variable name (type: character) of the code of the country
#' of origin and destination (e.g. ISO-3 codes from the variables \code{iso_o} and \code{iso_d}) in the
#' example datasets).
#'
#' The variables are grouped by using \code{iso_o} and \code{iso_d} to obtain estimates.
#'
#' Write this argument as \code{c(code origin, code destination)}.
#'
#' @param robust robust (type: logical) determines whether a robust
#' variance-covariance matrix should be used. By default is set to \code{TRUE}.
#'
#' If \code{robust = TRUE} the estimation results are consistent with the
#' Stata code provided at \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
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
#' When using panel data, a variable for the time may be included in the
#' dataset. Note that the variable for the time dimension should be of
#' type factor.
#'
#' The time variable can be used as a single dependent variable or interaction
#' term with other variables such as country identifiers by inserting it into
#' \code{regressors} or as an optional parameter.
#'
#' The function will remove zero flows and distances.
#'
#' @param ... additional arguments to be passed to \code{fixed_effects}.
#'
#' @references
#' For more information on fixed effects as well as informaton on gravity models,
#' theoretical foundations and suitable estimation methods in general see
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
#' \insertRef{Head2014}{gravity}
#'
#' \insertRef{Santos2006}{gravity}
#'
#' and the citations therein.
#'
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#'
#' For estimating gravity equations using panel data see
#'
#' \insertRef{Egger2003}{gravity}
#'
#' \insertRef{Gomez-Herrera2013}{gravity}
#'
#' and the references therein.
#'
#' @examples
#' \dontrun{
#' data(gravity_no_zeros)
#'
#' fixed_effects(dependent_variable = "flow",
#' regressors = c("distw", "rta"), codes = c("iso_o", "iso_d"),
#' robust = TRUE, data = gravity_no_zeros)
#'
#' fixed_effects(dependent_variable = "flow",
#' regressors = c("distw", "rta", "comcur", "contig"),
#' codes = c("iso_o", "iso_d"), robust = TRUE, data = gravity_no_zeros)
#' }
#'
#' \dontshow{
#' # examples for CRAN checks:
#' # executable in < 5 sec together with the examples above
#' # not shown to users
#'
#' data(gravity_no_zeros)
#' # choose exemplarily 10 biggest countries for check data
#' countries_chosen <- names(sort(table(gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
#' grav_small <- gravity_no_zeros[gravity_no_zeros$iso_o %in% countries_chosen,]
#' fixed_effects(dependent_variable = "flow", regressors = c("distw", "rta"),
#' codes = c("iso_o", "iso_d"), robust = TRUE, data = grav_small)
#' }
#'
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[sandwich]{vcovHC}}
#'
#' @export

fixed_effects <- function(dependent_variable, regressors, codes = c("iso_o", "iso_d"), robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(codes) | all(codes %in% colnames(data)) | length(codes) == 2)

  # Split input vectors -----------------------------------------------------
  distance <- regressors[1]
  additional_regressors <- regressors[-1]

  # Discarding unusable observations ----------------------------------------
  d <- data %>%
    filter_at(vars(!!sym(distance)), any_vars(!!sym(distance) > 0)) %>%
    filter_at(vars(!!sym(distance)), any_vars(is.finite(!!sym(distance)))) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(!!sym(dependent_variable) > 0)) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(is.finite(!!sym(dependent_variable))))

  # Transforming data, logging flows and distances --------------------------
  d <- data %>%
    mutate(
      dist_log = log(!!sym(distance)),
      y_log_fe = log(!!sym(dependent_variable))
    )

  # Model ----------------------------------------------------------------------
  vars <- paste(c("dist_log", additional_regressors, codes), collapse = " + ")
  form <- stats::as.formula(paste("y_log_fe", "~", vars))
  model_fe <- stats::lm(form, data = d)

  # Return ---------------------------------------------------------------------
  return_object <- robust_summary(model_fe, robust = robust)
  return_object$call <- form

  return(return_object)
}
