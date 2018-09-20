#' @title \insertCite{Eaton2001;textual}{gravity} Tobit model (EK Tobit)
#'
#' @description \code{ek_tobit} estimates gravity models in their additive form
#' by conducting a censored regression.
#'
#' It follows the \insertCite{Eaton2001;textual}{gravity} Tobit model where each country
#' is assigned specific censoring bounds.
#'
#' @details \code{ek_tobit} represents the \insertCite{Eaton2001;textual}{gravity} Tobit model.
#'
#' When taking the log of the gravity equation flows equal to zero
#' constitute a problem as their log is not defined. Therefore, in \code{ek_tobit} all values of
#' the dependent variable are redefined as intervals.
#'
#' The positive observations have both interval bounds equal to their original value.
#'
#' For zero flows the interval is left open. The right border of the interval is set to the
#' log of the minimum positive trade flow of the respective importing country.
#'
#' The defined data object of class \code{\link[survival]{Surv}} is then inserted
#' in \code{\link[survival]{survreg}} for the parameter estimation.
#'
#' \code{ek_tobit} is designed to be consistent with the Stata code provided at
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' For other Tobit functions, see \code{\link[gravity]{tobit}}
#' for a simple Tobit model where number \code{1} is added to all observations
#' and \code{\link[gravity]{et_tobit}} for the Eaton and Tamura (1994)
#' threshold Tobit model where instead of simply adding number \code{1}
#' to the data the threshold is estimated.
#'
#' The function is designed for cross-sectional data, but can be extended to panel data using the
#' \code{\link[survival]{survreg}} function.
#'
#' @param dependent_variable name (type: character) of the dependent variable in the dataset
#' \code{data} (e.g. trade flows).
#'
#' This variable is logged and then used as the dependent variable in the estimation.
#'
#' As the log of zero is not defined, all flows equal to zero are replaced by a left open interval
#' with the logged minimum trade flow of the respective importing country as right border.
#'
#' @param regressors name (type: character) of the regressors to include in the model.
#'
#' Include the distance variable in the dataset \code{data} containing a measure of
#' distance between all pairs of bilateral partners and bilateral variables that should
#' be taken as the independent variables in the estimation.
#'
#' Unilateral metric variables such as GDPs should be inserted via the argument \code{incomes}.
#'
#' Interaction terms can be added.
#'
#' Write this argument as \code{c(distance, contiguity, common curreny, ...)}.
#'
#' @param code_destination variable name (type: character) of the label of the country
#' of destination (e.g. ISO-3 code from the \code{iso_d} variable in the example datasets). The variables
#' are grouped by using \code{iso_d} to obtain estimates.
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
#' The function allows zero flows but will remove zero distances.
#'
#' @param ... additional arguments to be passed to \code{ek_tobit}.
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
#' ek_tobit(dependent_variable = "flow", regressors = c("distw", "rta","lgdp_o","lgdp_d"),
#' code_destination = "iso_d",
#' robust = TRUE, data = gravity_zeros)
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
#' ek_tobit(dependent_variable = "flow", regressors = c("distw", "rta", "lgdp_o", "lgdp_d"),
#' code_destination = "iso_d",
#' robust = TRUE, data = grav_small_zeros)
#' }
#'
#' @return
#' The function returns the summary of the estimated gravity model as a
#' \code{\link[survival]{survreg}}-object.
#'
#' @seealso \code{\link[survival]{Surv}}, \code{\link[survival]{survreg}}
#'
#' @export

ek_tobit <- function(dependent_variable, regressors, code_destination, robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(code_destination) | code_destination %in% colnames(data) | length(code_destination) == 1)

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
  d <- d %>%
    mutate(
      y_log_ek = ifelse(!!sym(dependent_variable) > 0, log(!!sym(dependent_variable)), NA)
    )

  # Minimum flows -----------------------------------------------------------
  d <- d %>%
    group_by(!!sym(code_destination)) %>%
    mutate(
      exportmin = min(!!sym("y_log_ek"), na.rm = TRUE)
    )

  # Transforming censored variables -----------------------------------------
  d <- d %>%
    mutate(
      flows_ek1 = ifelse(!!sym(dependent_variable) > 0, !!sym("y_log_ek"), -Inf)
    ) %>%
    mutate(
      flows_ek2 = ifelse(!!sym(dependent_variable) > 0, !!sym("flows_ek1"), !!sym("exportmin"))
    ) %>%
    ungroup()

  # Response variable -------------------------------------------------------
  f1 <- d %>% select(!!sym("flows_ek1")) %>% as_vector()
  f2 <- d %>% select(!!sym("flows_ek2")) %>% as_vector()

  y_cens_log_ek <- survival::Surv(f1, f2, type = "interval2") %>% as_vector()

  # Model -------------------------------------------------------------------
  vars <- paste(c("dist_log", additional_regressors), collapse = " + ")
  form <- stats::as.formula(paste("y_cens_log_ek", "~", vars))
  model_ek_tobit <- survival::survreg(form, data = d, dist = "gaussian", robust = robust)

  # Return ------------------------------------------------------------------
  return_object <- summary(model_ek_tobit)
  return_object$call <- form
  return(return_object)
}
