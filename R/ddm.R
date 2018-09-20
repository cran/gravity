#' @title Double Demeaning (DDM)
#'
#' @description \code{ddm} estimates gravity models via double demeaning the
#' left hand side and right hand side of the gravity equation.
#'
#' @details \code{ddm} is an estimation method for gravity models presented
#' in \insertCite{Head2014;textual}{gravity}.
#'
#' Country specific effects are subdued due double demeaning. Hence, unilateral income
#' proxies such as GDP cannot be considered as exogenous variables.
#'
#' Unilateral effect drop out due to double demeaning and therefore cannot be estimated.
#'
#' \code{ddm} is designed to be consistent with the Stata code provided at
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' As, to our knowledge at the moment, there is no explicit literature covering
#' the estimation of a gravity equation by \code{ddm} using panel data,
#' we do not recommend to apply this method in this case.
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
#' The function will remove zero flows and distances.
#'
#' @param ... additional arguments to be passed to \code{ddm}.
#'
#' @references
#' For more information on Double Demeaning as well as information on gravity
#' models, theoretical foundations and estimation methods in general see
#'
#' \insertRef{Head2014}{gravity}
#'
#' as well as
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
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#'
#' @examples
#' \dontrun{
#' data(gravity_no_zeros)
#'
#' ddm(dependent_variable = "flow", regressors = c("distw", "rta"),
#' codes = c("iso_o", "iso_d"),
#' robust = TRUE, data = gravity_no_zeros)
#'
#' ddm(dependent_variable = "flow", regressors = c("distw", "rta", "comcur", "contig"),
#' codes = c("iso_o", "iso_d"),
#' robust=TRUE, data=gravity_no_zeros)
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
#' ddm(dependent_variable = "flow", regressors = c("distw", "rta"),
#' codes = c("iso_o", "iso_d"),
#' robust = TRUE, data = grav_small)
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

ddm <- function(dependent_variable, regressors, codes, robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(codes) | all(codes %in% colnames(data)) | length(codes) == 2)

  # Split input vectors -----------------------------------------------------
  code_o <- codes[1]
  code_d <- codes[2]

  distance <- regressors[1]
  additional_regressors <- regressors[-1]

  # Discarding unusable observations ----------------------------------------
  d <- data %>%
    filter_at(vars(!!sym(distance)), any_vars(!!sym(distance) > 0)) %>%
    filter_at(vars(!!sym(distance)), any_vars(is.finite(!!sym(distance)))) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(!!sym(dependent_variable) > 0)) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(is.finite(!!sym(dependent_variable))))

  # Transforming data, logging distances ---------------------------------------
  d <- d %>%
    mutate(
      dist_log = log(!!sym(distance))
    )

  # Transforming data, logging flows -------------------------------------------
  d <- d %>%
    mutate(
      y_log = log(!!sym(dependent_variable))
    )

  # Substracting the means -----------------------------------------------------
  d <- d %>%
    mutate(
      y_log_ddm = !!sym("y_log"),
      dist_log_ddm = !!sym("dist_log")
    ) %>%
    group_by(!!sym(code_o), add = FALSE) %>%
    mutate(
      ym1 = mean(!!sym("y_log_ddm"), na.rm = TRUE),
      dm1 = mean(!!sym("dist_log_ddm"), na.rm = TRUE)
    ) %>%
    group_by(!!sym(code_d), add = FALSE) %>%
    mutate(
      ym2 = mean(!!sym("y_log_ddm"), na.rm = TRUE),
      dm2 = mean(!!sym("dist_log_ddm"), na.rm = TRUE)
    ) %>%
    group_by(!!sym(code_o), add = FALSE) %>%
    mutate(
      y_log_ddm = !!sym("y_log_ddm") - !!sym("ym1"),
      dist_log_ddm = !!sym("dist_log_ddm") - !!sym("dm1")
    ) %>%
    group_by(!!sym(code_d), add = FALSE) %>%
    mutate(
      y_log_ddm = !!sym("y_log_ddm") - !!sym("ym2"),
      dist_log_ddm = !!sym("dist_log_ddm") - !!sym("dm2")
    ) %>%
    ungroup() %>%
    mutate(
      y_log_ddm = !!sym("y_log_ddm") + mean(!!sym("y_log"), na.rm = TRUE),
      dist_log_ddm = !!sym("dist_log_ddm") + mean(!!sym("dist_log"), na.rm = TRUE)
    )

  # Substracting the means for the other independent variables -----------------
  d2 <- d %>%
    select(!!sym(code_o), !!sym(code_d), additional_regressors) %>%
    gather(!!sym("key"), !!sym("value"), -!!sym(code_o), -!!sym(code_d)) %>%
    mutate(key = paste0(!!sym("key"), "_ddm")) %>%
    group_by(!!sym(code_o), !!sym("key"), add = FALSE) %>%
    mutate(ddm = !!sym("value") - mean(!!sym("value"), na.rm = TRUE)) %>%
    group_by(!!sym(code_d), !!sym("key"), add = FALSE) %>%
    mutate(ddm = !!sym("ddm") - mean(!!sym("value"), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(value = !!sym("ddm") + mean(!!sym("value"), na.rm = TRUE)) %>%
    select(!!!syms(c(code_o, code_d, "key", "value"))) %>%
    spread(!!sym("key"), !!sym("value"))

  # Model ----------------------------------------------------------------------
  dmodel <- left_join(d, d2, by = c(code_o, code_d)) %>%
    select(!!sym("y_log_ddm"), ends_with("_ddm"))

  model_ddm <- stats::lm(y_log_ddm ~ . + 0, data = dmodel)

  # Return ---------------------------------------------------------------------
  if (robust == TRUE) {
    return_object <- robust_summary(model_ddm, robust = TRUE)
    return_object$call <- as.formula(model_ddm)
    return(return_object)
  }

  if (robust == FALSE) {
    return_object <- robust_summary(model_ddm, robust = FALSE)
    return_object$call <- as.formula(model_ddm)
    return(return_object)
  }
}
