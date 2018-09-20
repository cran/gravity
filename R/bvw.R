#' @title Bonus vetus OLS (BVW)
#'
#' @description \code{bvw} estimates gravity models via Bonus
#' vetus OLS with GDP-weights.
#'
#' @details Bonus vetus OLS is an estimation method for gravity models
#' developed by \insertCite{Baier2009,Baier2010;textual}{gravity} using GDP-weights to center a
#' Taylor-series (see the references for more information).
#'
#' The \code{bvw} function considers Multilateral Resistance terms and allows to
#' conduct comparative statics. Country specific effects are subdued due
#' to demeaning. Hence, unilateral variables apart from \code{inc_o}
#' and \code{inc_d} cannot be included in the estimation.
#'
#' \code{bvw} is designed to be consistent with the Stata code provided at
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' As, to our knowledge at the moment, there is no explicit literature covering
#' the estimation of a gravity equation by \code{bvw} using panel data,
#' we do not recommend to apply this method in this case.
#'
#' @param dependent_variable name (type: character) of the dependent variable in the dataset
#' \code{data} (e.g. trade flows).
#'
#' This dependent variable is divided by the
#' product of unilateral incomes (e.g.
#' GDPs or GNPs of the countries of interest, named \code{inc_o} and \code{inc_d} in the example datasets)
#' and logged afterwards.
#'
#' The transformed variable is then used as the dependent variable in the
#' estimation.
#'
#' @param regressors name (type: character) of the regressors to include in the model.
#'
#' Include the distance variable in the dataset \code{data} containing a measure of
#' distance between all pairs of bilateral partners and bilateral variables that should
#' be taken as the independent variables in the estimation.
#'
#' The distance is logged automatically when the function is executed.
#'
#' Unilateral metric variables such as GDPs should be inserted via the argument \code{incomes}.
#'
#' As country specific effects are subdued due to demeaning, no further unilateral variables
#' apart from unilateral incomes can be added.
#'
#' Write this argument as \code{c(distance, contiguity, common curreny, ...)}.
#'
#' @param incomes variable name (type: character) of the income of the country of
#' origin (e.g. \code{inc_o}) and destination (e.g. \code{inc_d}) in the dataset \code{data}.
#'
#' The dependent variable \code{dependent_variable} is divided by the product of the incomes.
#'
#' Write this argument as \code{c(income origin, income destination)}.
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
#' @param ... additional arguments to be passed to \code{bvw}.
#'
#' @references
#' For estimating gravity equations via Bonus Vetus OLS see
#'
#' \insertRef{Baier2009}{gravity}
#'
#' \insertRef{Baier2010}{gravity}
#'
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
#' \insertRef{Head2014}{gravity}
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
#' bvw(dependent_variable = "flow", regressors = c("distw", "rta"),
#' incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
#' robust = TRUE, data = gravity_no_zeros)
#'
#' bvw(dependent_variable = "flow", regressors = c("distw", "rta", "comcur", "contig"),
#' incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
#' robust = TRUE, data = gravity_no_zeros)
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
#' bvw(dependent_variable = "flow", regressors = c("distw", "rta"),
#' incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
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

bvw <- function(dependent_variable, regressors, incomes, codes, robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(incomes) | all(incomes %in% colnames(data)) | length(incomes) == 2)
  stopifnot(is.character(codes) | all(codes %in% colnames(data)) | length(codes) == 2)

  # Split input vectors -----------------------------------------------------
  inc_o <- incomes[1]
  inc_d <- incomes[2]

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
      y_log_bvw = log(
        !!sym(dependent_variable) / (!!sym(inc_o) * !!sym(inc_d))
      )
    )

  # GDP weights ----------------------------------------------------------------
  d <- d %>%
    group_by(!!sym(code_o)) %>%
    mutate(inc_world = sum(!!sym(inc_d), na.rm = TRUE)) %>%
    ungroup()

  # same for inc_o or inc_d as we have a squared dataset
  d <- d %>%
    mutate(
      theta_i = !!sym(inc_o) / !!sym("inc_world"),
      theta_j = !!sym(inc_d) / !!sym("inc_world")
    )

  # Multilateral resistance (MR) for distance ----------------------------------
  d <- d %>%
    group_by(!!sym(code_o), add = FALSE) %>%
    mutate(
      mr_dist_1 = sum(!!sym("theta_j") * !!sym("dist_log"), na.rm = TRUE)
    ) %>%
    group_by(!!sym(code_d), add = FALSE) %>%
    mutate(
      mr_dist_2 = sum(!!sym("theta_i") * !!sym("dist_log"), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(mr_dist_3 = sum(!!sym("theta_i") * !!sym("theta_j") * !!sym("dist_log"), na.rm = TRUE)) %>%
    mutate(dist_log_mr = !!sym("dist_log") - !!sym("mr_dist_1") - !!sym("mr_dist_2") + !!sym("mr_dist_3"))

  # Multilateral resistance (MR) for the other independent variables -----------
  d2 <- d %>%
    select(!!sym(code_o), !!sym(code_d), !!sym("theta_j"), !!sym("theta_i"), additional_regressors) %>%
    gather(!!sym("key"), !!sym("value"), -!!sym(code_o), -!!sym(code_d), -!!sym("theta_j"), -!!sym("theta_i")) %>%
    mutate(key = paste0(!!sym("key"), "_mr")) %>%
    group_by(!!sym(code_o), !!sym("key")) %>%
    mutate(mr1 = sum(!!sym("theta_j") * !!sym("value"), na.rm = TRUE)) %>%
    group_by(!!sym(code_d), !!sym("key")) %>%
    mutate(mr2 = sum(!!sym("theta_i") * !!sym("value"), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(mr3 = sum(!!sym("theta_i") * !!sym("theta_j") * !!sym("value"), na.rm = TRUE)) %>%
    mutate(value = !!sym("value") - !!sym("mr1") - !!sym("mr2") + !!sym("mr3")) %>%
    select(!!!syms(c(code_o, code_d, "key", "value"))) %>%
    spread(!!sym("key"), !!sym("value"))

  # Model ----------------------------------------------------------------------
  dmodel <- left_join(d, d2, by = c(code_o, code_d)) %>%
    select(!!sym("y_log_bvw"), ends_with("_mr"))

  model_bvw <- stats::lm(y_log_bvw ~ ., data = dmodel)

  # Return ---------------------------------------------------------------------
  if (robust == TRUE) {
    return_object <- robust_summary(model_bvw, robust = TRUE)
    return_object$call <- as.formula(model_bvw)
    return(return_object)
  }

  if (robust == FALSE) {
    return_object <- robust_summary(model_bvw, robust = FALSE)
    return_object$call <- as.formula(model_bvw)
    return(return_object)
  }
}
