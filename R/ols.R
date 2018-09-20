#' @title Ordinary Least Squares (OLS)
#'
#' @description \code{ols} estimates gravity models in their traditional form
#' via Ordinary Least Squares (ols). It does not consider Multilateral
#' Resistance terms.
#'
#' @details \code{ols} estimates gravity models in their traditional, additive,
#' form via Ordinary Least Squares using the \code{lm} function.
#' Multilateral Resistance terms are not considered by this function.
#'
#' As the coefficients for the country's incomes were often found to be close to
#' unitary and unitary income elasticities are in line with some theoretical
#' foundations on international trade, it is sometimes assumed that the income
#' elasticities are equal to unity.
#'
#' In order to allow for the estimation with
#' and without the assumption of unitary income elasticities, the option
#' \code{uie} is built into \code{ols} with the default set to \code{FALSE}.
#'
#' \code{ols} estimation can be used for both, cross-sectional and
#' panel data. Nonetheless, the function is designed to be consistent with the
#' Stata code for cross-sectional data provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' The function \code{ols} was therefore tested for cross-sectional data. For the use with panel data
#' no tests were performed.
#'
#' Therefore, it is up to the user to ensure that the functions can be applied
#' to panel data.
#'
#' Depending on the panel dataset and the variables -
#' specifically the type of fixed effects -
#' included in the model, it may easily occur that the model is not computable.
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
#' If \code{uie = TRUE} the dependent variable is divided by the product of
#' unilateral incomes (e.g. GDP variables \code{inc_o} and \code{inc_d} in the example datasets) of the
#' countries of interest and logged afterwards.
#'
#' If \code{uie=FALSE} the dependent variable is logged directly. The transformed variable is then used as
#' the dependent variable and the logged income variables are used as independent variables in the
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
#' Unilateral metric variables such as GDPs can be added but those variables have to be logged first.
#'
#' Interaction terms can be added.
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
#' @param uie Unitary Income Elasticities (type: logic) determines whether the
#' parameters are to be estimated assuming unitary income elasticities. The default value is set
#' to \code{FALSE}.
#'
#' If \code{uie} is set \code{TRUE}, the flows in the dependent variable \code{y} are divided
#' by the product of the country pairs' incomes before the estimation.
#'
#' If \code{uie} is set to \code{FALSE}, the income variables are logged and taken as independent
#' variables in the estimation. The variable names for the incomes should be included (e.g. \code{inc_o}
#' and \code{inc_d} in the example datasets).
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
#' The function will remove zero flows and distances.
#'
#' @param ... additional arguments to be passed to \code{ols}.
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
#' ols(dependent_variable = "flow", regressors = c("distw", "rta", "contig", "comcur"),
#' incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
#' uie = TRUE, robust = TRUE, data = gravity_no_zeros)
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
#' ols(dependent_variable = "flow", regressors = c("distw", "rta"),
#' incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
#' uie = FALSE, robust = TRUE, data = grav_small)
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

ols <- function(dependent_variable, regressors, incomes, codes, uie = FALSE, robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(uie))
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

  # Income elasticities --------------------------------------------------------
  if (uie == TRUE) {
    # Transforming data, logging flows -----------------------------------------
    d <- d %>%
      mutate(
        y_log_ols = log(
          !!sym(dependent_variable) / (!!sym(inc_o) * !!sym(inc_d))
        )
      ) %>%
      select(
        !!sym("y_log_ols"), !!sym("dist_log"), !!sym("additional_regressors")
      )

    # Model --------------------------------------------------------------------
    model_ols <- stats::lm(y_log_ols ~ ., data = d)
  }

  if (uie == FALSE) {
    # Transforming data, logging flows -----------------------------------------
    d <- d %>%
      mutate(
        y_log_ols = log(!!sym(dependent_variable)),
        inc_o_log = log(!!sym(inc_o)),
        inc_d_log = log(!!sym(inc_d))
      ) %>%
      select(
        !!sym("y_log_ols"), !!sym("inc_o_log"), !!sym("inc_d_log"), !!sym("dist_log"), !!sym("additional_regressors")
      )

    # Model --------------------------------------------------------------------
    model_ols <- stats::lm(y_log_ols ~ ., data = d)
  }

  # Return ---------------------------------------------------------------------
  if (robust == TRUE) {
    return_object_1 <- robust_summary(model_ols, robust = TRUE)
    return_object_1$call <- as.formula(model_ols)
    return(return_object_1)
  }

  if (robust == FALSE) {
    return_object_1 <- robust_summary(model_ols, robust = FALSE)
    return_object_1$call <- as.formula(model_ols)
    return(return_object_1)
  }
}
