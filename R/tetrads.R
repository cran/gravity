#' @title Tetrads
#'
#' @description \code{tetrads} estimates gravity models
#' by taking the ratio of the ratio of flows.
#'
#' @details \code{tetrads} is an estimation method for gravity models
#' developed by \insertCite{Head2010;textual}{gravity}.
#'
#' The function \code{tetrads} utilizes the multiplicative form of the
#' gravity equation. After choosing a reference exporter \code{A} and
#' importer \code{B} one can eliminate importer and exporter fixed effects
#' by taking the ratio of ratios.
#'
#' Only those exporters trading with the
#' reference importer and importers trading with the reference exporter will
#' remain for the estimation. Therefore, reference countries should
#' preferably be countries which trade with every other country in the dataset.
#'
#' After restricting the data in this way, \code{tetrads} estimates the gravity
#' equation in its additive form by OLS.
#'
#' By taking the ratio of ratios, all monadic effects diminish, hence no
#' unilateral variables such as GDP can be included as independent variables.
#'
#' \code{tetrads} estimation can be used for both, cross-sectional as well as
#' panel data. Nonetheless, the function is designed to be consistent with the
#' Stata code for cross-sectional data provided on the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' The function \code{tetrads} was therefore tested for cross-sectional data.
#'
#' If tetrads is used for panel data, the user may have to drop distance as an
#' independent variable as time-invariant effects drop.
#'
#' For applying \code{tetrads} to panel data see \insertCite{Head2010;textual}{gravity}.
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
#' Unilateral effects drop as the ratio of ratios is taken.
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
#' @param reference_countries reference exporting and importing country, default is set to
#' \code{c("JPN", "USA")}
#'
#' Write this argument as \code{c(importing country, exporting country)}.
#'
#' @param multiway (type: logic) In case \code{multiway = TRUE}, the
#' \code{\link[multiwayvcov]{cluster.vcov}} function is used for estimation following
#' \insertCite{Cameron2011;textual}{gravity} multi-way clustering of
#' variance-covariance matrices.
#'
#' The default value is set to \code{TRUE}.
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
#' @param ... additional arguments to be passed to functions used by
#' \code{tetrads}.
#'
#' @references
#' For information on \code{tetrads} see
#'
#' \insertRef{Cameron2011}{gravity}
#'
#' \insertRef{Head2010}{gravity}
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
#' tetrads(dependent_variable = "flow", regressors = c("distw", "rta"),
#' codes = c("iso_o", "iso_d"), reference_countries = c("JPN", "USA"),
#' multiway = TRUE, data = gravity_no_zeros)
#'
#' tetrads(dependent_variable = "flow", regressors = c("distw", "rta", "comcur", "contig"),
#' codes = c("iso_o", "iso_d"), reference_countries = c("JPN", "USA"),
#' multiway = FALSE, data = gravity_no_zeros)
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
#' tetrads(dependent_variable ="flow",
#' regressors = c("distw", "rta"),
#' codes = c("iso_o", "iso_d"),
#' reference_countries = c(countries_chosen[1], countries_chosen[2]),
#' multiway = FALSE,
#' data = grav_small)
#' }
#'
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[multiwayvcov]{cluster.vcov}}
#'
#' @export

tetrads <- function(dependent_variable, regressors, codes, reference_countries = c("JPN", "USA"), multiway = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(multiway))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(reference_countries) | all(reference_countries %in% colnames(data)) | length(reference_countries) == 2)

  # Split input vectors -----------------------------------------------------
  code_o <- codes[1]
  code_d <- codes[2]

  filter_o <- reference_countries[1]
  filter_d <- reference_countries[2]

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
      y_log_tetrads = log(!!sym(dependent_variable))
    )

  # Truncating dataset to reference importer and exporter partners -------------
  d2_filter_d <- d %>%
    filter_at(vars(!!sym(code_d)), any_vars(!!sym(code_d) == filter_d)) %>%
    select(!!sym(code_o)) %>%
    distinct() %>%
    as_vector()

  d2 <- d %>%
    filter_at(vars(!!sym(code_o)), any_vars(!!sym(code_o) %in% d2_filter_d))

  d2_filter_o <- d2 %>%
    filter_at(vars(!!sym(code_o)), any_vars(!!sym(code_o) == filter_o)) %>%
    select(!!sym(code_d)) %>%
    distinct() %>%
    as_vector()

  d2 <- d2 %>%
    filter_at(vars(!!sym(code_d)), any_vars(!!sym(code_d) %in% d2_filter_o))

  # Taking ratios, ratk --------------------------------------------------------
  d2 <- left_join(
    d2,
    d2 %>%
      filter_at(vars(!!sym(code_d)), any_vars(!!sym(code_d) == filter_d)) %>%
      select(!!sym("code_o"), y_log_tetrads_d = !!sym("y_log_tetrads"), dist_log_d = !!sym("dist_log")),
    by = code_o
  ) %>%
    mutate(
      lXinratk = !!sym("y_log_tetrads") - !!sym("y_log_tetrads_d"),
      ldistratk = !!sym("dist_log") - !!sym("dist_log_d")
    ) %>%
    select(-!!sym("y_log_tetrads_d"), -!!sym("dist_log_d"))

  # Taking ratios, ratk, for the other independent variables -------------------
  d2 <- d2 %>%
    select(c(!!sym("code_o"), !!sym("code_d"), regressors)) %>%
    gather(!!sym("key"), !!sym("value"), -!!sym("code_o"), -!!sym("code_d")) %>%
    left_join(
      d2 %>%
        filter_at(vars(!!sym(code_d)), any_vars(!!sym(code_d) == filter_d)) %>%
        select(c(!!sym("code_o"), regressors)) %>%
        gather(!!sym("key"), !!sym("value"), -!!sym("code_o")),
      by = c(code_o, "key")
    ) %>%
    mutate(value = !!sym("value.x") - !!sym("value.y")) %>%
    select(c(!!sym("code_o"), !!sym("code_d"), "key", "value")) %>%
    spread(!!sym("key"), !!sym("value")) %>%
    left_join(d2 %>% select(-one_of(regressors)), by = c(code_o, code_d))

  # Taking the ratio of ratios, rat --------------------------------------------
  d2 <- left_join(
    d2,
    d2 %>%
      filter_at(vars(!!sym(code_o)), any_vars(!!sym(code_o) == filter_o)) %>%
      select(!!sym("code_d"), lXinratk_o = !!sym("lXinratk"), ldistratk_o = !!sym("ldistratk")),
    by = code_d
  ) %>%
    mutate(
      y_log_rat = !!sym("lXinratk") - !!sym("lXinratk_o"),
      dist_log_rat = !!sym("ldistratk") - !!sym("ldistratk_o")
    ) %>%
    select(-!!sym("lXinratk_o"), -!!sym("ldistratk_o"))

  # Taking the ratio of ratios, rat, for the other independent variables -------
  d2 <- d2 %>%
    select(c(!!sym("code_o"), !!sym("code_d"), additional_regressors)) %>%
    gather(!!sym("key"), !!sym("value"), -!!sym("code_o"), -!!sym("code_d")) %>%
    left_join(
      d2 %>%
        filter_at(vars(!!sym(code_o)), any_vars(!!sym(code_o) == filter_o)) %>%
        select(c(!!sym("code_d"), additional_regressors)) %>%
        gather(!!sym("key"), !!sym("value"), -!!sym("code_d")),
      by = c(code_d, "key")
    ) %>%
    mutate(
      key = paste0(!!sym("key"), "_rat"),
      value = !!sym("value.x") - !!sym("value.y")
    ) %>%
    select(c(!!sym("code_o"), !!sym("code_d"), "key", "value")) %>%
    spread(!!sym("key"), !!sym("value")) %>%
    left_join(d2 %>% select(-one_of(regressors)), by = c(code_o, code_d)) %>%
    select(ends_with("_rat"), codes)

  # Model ----------------------------------------------------------------------
  additional_regressors <- paste0(additional_regressors, "_rat")
  form <- stats::as.formula(
    sprintf("y_log_rat ~ dist_log_rat + %s", paste(additional_regressors, collapse = " + "))
  )
  form2 <- as.formula(sprintf("~ %s + %s", code_o, code_d))
  model_tetrads <- stats::lm(form, data = d2)
  model_tetrads_vcov <- multiwayvcov::cluster.vcov(model = model_tetrads, cluster = form2)
  model_tetrads_robust <- lmtest::coeftest(x = model_tetrads, vcov = model_tetrads_vcov)

  # Return ---------------------------------------------------------------------
  if (multiway == TRUE) {
    summary_ted <- robust_summary(model_tetrads, robust = TRUE)
    summary_ted$coefficients <- model_tetrads_robust[1:length(rownames(model_tetrads_robust)), ]
    return_object <- summary_ted
    return_object$call <- form
    return(return_object)
  }

  if (multiway == FALSE) {
    return_object <- robust_summary(model_tetrads, robust = FALSE)
    return_object$call <- form
    return(return_object)
  }
}
