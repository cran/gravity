#' @title Structural Iterated Least Squares (SILS)
#'
#' @description \code{sils} estimates gravity models via
#' Structural Iterated Least Squares and an explicit inclusion
#' of the Multilateral Resistance terms.
#'
#' @details \code{sils} is an estimation method for gravity models
#' developed by \insertCite{Head2014;textual}{gravity}.
#'
#' The function \code{sils} utilizes the relationship between the Multilateral
#' Resistance terms and the transaction costs. The parameters are estimated by
#' an iterative procedure. The function executes loops until the parameters
#' stop changing significantly.
#'
#' \code{sils} is designed to be consistent with the Stata code provided at
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' As, to our knowledge at the moment, there is no explicit literature covering
#' the estimation of a gravity equation by \code{sils} using panel data,
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
#' @param maxloop maximum number of outer loop iterations.
#' The default is set to 100. There will be a warning if the iterations
#' did not converge.
#'
#' @param decimal_places number of decimal places that should not change after a new
#' iteration for the estimation to stop. The default is set to 4.
#'
#' @param robust robust (type: logic) determines whether a robust
#' variance-covariance matrix should be used. The default is set to \code{TRUE}.
#'
#' If set \code{TRUE} the estimation results are consistent with the
#' Stata code provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' @param verbose (type: logic) determines whether the estimated coefficients
#' of each iteration should be printed in the console. The default is set
#' to \code{FALSE}.
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
#' @param ... additional arguments to be passed to functions used by
#' \code{sils}.
#'
#' @references
#' For information on \code{sils} as well as more information on gravity
#' models, theoretical foundations and suitable estimation methods in general see
#'
#' \insertRef{Head2014}{gravity}
#'
#' \insertRef{Anderson2001}{gravity}
#'
#' For more information on gravity models, theoretical foundations and
#' estimation methods in general see
#'
#' \insertRef{Anderson1979}{gravity}
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
#' sils(dependent_variable = "flow", regressors = c("distw", "rta"),
#' incomes = c("gdp_o", "gdp_d"),
#' maxloop = 100, dec_places = 4, robust = TRUE, verbose = FALSE,
#' data = gravity_no_zeros)
#'
#' sils(dependent_variable = "flow", regressors = c("distw", "rta", "comcur", "contig"),
#' incomes = c("gdp_o", "gdp_d"),
#' maxloop = 100, dec_places = 4, robust = TRUE, verbose = FALSE,
#' data = gravity_no_zeros)
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
#' sils(dependent_variable = "flow", regressors = c("distw", "rta"),
#' incomes = c("gdp_o", "gdp_d"), maxloop = 100,
#' dec_places = 4, robust = TRUE, verbose = TRUE, data = grav_small)
#' }
#'
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object. It furthermore returns the resulting coefficients for each
#' iteration.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[sandwich]{vcovHC}}
#'
#' @export
#'
sils <- function(dependent_variable, regressors, incomes, maxloop=50, decimal_places=4, robust=TRUE, verbose=FALSE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(incomes) | all(incomes %in% colnames(data)) | length(incomes) == 2)
  stopifnot(maxloop > 0)
  stopifnot(decimal_places >= 1)

  # Split input vectors -----------------------------------------------------
  inc_o <- incomes[1]
  inc_d <- incomes[2]

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

  # Setting starting values for the first iteration ----------------------------
  d <- d %>%
    mutate(
      P_i = 1,
      P_j = 1
    )

  loop <- 0
  dec_point <- 1 * 10^-decimal_places
  beta_distance <- 1
  beta_distance_old <- 0
  coef_dist <- 1

  beta <- vector(length = length(additional_regressors))
  names(beta) <- additional_regressors
  for (j in 1:length(additional_regressors)) {
    beta[j] <- 1
  }
  beta_old <- vector(length = length(additional_regressors))
  names(beta_old) <- additional_regressors

  for (j in 1:length(additional_regressors)) {
    beta_old[j] <- 0
  }

  coef_additional_regressors <- data.frame(matrix(nrow = 1, ncol = length(additional_regressors)))
  coef_additional_regressors[1, ] <- 1
  names(coef_additional_regressors) <- additional_regressors

  # Begin iterations -----------------------------------------------------------
  while (loop <= maxloop &
    abs(beta_distance - beta_distance_old) > dec_point &
    prod(abs(beta - beta_old) > dec_point) == 1) {

    # Updating betas -----------------------------------------------------------
    beta_distance_old <- beta_distance
    for (j in 1:length(additional_regressors)) {
      beta_old[j] <- beta[j]
    }

    # Updating transaction costs -----------------------------------------------
    costs <- data.frame(matrix(nrow = nrow(d), ncol = length(additional_regressors)))
    for (j in 1:length(additional_regressors)) {
      costs[, j] <- beta[j] * d[additional_regressors[j]][, 1]
    }
    costs <- apply(X = costs, MARGIN = 1, FUN = sum)

    d <- d %>%
      ungroup() %>%
      mutate(
        bd = beta_distance,
        co = costs,
        t_ij = exp(!!sym("bd") * !!sym("dist_log") + !!sym("co"))
      )

    # Contraction mapping ------------------------------------------------------
    d <- d %>%
      mutate(
        P_j_old = 0,
        P_i_old = 0
      )

    j <- 1

    while (j <= maxloop &
      sum(abs(d$P_j - d$P_j_old)) > dec_point &
      sum(abs(d$P_i - d$P_i_old)) > dec_point) {
      d <- d %>%
        mutate(
          P_j_old = !!sym("P_j"),
          P_i_old = !!sym("P_i")
        )

      # Inward MR --------------------------------------------------------------
      d <- d %>%
        group_by(!!sym("iso_d")) %>%
        mutate(
          P_j = sum(!!sym("t_ij") * !!sym(inc_o) / !!sym("P_i"))
        )

      # Outward MR -------------------------------------------------------------
      d <- d %>%
        group_by(!!sym("iso_o")) %>%
        mutate(
          P_i = sum(!!sym("t_ij") * !!sym(inc_d) / !!sym("P_j"))
        )

      j <- j + 1

      if (j == maxloop) {
        warning("The inner iteration did not converge before the inner loop reached maxloop=", maxloop, " iterations")
      }
    }

    # Model --------------------------------------------------------------------
    d <- d %>%
      mutate(
        y_log_sils = log(!!sym(dependent_variable)) -
          log((!!sym(inc_o) * !!sym(inc_d)) / (!!sym("P_i") * !!sym("P_j")))
      )

    vars <- paste(c("dist_log", additional_regressors), collapse = " + ")
    form <- stats::as.formula(paste("y_log_sils", "~", vars))

    model_sils <- stats::lm(form, data = d)

    # Updating coefficients ----------------------------------------------------
    beta_distance <- stats::coef(model_sils)[2]
    for (j in 1:length(additional_regressors)) {
      beta[j] <- stats::coef(model_sils)[j + 2]
    }

    coef_dist <- c(coef_dist, beta_distance)
    coef_additional_regressors <- rbind(coef_additional_regressors, rep(0, times = length(additional_regressors)))
    for (j in 1:length(additional_regressors)) {
      coef_additional_regressors[additional_regressors[j]][loop + 2, ] <- beta[j]
    }

    # Coefficients -------------------------------------------------------------
    coef_sils <- cbind(
      loop = c(1:loop), dist = as.numeric(coef_dist)[2:(loop + 1)],
      coef_additional_regressors[2:(loop + 1), ]
    )

    if (verbose == TRUE) {
      cat("This is round", loop, "\n")
      cat("The coefficients are", beta_distance, beta, "\n")
    }

    loop <- loop + 1

    if (loop == maxloop) {
      warning("The outer iteration did not converge before the outer loop reached maxloop=", maxloop, " iterations")
    }
  }

  # Return ---------------------------------------------------------------------
  if (robust == TRUE) {
    return_object <- robust_summary(model_sils, robust = TRUE)
    return_object$call <- form
    return(return_object)
  }

  if (robust == FALSE) {
    return_object <- robust_summary(model_sils, robust = FALSE)
    return_object$call <- form
    return(return_object)
  }
}
