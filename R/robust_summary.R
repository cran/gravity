#' @title Robust summary statistics
#'
#' @description Internal function for nice summary of estimates
#' adapted after a function written by Isidore Beautrelet
#' \url{https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R}
#'
#' @param object \code{lm} object
#'
#' @param robust robust or non-robust estimation
#'
#' @return summary \code{lm} object
#'
#' @keywords internal
#'
#' @noRd

robust_summary <- function(object, robust = FALSE, ...) {
  # Check -------------------------------------------------------------------
  qr_lm <- function(x, ...) {
    if (is.null(r <- x$qr)) {
      stop("lm object does not have a proper 'qr' component.\n 
           Rank zero or should not have used lm(.., qr=FALSE).")
    }
    r
  }

  # Robust standard errors --------------------------------------------------
  if (robust == TRUE) {
    # Save variables that are necessary to calculate robust sd
    X <- stats::model.matrix(object)
    u2 <- stats::residuals(object)^2
    XDX <- t(X) %*% (X * u2)

    # Inverse(X'X)
    XX1 <- solve(t(X) %*% X, tol = 1e-100)

    # Sandwich Variance calculation (Bread x meat x Bread)
    varcovar <- XX1 %*% XDX %*% XX1

    # Adjust degrees of freedom
    dfc_r <- sqrt(nrow(X)) / sqrt(nrow(X) - ncol(X))

    # Standard errors of the coefficient estimates are the
    # square roots of the diagonal elements
    rstdh <- dfc_r * sqrt(diag(varcovar))
  }

  # Clustered standard errors -----------------------------------------------
  z <- object
  p <- z$rank
  rdf <- z$df.residual

  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss / rdf
    ans <- z[c(
      "call", "terms",
      if (!is.null(z$weights)) "weights"
    )]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(stats::coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(
      NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    )
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }

  if (is.null(z$terms)) {
    stop("invalid 'lm' object:  no 'terms' component")
  }

  if (!inherits(object, "lm")) {
    warning("calling summary.lm(<fake-lm-object>) ...")
  }

  Qr <- qr_lm(object)
  n <- NROW(Qr$qr)

  if (is.na(z$df.residual) || n - p != z$df.residual) {
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  }

  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights

  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept")) {
      sum((f - mean(f))^2)
    } else {
      sum(f^2)
    }
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f / sum(w))
      sum(w * (f - m)^2)
    }
    else {
      sum(w * f^2)
    }
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }

  resvar <- rss / rdf

  if (is.finite(resvar) && resvar < (mean(f)^2 + stats::var(f)) * 1e-30) {
    warning("essentially perfect fit: summary may be unreliable")
  }

  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)

  if (robust == TRUE) se <- rstdh

  est <- z$coefficients[Qr$pivot[p1]]
  tval <- est / se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  pval <- 2 * stats::pt(abs(tval), rdf, lower.tail = FALSE)
  ans$coefficients <- cbind(est, se, tval, pval)
  dimnames(ans$coefficients) <- list(
    names(z$coefficients)[Qr$pivot[p1]],
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  )
  ans$aliased <- is.na(stats::coef(object))
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf, NCOL(Qr$qr))

  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) {
      1L
    } else {
      0L
    }
    ans$r.squared <- mss / (mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int) / rdf)
    ans$fstatistic <- c(
      value = (mss / (p - df.int)) / resvar,
      numdf = p - df.int, dendf = rdf
    )
    if (robust == TRUE) {
      if (df.int == 0) {
        pos_coef <- 1:length(z$coefficients)
      } else {
        pos_coef <- match(
          names(z$coefficients)[-match(
            "(Intercept)",
            names(z$coefficients)
          )],
          names(z$coefficients)
        )
      }
      P_m <- matrix(z$coefficients[pos_coef])

      R_m <- diag(
        1,
        length(pos_coef),
        length(pos_coef)
      )

      ans$fstatistic <- c(
        value = t(R_m %*% P_m) %*%
          (solve(varcovar[pos_coef, pos_coef], tol = 1e-100)) %*%
          (R_m %*% P_m) / (p - df.int),
        numdf = p - df.int, dendf = rdf
      )
    }
  }
  else {
    ans$r.squared <- 0
    ans$adj.r.squared <- 0
  }

  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 1)]

  if (!is.null(z$na.action)) {
    ans$na.action <- z$na.action
  }

  # Output ------------------------------------------------------------------
  class(ans) <- "summary.lm"
  ans
}
