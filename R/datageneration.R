#' Generate artificial balanced panel data by simulation
#'
#' The `datageneration()` function generates artificial balanced panel data from
#' the data generating process considered in the Monte Carlo experiment
#' in Yanagi (2023).
#'
#' @param N The number of cross-sectional units
#' @param S The length of time series
#'
#' @returns A data.frame that contains the following elements.
#' \item{id}{The unit index}
#' \item{period}{The time period index}
#' \item{Y}{The outcome}
#' \item{D}{The treatment}
#' \item{X}{The unit-specific covariate}
#' \item{once}{Whether some treatment realization has been received at least once so far}
#' \item{event}{The earliest time period of receiving some treatment realization so far (the event date)}
#' \item{number}{The number of treatment adoptions so far}
#'
#' @examples
#' set.seed(1)
#' data <- datageneration(N = 1000, S = 4)
#'
#' @references Yanagi, T., 2023.
#' An effective treatment approach to
#' difference-in-differences with general treatment patterns.
#' arXiv:2212.13226.
#'
#' @export
#'
datageneration <- function(N, S) {

  # Unit specific covariate
  X <- stats::rnorm(N)

  # Individual effect
  alpha <- stats::rnorm(N)

  # Error terms
  v  <- matrix(stats::rnorm(N * S), nrow = N, ncol = S)
  u  <- matrix(stats::rnorm(N * S), nrow = N, ncol = S)
  xi <- matrix(stats::rnorm(N * S), nrow = N, ncol = S)

  # Parameters
  pi1 <- -1
  pi2 <-  1
  tau <- matrix(NA, nrow = S, ncol = S)
  eta <- gamma <- lambda <- rep(NA, S)
  for (t in 1:S) {
    for (e in 1:S) {
      tau[t, e] <- (t + S - e) / S
    }
    eta[t]    <- t
    gamma[t]  <- t
    lambda[t] <- t / S
  }

  # Binary treatment
  D <- matrix(0, nrow = N, ncol = S)
  for (t in 2:S) {
    D[, t] <- ifelse(pi1 + X * pi2 + alpha + lambda[t] >= u[, t], 1, 0)
  }

  # The once specification
  once_func <- function(D, t) {
    max(D[1:t] != 0)
  }

  # The event specification
  event_func <- function(D, t) {
    if (any(D[1:t] != 0)) {
      min(which(D[1:t] != 0))
    } else {
      0
    }
  }

  # The number specification
  number_func <- function(D, t) {
    sum(D[1:t] != 0)
  }

  # Each realization
  once <- event <- number <- matrix(NA, nrow = N, ncol = S)
  for (t in 1:S) {
    once[, t]   <- apply(D, MARGIN = 1, FUN = once_func,   t = t)
    event[, t]  <- apply(D, MARGIN = 1, FUN = event_func,  t = t)
    number[, t] <- apply(D, MARGIN = 1, FUN = number_func, t = t)
  }

  # Outcome
  Y <- Y0 <- matrix(NA, nrow = N, ncol = S)
  for (t in 1:S) {
    effect <- rep(0, N)
    for (e in 1:t) {
      effect <- effect + tau[t, e] * ifelse(event[, t] == e, 1, 0)
    }
    Y0[, t] <- X * gamma[t] + alpha + eta[t] + v[, t]
    Y[, t]  <- Y0[, t] + effect + xi[, t]
  }

  # Return a data.frame
  data <- NULL
  for (t in 1:S) {
    data <- rbind(data,
                  cbind(1:N,
                        t,
                        Y[, t],
                        D[, t],
                        X,
                        once[, t],
                        event[, t],
                        number[, t]
                  )
    )
  }

  colnames(data) <- c("id", "period", "Y", "D", "X", "once", "event", "number")
  rownames(data) <- NULL

  data <- as.data.frame(data) %>%
    dplyr::arrange(id, period)

  return(data)

}
