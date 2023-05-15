#' Error handling
#'
#' @param yname The name of the outcome variable
#' @param dname The name of the treatment variable
#' @param tname The name of the time period
#' @param idname The name of the unit index
#' @param xformla A formula for the unit-specific covariates to be included
#' @param data A data.frame of panel data (long format)
#' @param specification A character for the effective treatment specification
#' @param alp The significance level
#' @param nboot The number of bootstrap repetitions
#'
#' @returns NULL if there are no errors
#'
#' @noRd
#'
errorhandling <- function(yname,
                          dname,
                          tname,
                          idname,
                          xformla,
                          data,
                          specification,
                          alp,
                          nboot) {

  # yname, tname, idname, dname ------------------------------------------------

  if (!is.character(yname)  | length(yname)  != 1 |
      !is.character(dname)  | length(dname)  != 1 |
      !is.character(tname)  | length(tname)  != 1 |
      !is.character(idname) | length(idname) != 1 ) {

    # Character
    stop(paste("yname, dname, tname, and idname must be characters."))

  }

  # xformla --------------------------------------------------------------------

  if (!rlang::is_formula(xformla)) {

    # Formula
    stop(paste("xformla must be a formula of the form `xformla = ~ X1 + X2`."))

  }

  # data -----------------------------------------------------------------------

  if (!is.data.frame(data)) {

    # data.frame
    stop(paste("data must be a data.frame."))

  }

  # specification --------------------------------------------------------------

  if (specification != "once"   &
      specification != "event"  &
      specification != "number" &
      specification != "aggregate") {

    # Inappropriate specification
    stop(paste("specification must be one of 'once', 'event', 'number', and 'aggregate'."))

  }

  # alp ------------------------------------------------------------------------

  if (!is.numeric(alp)) {

    # Numeric
    stop(paste("alp must be a positive number between 0 and 1."))

  } else {

    # Scalar
    if (length(alp) != 1) {
      stop(paste("alp must be a positive number between 0 and 1."))
    }

    # Positive value between 0 and 1
    if (alp <= 0 | alp >= 1) {
      stop(paste("alp must be a positive number between 0 and 1."))
    }

  }

  # nboot ----------------------------------------------------------------------

  # Scalar
  if (length(nboot) != 1) {
    stop(paste("nboot must be a positive number."))
  }

  # Positive number
  if (nboot <= 0) {
    stop(paste("nboot must be a positive number."))
  }

  # Return ---------------------------------------------------------------------

  return(NULL)

}
