#' An effective treatment approach to DiD with general treatment patterns
#'
#' The `didet()` function takes the effective treatment approach
#' to difference-in-differences estimation with general treatment patterns.
#' The method is developed by Yanagi (2023).
#' For more details, see the package vignette with `vignette("didet")` and
#' `vignette("review", package = "didet")`.
#'
#' @param yname The name of the outcome variable
#' @param dname The name of the treatment variable
#' @param tname The name of the time period
#' @param idname The name of the unit index
#' @param xformla A formula for the unit-specific covariates to be included.
#' It should be of the form `xformla = ~ X1 + X2`.
#' Default is `xformla = ~ 1`.
#' @param data A data.frame of balanced panel data (long format)
#' @param specification A character specifying the effective treatment function.
#' Options are "once", "event", "number", and "aggregate".
#' Default is "event".
#' @param alp The significance level.
#' Default is 0.05.
#' @param nboot The number of bootstrap repetitions.
#' Default is 1000.
#'
#' @returns A list that contains the following elements.
#' \item{ATEM}{A data.frame that collects results for ATEM(t,s,e)}
#' \item{mover}{A data.frame that collects results for the probability of mover}
#' \item{stayer}{A data.frame that collects results for the probability of stayer}
#' \item{figure}{A list that contains the ggplot2 figures for ATEM(t,s,e)}
#'
#' @examples
#' set.seed(1)
#' data <- datageneration(N = 1000, S = 4)
#' did_est <- didet(yname = "Y",
#'                  dname = "D",
#'                  tname = "period",
#'                  idname = "id",
#'                  xformla = ~ X,
#'                  data = data,
#'                  specification = "event",
#'                  alp = 0.05,
#'                  nboot = 1000)
#'
#' @references Yanagi, T., 2023.
#' An effective treatment approach to
#' difference-in-differences with general treatment patterns.
#' arXiv:2212.13226.
#'
#' @importFrom data.table :=
#' @importFrom dplyr arrange distinct filter group_by mutate pull
#' row_number sym ungroup
#' @importFrom ggplot2 aes guides guide_legend element_line element_rect
#' element_text geom_errorbar geom_hline geom_point
#' ggplot labs position_dodge
#' scale_color_manual scale_shape_manual
#' scale_x_discrete scale_y_continuous
#' theme theme_classic theme_set unit
#' @importFrom magrittr %>%
#' @importFrom rlang UQ
#' @importFrom utils globalVariables
#'
#' @export
#'
didet <- function(yname,
                  dname,
                  tname,
                  idname,
                  xformla = ~ 1,
                  data,
                  specification = "event",
                  alp = 0.05,
                  nboot = 1000) {

  #-----------------------------------------------------------------------------
  # Error handling
  #-----------------------------------------------------------------------------

  error <- errorhandling(yname = yname,
                         dname = dname,
                         tname = tname,
                         idname = idname,
                         xformla = xformla,
                         data = data,
                         specification = specification,
                         alp = alp,
                         nboot = nboot)

  #-----------------------------------------------------------------------------
  # Effective treatment functions
  #-----------------------------------------------------------------------------

  # Arguments of the function defined below
  # D: The T-dimensional treatment vector
  # t: A time period (1 <= t <= T)

  if (specification == "once" | specification == "aggregate") {

    et_func <- function(D, t) {
      max(D[1:t] != 0)
    }

  } else if (specification == "event") {

    et_func <- function(D, t) {
      if (any(D[1:t] != 0)) {
        min(which(D[1:t] != 0))
      } else {
        0
      }
    }

  } else if (specification == "number") {

    et_func <- function(D, t) {
      sum(D[1:t] != 0)
    }

  }

  #-----------------------------------------------------------------------------
  # Functions for movers and stayers
  #-----------------------------------------------------------------------------

  # Arguments of the functions defined below
  # E: The T-dimensional effective treatment vector
  # t: A time period (1 <= t <= T)
  # s: A time period (s < t)
  # e: An effective treatment intensity

  mover_func <- function(E, t, s, e) {
    ifelse(E[t] == e & E[s] == 0, 1, 0)
  }

  stayer_func <- function(E, t, s) {
    ifelse(E[t] == 0 & E[s] == 0, 1, 0)
  }

  #-----------------------------------------------------------------------------
  # Variable definitions
  #-----------------------------------------------------------------------------

  # Arrange
  data <- data %>%
    arrange(UQ(sym(idname)), UQ(sym(tname)))

  # Original time periods
  origin_period <- data %>%
    pull(UQ(sym(tname))) %>%
    unique()

  # Add unit index (i = 1, ..., N) and time index (t = 1, ..., T)
  data <- data %>%
    group_by(UQ(sym(tname))) %>%
    mutate(unit = row_number(UQ(sym(idname)))) %>%
    ungroup() %>%
    group_by(unit) %>%
    mutate(period = row_number(UQ(sym(tname)))) %>%
    ungroup() %>%
    arrange(unit, period)

  # The number of cross-sectional units
  N <- data$unit %>%
    unique() %>%
    length()

  # The length of the time series
  S <- data$period %>%
    unique() %>%
    length()

  # Add the outcome, treatment, and effective treatment
  data <- data %>%
    mutate(Y = UQ(sym(yname)),
           D = UQ(sym(dname))) %>%
    group_by(unit) %>%
    mutate(E = purrr::map_dbl(.x = 1:S,
                              .f = et_func,
                              D = D)) %>%
    ungroup()

  #-----------------------------------------------------------------------------
  # Set of (t, s, r, e) at which we estimate ATEM(t, s, r, e)
  # Definition: ATEM(t, s, r, e) = E[ Y_{ir} - Y_{ir}^*(0) | M_{i,t,s,e} = 1 ]
  #-----------------------------------------------------------------------------

  # A1: Set of (t, s, r, e)
  A1 <- NULL

  if (specification == "once" | specification == "aggregate") {

    # Effective treatment intensity
    e <- 1

    for (t in 2:S) {

      # Periods s and r
      s <- 1
      r <- NA

      # Add the indicators for movers and stayers to data
      data <- data %>%
        group_by(unit) %>%
        mutate(ind_mover  =  mover_func(E = E, t = t, s = s, e = e),
               ind_stayer = stayer_func(E = E, t = t, s = s)) %>%
        ungroup() %>%
        mutate("M_{t}_{s}_{e}" := ind_mover,
               "S_{t}_{s}"     := ind_stayer)

      # Message for the overlap condition
      if (sum(data$ind_mover) == 0 | sum(data$ind_stayer) == 0) {
        message(paste0("Movers or stayers do not exist for t = ",
                       origin_period[t], ", s = ", origin_period[s],
                       ", e = ", e, "."))

        next
      }

      # Construct A1
      A1 <- rbind(A1, c(t, s, r, e))

    }

  } else if (specification == "event") {

    # Support of event date, except for "0" and "1"
    supp_e <- data$E %>%
      unique() %>%
      sort() %>%
      setdiff(0:1)

    for (e in supp_e) {

      # Period s
      s <- e - 1

      for (t in 2:S) {

        # Period r
        if (t >= e) {
          r <- NA
        } else {
          r <- t
          t <- e
        }

        # Add the indicators for movers and stayers to data
        data <- data %>%
          group_by(unit) %>%
          mutate(ind_mover  =  mover_func(E = E, t = t, s = s, e = e),
                 ind_stayer = stayer_func(E = E, t = t, s = s)) %>%
          ungroup() %>%
          mutate("M_{t}_{s}_{e}" := ind_mover,
                 "S_{t}_{s}"     := ind_stayer)

        # Message for the overlap condition
        if (sum(data$ind_mover) == 0 | sum(data$ind_stayer) == 0) {
          message(paste0("Movers or stayers do not exist for t = ",
                         origin_period[t], ", s = ", origin_period[s],
                         ", e = ", e, "."))

          next
        }

        # Construct A1
        A1 <- rbind(A1, cbind(t, s, r, e))

      }
    }

  } else if (specification == "number") {

    # Support of number, except for "0"
    supp_e <- data$E %>%
      unique() %>%
      sort() %>%
      setdiff(0)

    # Periods s and r
    s <- 1
    r <- NA

    for (e in supp_e) {

      if (e == S) {
        break
      }

      for (t in (e + 1):S) {

        # Add the indicators for movers and stayers to data
        data <- data %>%
          group_by(unit) %>%
          mutate(ind_mover  =  mover_func(E = E, t = t, s = s, e = e),
                 ind_stayer = stayer_func(E = E, t = t, s = s)) %>%
          ungroup() %>%
          mutate("M_{t}_{s}_{e}" := ind_mover,
                 "S_{t}_{s}"     := ind_stayer)

        # Message for the overlap condition
        if (sum(data$ind_mover) == 0 | sum(data$ind_stayer) == 0) {
          message(paste0("Movers or stayers do not exist for t = ",
                         origin_period[t], ", s = ", origin_period[s],
                         ", e = ", e, "."))

          next
        }

        # Construct A1
        A1 <- rbind(A1, c(t, s, r, e))

      }
    }

  }

  # Name
  colnames(A1) <- c("t", "s", "r", "e")

  # Convert to data.frame, arrange, and unique rows
  A1 <- as.data.frame(A1) %>%
    arrange(e, t, s, r) %>%
    distinct()

  #-----------------------------------------------------------------------------
  # DiD estimation of ATEM(t, s, r, e)
  #-----------------------------------------------------------------------------

  # For the ATEM estimates
  ATEM_t_s_r_e <- NULL

  # For the influence functions
  inf_func_ATEM_i_t_s_r_e <- NULL

  # A2: Resulting set of (t, s, r, e)
  A2 <- NULL

  for (a in 1:nrow(A1)) {

    # Specify (t, s, r, e)
    t <- A1$t[a]
    s <- A1$s[a]
    r <- A1$r[a]
    e <- A1$e[a]

    # Data for movers and stayers
    var1 <- sym(paste0("M_", t, "_", s, "_", e))
    var2 <- sym(paste0("S_", t, "_", s))
    DiDdata <- data %>%
      filter(UQ(var1) == 1 | UQ(var2) == 1) %>%
      mutate(mover = UQ(var1))

    # Unit IDs used for the DiD estimation
    unit_id <- unique(DiDdata$unit)

    # Outcomes for the DRDID estimation
    if (is.na(r)) {

      y1 <- DiDdata %>%
        filter(period == t) %>%
        pull(Y)

      y0 <- DiDdata %>%
        filter(period == s) %>%
        pull(Y)

    } else {

      y1 <- DiDdata %>%
        filter(period == r) %>%
        pull(Y)

      y0 <- DiDdata %>%
        filter(period == r - 1) %>%
        pull(Y)

    }

    # Treatment
    D <- DiDdata %>%
      filter(period == t) %>%
      pull(mover)

    # Matrix of covariates
    DiDdata0 <- DiDdata %>%
      filter(period == 1)

    covariates <- stats::model.matrix(xformla, data = DiDdata0)

    # Estimate ATEM(t, s, e)
    DiDest <- try(
      DRDID::drdid_panel(y1 = y1,
                         y0 = y0,
                         D = D,
                         covariates = covariates,
                         inffunc = TRUE)
    )

    if (methods::is(DiDest) != "drdid") {

      next

    } else {

      # The ATEM estimate
      ATEM_t_s_r_e <- c(ATEM_t_s_r_e, DiDest$ATT)

      # Influence function
      inf_func_temp <- rep(0, N)
      inf_func_temp[unit_id]  <- DiDest$att.inf.func * N / length(unit_id)
      inf_func_ATEM_i_t_s_r_e <- cbind(inf_func_ATEM_i_t_s_r_e,
                                       inf_func_temp)

      # Construct A2: Set of (t, s, r, e)
      A2 <- rbind(A2, c(t, s, r, e))

    }
  }

  # Name
  colnames(A2) <- c("t", "s", "r", "e")

  # Convert to data.frame and arrange
  A2 <- as.data.frame(A2) %>%
    arrange(e, t, s, r)

  #-----------------------------------------------------------------------------
  # End if 'specification == "aggregate"'
  #-----------------------------------------------------------------------------

  if (specification == "aggregate") {

    # Aggregated ATEM estimate
    aggATEM <- mean(ATEM_t_s_r_e)

    # Influence function for aggregated ATEM estimation
    inf_func_aggATEM_i <- rowMeans(inf_func_ATEM_i_t_s_r_e)

    # Bootstrap weight (Mammen)
    kappa <- (sqrt(5) + 1) / 2

    V_i_b <- sample(x = c(1 - kappa, kappa),
                    size = N * nboot,
                    prob = c(kappa/sqrt(5), 1 - kappa/sqrt(5)),
                    replace = TRUE) %>%
      matrix(nrow = N, ncol = nboot)

    # Bootstrap aggregated ATEM estimates and residuals
    aggATEM_b <- aggATEM + t(inf_func_aggATEM_i) %*% V_i_b / N

    R_aggATEM_b <- aggATEM_b - aggATEM

    # Bootstrap standard error
    SE_aggATEM <- apply(R_aggATEM_b,
                        MARGIN = 1,
                        FUN = stats::IQR,
                        na.rm = TRUE) /
      (stats::qnorm(0.75) - stats::qnorm(0.25))

    # Bootstrapped (1 - alp) critical value
    max_t_stat_aggATEM_b <- apply(abs(R_aggATEM_b) / SE_aggATEM,
                                  MARGIN = 2,
                                  FUN = max,
                                  na.rm = TRUE)

    c_aggATEM <- stats::quantile(max_t_stat_aggATEM_b,
                                 probs = 1 - alp)

    # Bootstrap (1 - alp) confidence interval
    CIL_aggATEM <- aggATEM - c_aggATEM * SE_aggATEM
    CIU_aggATEM <- aggATEM + c_aggATEM * SE_aggATEM

    # Return
    ATEM <- data.frame(t   = NA,
                       s   = NA,
                       r   = NA,
                       e   = 1,
                       est = aggATEM,
                       SE  = SE_aggATEM,
                       CIL = CIL_aggATEM,
                       CIU = CIU_aggATEM)
    rownames(ATEM) <- NULL

    return(list(ATEM   = ATEM,
                mover  = NA,
                stayer = NA,
                figure = NA))

  }

  #-----------------------------------------------------------------------------
  # Estimation of the probabilities for movers and stayers
  #-----------------------------------------------------------------------------

  # For the probabilities and standard errors
  Pr_M_t_s_e <- SE_Pr_M_t_s_e <- NULL
  Pr_S_t_s   <- SE_Pr_S_t_s   <- NULL

  for (a in 1:nrow(A2)) {

    # Periods (t, s, e)
    t <- A2$t[a]
    s <- A2$s[a]
    e <- A2$e[a]

    # Variable definitions
    var1 <- paste0("M_", t, "_", s, "_", e)
    var2 <- paste0("S_", t, "_", s)

    mover <- data %>%
      filter(period == t) %>%
      pull(UQ(var1))

    stayer <- data %>%
      filter(period == t) %>%
      pull(UQ(var2))

    # Estimated probabilities
    Pr_M_t_s_e <- c(Pr_M_t_s_e, mean(mover))
    Pr_S_t_s   <- c(Pr_S_t_s,   mean(stayer))

    # Standard errors
    SE_Pr_M_t_s_e <- c(SE_Pr_M_t_s_e, sqrt(mean(mover)  * mean(1 - mover)  / N))
    SE_Pr_S_t_s   <- c(SE_Pr_S_t_s,   sqrt(mean(stayer) * mean(1 - stayer) / N))

  }

  # (1 - alp) point-wise confidence intervals
  CIL_Pr_M_t_s_e <- Pr_M_t_s_e + stats::qnorm(alp / 2) * SE_Pr_M_t_s_e
  CIU_Pr_M_t_s_e <- Pr_M_t_s_e - stats::qnorm(alp / 2) * SE_Pr_M_t_s_e

  CIL_Pr_S_t_s   <- Pr_S_t_s   + stats::qnorm(alp / 2) * SE_Pr_S_t_s
  CIU_Pr_S_t_s   <- Pr_S_t_s   - stats::qnorm(alp / 2) * SE_Pr_S_t_s

  #-----------------------------------------------------------------------------
  # Multiplier bootstrap inference
  #-----------------------------------------------------------------------------

  # Bootstrap weight (Mammen)
  kappa <- (sqrt(5) + 1) / 2

  V_i_b <- sample(x = c(1 - kappa, kappa),
                  size = N * nboot,
                  prob = c(kappa/sqrt(5), 1 - kappa/sqrt(5)),
                  replace = TRUE) %>%
    matrix(nrow = N, ncol = nboot)

  # Bootstrap ATEM estimates and residuals
  ATEM_t_s_r_e_b <- ATEM_t_s_r_e + t(inf_func_ATEM_i_t_s_r_e) %*% V_i_b / N

  R_ATEM_t_s_r_e_b <- ATEM_t_s_r_e_b - ATEM_t_s_r_e

  # Bootstrap standard error
  SE_ATEM_t_s_r_e <- apply(R_ATEM_t_s_r_e_b,
                           MARGIN = 1,
                           FUN = stats::IQR,
                           na.rm = TRUE) /
    (stats::qnorm(0.75) - stats::qnorm(0.25))

  # Bootstrapped (1 - alp) critical value
  max_t_stat_ATEM_t_s_r_e_b <- apply(abs(R_ATEM_t_s_r_e_b) / SE_ATEM_t_s_r_e,
                                     MARGIN = 2,
                                     FUN = max,
                                     na.rm = TRUE)

  c_ATEM_t_s_r_e <- stats::quantile(max_t_stat_ATEM_t_s_r_e_b,
                                    probs = 1 - alp)

  # Bootstrap (1 - alp) uniform confidence band
  CIL_ATEM_t_s_r_e <- ATEM_t_s_r_e - c_ATEM_t_s_r_e * SE_ATEM_t_s_r_e
  CIU_ATEM_t_s_r_e <- ATEM_t_s_r_e + c_ATEM_t_s_r_e * SE_ATEM_t_s_r_e

  #-----------------------------------------------------------------------------
  # Collect the estimation results
  #-----------------------------------------------------------------------------

  # Modifications
  A2$t <- origin_period[A2$t]
  A2$s <- origin_period[A2$s]
  A2$r <- origin_period[A2$r]

  if (specification == "event") {
    A2$e <- origin_period[A2$e]
  }

  # The ATEM results
  ATEM <- data.frame(t   = A2$t,
                     s   = A2$s,
                     r   = A2$r,
                     e   = A2$e,
                     est = ATEM_t_s_r_e,
                     SE  = SE_ATEM_t_s_r_e,
                     CIL = CIL_ATEM_t_s_r_e,
                     CIU = CIU_ATEM_t_s_r_e) %>%
    distinct()

  # The mover results
  mover <- data.frame(t   = A2$t,
                      s   = A2$s,
                      e   = A2$e,
                      est = Pr_M_t_s_e,
                      SE  = SE_Pr_M_t_s_e,
                      CIL = CIL_Pr_M_t_s_e,
                      CIU = CIU_Pr_M_t_s_e) %>%
    distinct()

  # The stayer results
  stayer <- data.frame(t   = A2$t,
                       s   = A2$s,
                       est = Pr_S_t_s,
                       SE  = SE_Pr_S_t_s,
                       CIL = CIL_Pr_S_t_s,
                       CIU = CIU_Pr_S_t_s) %>%
    distinct()

  #-----------------------------------------------------------------------------
  # Make figures using ggplot 2
  #-----------------------------------------------------------------------------

  # Themes used for ggplot2
  theme_set(
    theme_classic() +
      theme(
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(angle = 0, vjust = 1),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "grey"),
        plot.title = element_text(hjust = 0),
        text = element_text(size = 20)
      )
  )

  # A list that contains the figures
  figure <- list()

  # Support of e, except for "0"
  supp_e <- ATEM$e %>%
    unique() %>%
    sort()

  # Support of t, except for "1"
  supp_t <- ATEM$t %>%
    unique() %>%
    sort() %>%
    as.character()

  # Limits of y-axis
  min_y <- min(ATEM$CIL, -1)
  max_y <- max(ATEM$CIU,  1)

  # Make ATEM figure for each effective treatment intensity
  for (e0 in supp_e) {

    if (specification == "once") {

      # A data.frame used for ggplot2
      df_ATEM <- ATEM %>%
        filter(e == e0) %>%
        mutate(t = as.character(t))

      # Title
      mytitle <- "The once specification"

      # Make a figure
      fig <- ggplot(df_ATEM, aes(x = t, y = est, color = "#F8766D")) +
        theme(legend.position = "none")

    } else if (specification == "event") {

      # A data.frame used for ggplot2
      df_ATEM <- ATEM %>%
        filter(e == e0) %>%
        mutate(prepost = ifelse(is.na(r), "Post", "Pre")) %>%
        mutate(t = ifelse(is.na(r), t, r)) %>%
        mutate(t = as.character(t))

      # Title
      mytitle <- paste0("The event specification (event date = ", e0, ")")

      # Make a figure
      fig <- ggplot(df_ATEM, aes(x = t,
                                 y = est,
                                 color    = prepost,
                                 shape    = prepost)) +
        scale_color_manual(name = "",
                           breaks = c("Pre", "Post"),
                           values = c("#00BFC4", "#F8766D"),
                           limits = c("Pre", "Post")) +
        scale_shape_manual(name = "",
                           breaks = c("Pre", "Post"),
                           values = c(15, 16),
                           limits = c("Pre", "Post")) +
        theme(legend.key.width  = unit(2, "cm")) +
        guides(shape = guide_legend(override.aes = list(size = 7)))

    } else if (specification == "number") {

      # A data.frame used for ggplot2
      df_ATEM <- ATEM %>%
        filter(e == e0) %>%
        mutate(t = as.character(t))

      # Title
      mytitle <- paste0("The number specification (number = ", e0, ")")

      # Make a figure
      fig <- ggplot(df_ATEM, aes(x = t, y = est, color = "#F8766D")) +
        theme(legend.position = "none")

    }

    # Modifications
    fig <- fig +
      labs(
        title = mytitle,
        x = "Period",
        y = "",
        color = ""
      ) +
      geom_errorbar(
        aes(
          ymin = est - est + CIL,
          ymax = est - est + CIU
        ),
        position = position_dodge(0.5),
        width = 0.1,
        linewidth = 2
      ) +
      geom_point(
        position = position_dodge(0.5),
        size = 6
      ) +
      scale_x_discrete(
        breaks = supp_t,
        limits = supp_t
      ) +
      scale_y_continuous(
        limits = c(min_y, max_y)
      )  +
      geom_hline(
        yintercept = 0,
        color = "black",
        linetype = "dotted",
        linewidth = 1.5
      )

    # Record
    figure[[paste0("e", e0)]] <- fig

  }

  #-----------------------------------------------------------------------------
  # return
  #-----------------------------------------------------------------------------

  return(list(ATEM   = ATEM,
              mover  = mover,
              stayer = stayer,
              figure = figure))

}
