# Artificial data --------------------------------------------------------------

set.seed(1)
data <- didet::datageneration(N = 3000, S = 4)

# Variable Definitions ---------------------------------------------------------

yname <- "Y"
tname <- "period"
idname <- "id"
dname <- "D"
xformla <- ~ X
specification <- "event"
alp <- 0.05
nboot <- 1000

# yname, dname, tname, idname --------------------------------------------------

expect_error(didet::didet(yname = TRUE,
                          tname = "period",
                          idname = "id",
                          dname = "D",
                          xformla = ~ X,
                          data = data,
                          specification = "event",
                          alp = 0.05,
                          nboot = 1000),
             "yname, dname, tname, and idname must be characters.")

expect_error(didet::didet(yname = "Y",
                          tname = NA,
                          idname = "id",
                          dname = "D",
                          xformla = ~ X,
                          data = data,
                          specification = "event",
                          alp = 0.05,
                          nboot = 1000),
             "yname, dname, tname, and idname must be characters.")

expect_error(didet::didet(yname = "Y",
                          tname = "period",
                          idname = 10,
                          dname = "D",
                          xformla = ~ X,
                          data = data,
                          specification = "event",
                          alp = 0.05,
                          nboot = 1000),
             "yname, dname, tname, and idname must be characters.")

expect_error(didet::didet(yname = "Y",
                          tname = "period",
                          idname = "id",
                          dname = Inf,
                          xformla = ~ X,
                          data = data,
                          specification = "event",
                          alp = 0.05,
                          nboot = 1000),
             "yname, dname, tname, and idname must be characters.")

# data -------------------------------------------------------------------------

expect_error(didet::didet(yname = "Y",
                          tname = "period",
                          idname = "id",
                          dname = "D",
                          xformla = ~ X,
                          data = as.matrix(data),
                          specification = "event",
                          alp = 0.05,
                          nboot = 1000),
             "data must be a data.frame.")

# specification ----------------------------------------------------------------

expect_error(didet::didet(yname = "Y",
                          tname = "period",
                          idname = "id",
                          dname = "D",
                          xformla = ~ X,
                          data = data,
                          specification = "DID",
                          alp = 0.05,
                          nboot = 1000),
             "specification must be one of 'once', 'event', 'number', and 'aggregate'.")

# alp --------------------------------------------------------------------------

expect_error(didet::didet(yname = "Y",
                          tname = "period",
                          idname = "id",
                          dname = "D",
                          xformla = ~ X,
                          data = data,
                          specification = "event",
                          alp = 5,
                          nboot = 1000),
             "alp must be a positive number between 0 and 1.")

# nboot ------------------------------------------------------------------------

expect_error(didet::didet(yname = "Y",
                          tname = "period",
                          idname = "id",
                          dname = "D",
                          xformla = ~ X,
                          data = data,
                          specification = "event",
                          alp = 0.05,
                          nboot = FALSE),
             "nboot must be a positive number.")
