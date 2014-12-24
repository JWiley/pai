#' Create a varying index
#'
#' This is a simple utility to create a fixed length index.
#'
#' @param i Where the index should start
#' @param end Where the index should end (essentially length)
#' @return A vector of the index
#' @export
#' @keywords utils
#' @examples
#'
#' indexer(1, 10)
#' indexer(2, 10)
#' indexer(3, 5)
indexer <- function(i, end) {
  ((i - 1) * end + 1):(i * end)
}

#' Calculate Personalized Advantage Index (PAI)
#'
#' Calculates the PAI for each observation in the dataset using simple
#' linear regression, but with cross-validation and optionally bootstrapping
#' in order to construct confidence intervals.
#'
#' @param formula An R formula, but without dependent variable
#' @param DV a character string indicating the name of the dependent variable
#' @param TreatVar a character string indicating the name of the treatment variable
#'   (note that values must be coded as 0/1 for the treatment variable.
#' @param dat A dataset to use
#' @param boot A logical value whether to bootstrap or not
#' @param k The number of bootstrap resamples to take
#' @param holdouts The holdouts to be used.  If a vector from 1 to the number of
#'   cases available in the dataset, it is equivalent to leave-one-out.
#'   Must be constructed manually.
#' @param yhigherisbetter A logical value whether higher scores indicate better
#'   treatment response on the dependent variable
#' @return A three dimensional array with cases on the first dimension,
#'   multiple estimates calculated on the second dimension, and bootstrap resamples
#'   on the third dimension.
#' @import Rcpp RcppEigen
#' @export
#' @examples
#' \dontrun{
#' # builtin dataset with 32 cases
#' dat <- mtcars
#' dat$cyl <- factor(dat$cyl)
#' formula <- ~ cyl + am * (mpg + hp + drat + wt)
#' DV <- "disp"
#' TreatVar <- "am"
#' k <- 50 # use a few just for speed
#' holdouts <- 1:32 # for leave-one-out
#'
#' # for 10 fold cross-validation
#' holdouts <- by(sample(1:32), cut(1:32, ceiling(seq(from = 0, to = 32, length.out = 11))), function(x) x)
#' # remove names
#' names(holdouts) <- NULL
#'
#' m <- lm_pai_cvboot_fit(formula, DV, TreatVar, dat, boot = TRUE, k = k, holdouts = holdouts)
#'
#' # see result from first bootstrap resample
#' m[, , 1]
#' }
lm_pai_cvboot_fit <- function(formula, DV, TreatVar, dat, boot = TRUE, k, holdouts, yhigherisbetter = TRUE) {

  stopifnot(all(unique(dat[, TreatVar]) %in% c(0, 1)))

  res <- lapply(holdouts, function(n) {
    nholdout <- length(n)
    ntrain <- nrow(dat) - nholdout

    if (boot) {
      boot.index <- sample(1:ntrain, size = ntrain * k, replace = TRUE)
    } else {
      boot.index <- 1:ntrain
      k <- 1
    }

    X <- model.matrix(formula, data = dat[-n, ])
    Y <- dat[-n, DV]

    preddat <- dat[c(n, n, n), ]
    # first set is for zero, second set is for one, third set is for true
    preddat[, TreatVar] <- c(rep(c(0, 1), each = nholdout), dat[n, TreatVar])

    predX <- model.matrix(formula, data = preddat)

    lapply(lapply(1:k, indexer, end = ntrain), function(i) {
      b <- fastLmPure(X[boot.index[i], ], Y[boot.index[i]])$coefficients
      yhat <- as.vector(predX %*% b)
      out <- cbind(Y = dat[n, DV],
                   Treatment = dat[n, TreatVar],
                   Yhat.reality = yhat[(2 * nholdout + 1):(3 * nholdout)],
                   Yhat.0 = yhat[1:nholdout],
                   Yhat.1 = yhat[(nholdout + 1):(2 * nholdout)])
      colnames(out)[2] <- TreatVar
      out <- cbind(out,
                   Yhat.PAI = out[, "Yhat.0"] - out[, "Yhat.1"],
                   Residual = out[, "Y"] - out[, "Yhat.reality"])

      if (yhigherisbetter) {
        out <- cbind(out, Lucky = as.integer(ifelse(out[, "Yhat.0"] < out[, "Yhat.1"], 1, 0) == out[, TreatVar]))
      } else {
        out <- cbind(out, Lucky = as.integer(ifelse(out[, "Yhat.0"] < out[, "Yhat.1"], 0, 1) == out[, TreatVar]))
      }

      ties <- which(out[, "Yhat.PAI"] == 0)
      if (length(ties)) {
        out[ties, "Lucky"] <- 9
      }
      return(out)
    })
  })

  bigres <- array(NA_real_, dim = c(nrow(dat), 8L, k),
    dimnames = list(NULL, colnames(res[[1]][[1]]), NULL))

  for (i in 1:k) {
    for (j in 1:length(holdouts)) {
      bigres[holdouts[[j]], ,i] <- as.vector(res[[j]][[i]])
    }
  }

  class(bigres) <- c("pai", "array")

  return(bigres)
}

#' Calculate PAI on multiply imputed data in parallel
#'
#' This is a simple function that wraps \code{lm_pai_cvboot}.
#'
#' @param formula An R formula, but without dependent variable
#' @param DV a character string indicating the name of the dependent variable
#' @param TreatVar a character string indicating the name of the treatment variable
#'   (note that values must be coded as 0/1 for the treatment variable.
#' @param data A list of imputed datasets or for a single dataset, a list of one
#' @param nboot The number of bootstrap resamples to take
#' @param holdouts The holdouts to be used.  This should be a numeric list with as many elements
#'   desired for holdouts.  For example, if the list is the same length
#'   as the number of cases and each element is one case (i.e., 1, 2, \ldots, N),
#'   then it is equivalent to leave-one-out.  If there are ten elements, each containing
#'   approximately 1/10th of the sample, it is 10-fold cross validation, etc.
#'   If a single character string is passed, K-fold CV is assumed and holdouts
#'   are automatically constructed.
#' @param cl An existing cluster to use (optional)
#' @param cores The number of cores to use when creating a cluster, if an existing
#'   cluster is not passed (optional).  If left blank and no cluster passed,
#'   defaults to the number of cores available.
#' @param seed A seed to make the bootstrap reproducible
#' @param boot A logical value whether to bootstrap or not
#' @param yhigherisbetter A logical value whether higher scores indicate better
#'   treatment response on the dependent variable
#' @return A list containing a four dimensional array of the results, as well as
#'   the input parameters for the model.  For the array, the first dimension is cases,
#'   the second dimension is multiple estimates calculated, the third dimension is
#'   bootstrap resamples, and the fourth dimension is imputations.
#' @export
#' @import mice
#' @examples
#' \dontrun{
#' # builtin dataset with 32 cases
#' dat <- mtcars
#' dat$cyl <- factor(dat$cyl)
#' formula <- ~ cyl + am * (mpg + hp + drat + wt)
#' DV <- "disp"
#' TreatVar <- "am"
#' k <- 50 # use a few just for speed
#' holdouts <- 1:32 # for leave-one-out
#'
#' m <- mi_lm_pai_cvboot(formula, DV, TreatVar, list(dat, dat),
#'   nboot = k, holdouts = holdouts, cores = 2)
#' }
mi_lm_pai_cvboot <- function(formula, DV, TreatVar, data, nboot = 500, holdouts,
                             cl, cores, seed = 1234, boot = TRUE, yhigherisbetter = TRUE) {

  N <- length(data)
  ncases <- nrow(data[[1]])

  if (length(holdouts) == 1 && is.character(holdouts)) {
    set.seed(seed)
    holdouts <- by(sample(1:ncases), cut(1:ncases,
      ceiling(seq(from = 0, to = ncases, length.out = as.numeric(holdouts) + 1))),
      function(x) x)
  }

  if (missing(cl)) {
    if (missing(cores)) {
      cores <- detectCores()
    }
    cl <- makeCluster(cores)
    on.exit(stopCluster(cl))
  }

  clusterEvalQ(cl, {
    library(pai)
  })

  env <- environment()
  clusterExport(cl, c("formula", "DV", "TreatVar", "data", "nboot", "holdouts", "yhigherisbetter"),
    envir = env)

  clusterSetRNGStream(cl, seed)

  out <- parLapply(cl, 1:N, function(i) {
    lm_pai_cvboot_fit(formula = formula,
      DV = DV, TreatVar = TreatVar, dat = data[[i]],
      boot = boot, k = nboot, holdouts = holdouts,
      yhigherisbetter = yhigherisbetter)
  })

  finalout <- array(NA_real_, c(dim(out[[1]]), N),
    dimnames = list(NULL, colnames(out[[1]]), NULL, NULL))

  for (i in 1:N) {
    finalout[, , , i] <- out[[i]]
  }

  if (boot) {
    nobootout <- parLapply(cl, 1:N, function(i) {
    lm_pai_cvboot_fit(formula = formula,
      DV = DV, TreatVar = TreatVar, dat = data[[i]],
      boot = FALSE, k = 1, holdouts = holdouts,
      yhigherisbetter = yhigherisbetter)
    })

    nobootfinalout <- array(NA_real_, c(dim(nobootout[[1]]), N),
    dimnames = list(NULL, colnames(nobootout[[1]]), NULL, NULL))

    for (i in 1:N) {
      nobootfinalout[, , , i] <- nobootout[[i]]
    }

    res <- list(bootresults = finalout, nobootresults = nobootfinalout,
                DV = DV, TreatVar = TreatVar, nboot = nboot,
                holdouts = holdouts, seed = seed, yhigherisbetter = yhigherisbetter)

    class(res) <- c("mibootpai", "mipai", "list")

  } else {
    nobootout <- out
    nobootfinalout <- finalout

    res <- list(nobootresults = nobootfinalout,
                DV = DV, TreatVar = TreatVar, nboot = nboot,
                holdouts = holdouts, seed = seed, yhigherisbetter = yhigherisbetter)

    class(res) <- c("mipai", "list")
  }

  return(res)
}

#' Calculate Summary Statistics for PAI Scores
#'
#' @param object an mipai or mibootpai class object returned from
#'   \code{mi_lm_pai_cvboot}.
#' @param plot logical whether to plot a histogram of the PAI scores.
#' @param dots Additional arguments, not currently used.
#'   a label for the results.
#' @return A list of data frames with an element of summary statistics
#'   for the PAI (PAI) and observed results for patients assigned to
#'   their optimal and nonoptimal treatment (ObservedOutcome).
#' @export
#' @import boot pander
#' @examples
#' \dontrun{
#' # builtin dataset with 32 cases
#' dat <- mtcars
#' dat$cyl <- factor(dat$cyl)
#'
#' m <- mi_lm_pai_cvboot(~ cyl + am * (mpg + hp + drat + wt),
#'   "disp", "am", list(dat, dat),
#'   nboot = 50, holdouts = "10", cores = 2)
#'
#' summary(m)
#' }
summary.mipai <- function(object, plot = TRUE, ...) {
  stopifnot(inherits(object, "mipai"))
  pais <- apply(object$nobootresults[, "Yhat.PAI", , ], 1, function(x) {
    mean(abs(x), na.rm = TRUE)
  })

  output <- data.frame(
    Mean = mean(pais, na.rm=TRUE),
    Median = median(pais, na.rm=TRUE),
    Min = min(pais, na.rm=TRUE),
    Max = max(pais, na.rm=TRUE),
    SD = sd(pais, na.rm=TRUE)
  )

  if (inherits(object, "mibootpai")) {
    pai.means <- matrix(NA_real_, nrow = dim(object$bootresults)[3],
      ncol = dim(object$bootresults)[4])

    for (i in 1:dim(object$bootresults)[3]) {
      for (j in 1:dim(object$bootresults)[4]) {
        x <- object$bootresults[, , i, j]
        pai.means[i, j] <- mean(abs(x[, "Yhat.PAI"]), na.rm=TRUE)
      }
    }

    bres <- boot(data.frame(X = 1:dim(object$bootresults)[1]),
      function(d, i) 1,
      R = prod(dim(object$bootresults)[3:4]))

    bres$data <- rep(NA_real_, dim(object$bootresults)[1])

    bres$t0 <- mean(pais, na.rm = TRUE)
    bres$t <- cbind(as.vector(pai.means))

    ci.bca <- tryCatch(boot.ci(bres, index = 1, type = "bca"),
                  error = function(x) x)

    ci.perc <- tryCatch(boot.ci(bres, index = 1, type = "perc"),
                   error = function(x) x)

    if (inherits(ci.bca, "error")) {
      LL.bca <- NA
      UL.bca <- NA
      message(ci.bca)
    } else {
      LL.bca <- ci.bca$bca[1, 4]
      UL.bca <- ci.bca$bca[1, 5]
    }

    if (inherits(ci.perc, "error")) {
      LL.perc <- NA
      UL.perc <- NA
      message(ci.perc)
    } else {
      LL.perc <- ci.perc$percent[1, 4]
      UL.perc <- ci.perc$percent[1, 5]
    }

    output <- cbind(output,
      LLPerc = LL.perc, ULPerc = UL.perc,
      LLBCa = LL.bca, ULBCa = UL.bca)
  }

  obs.output <- observed_outcome(object)

  if (plot) {
    hist(pais, main = "Histogram of PAI Scores", xlab = "PAI Score")
  }

  pander(output, digits = 3, split.tables = 200,
    caption = "Summary Statistics for the Average PAI")

  pander(obs.output, digits = 3, split.tables = 200,
    caption = "Observed Outcome for Patients Receiving Optimal and NonOptimal Treatment")

  return(invisible(list(PAI = output, ObservedOutcome = obs.output)))
}



#' Calculate the Average Outcome for Those in Optimal and Non-optimal Treatment
#'
#' This function claculates the average outcome for cases assigned to their optimal and nonoptimal treatment.  If bootstrapping was performed, confidence intervals based on the bootstrap.  This is done by calculating the optimal treatment based on predicted PAI scores from each bootstrap replication, and using these to select cases and calculate the average (observed) outcome.
#'
#' @param object an mipai or mibootpai class object returned from
#'   \code{mi_lm_pai_cvboot}.
#' @return A data frame with summary statistics.
#' @export
#' @import boot
#' @examples
#' \dontrun{
#' # builtin dataset with 32 cases
#' dat <- mtcars
#' dat$cyl <- factor(dat$cyl)
#'
#' m <- mi_lm_pai_cvboot(~ cyl + am * (mpg + hp + drat + wt),
#'   "disp", "am", list(dat, dat),
#'   nboot = 50, holdouts = "10", cores = 2)
#'
#' observed_outcome(m)
#' }
observed_outcome <- function(object) {
  stopifnot(inherits(object, "mipai"))

  # average by optimal/nonoptimal across all results
  lucky.index <- which(as.vector(object$nobootresults[, "Lucky", , ]) == 1)
  unlucky.index <- which(as.vector(object$nobootresults[, "Lucky", , ]) == 0)

  noboot.mean <- c(
    Optimal = mean(
      as.vector(object$nobootresults[, "Y", , ])[lucky.index], na.rm = TRUE),
    NonOptimal = mean(
      as.vector(object$nobootres[, "Y", , ])[unlucky.index], na.rm = TRUE))

  noboot.sd <- c(
    Optimal = sd(
      as.vector(object$nobootresults[, "Y", , ])[lucky.index], na.rm = TRUE),
    NonOptimal = sd(
      as.vector(object$nobootres[, "Y", , ])[unlucky.index], na.rm = TRUE))

  noboot.median <- c(
    Optimal = median(
      as.vector(object$nobootresults[, "Y", , ])[lucky.index], na.rm = TRUE),
    NonOptimal = median(
      as.vector(object$nobootres[, "Y", , ])[unlucky.index], na.rm = TRUE))

  output <- data.frame(
    Mean = noboot.mean,
    Median = noboot.median,
    SD = noboot.sd)

  if (inherits(object, "mibootpai")) {
    lucky.means <- matrix(NA_real_, nrow = dim(object$bootresults)[3],
      ncol = dim(object$bootresults)[4])

    for (i in 1:dim(object$bootresults)[3]) {
      for (j in 1:dim(object$bootresults)[4]) {
        x <- object$bootresults[, , i, j]
        lucky.means[i, j] <- mean(x[x[, "Lucky"] == 1, "Y"], na.rm=TRUE)
      }
    }

    unlucky.means <- matrix(NA_real_, nrow = dim(object$bootresults)[3],
      ncol = dim(object$bootresults)[4])

    for (i in 1:dim(object$bootresults)[3]) {
      for (j in 1:dim(object$bootresults)[4]) {
        x <- object$bootresults[, , i, j]
        unlucky.means[i, j] <- mean(x[x[, "Lucky"] == 0, "Y"], na.rm=TRUE)
      }
    }

    bres <- boot(1:2, function(d, i) 1,
      R = prod(dim(object$bootresults)[3:4]))

    bres$data <- rep(NA_real_, 2)

    bres$t0 <- noboot.mean
    bres$t <- cbind(MeanOptimal = as.vector(lucky.means),
                    MeanNonOptimal = as.vector(unlucky.means))

    ci.bca <- list(
      Optimal = boot.ci(bres, index = 1, type = "bca"),
      NonOptimal = boot.ci(bres, index = 2, type = "bca"))

    ci.perc <- list(
      Optimal = boot.ci(bres, index = 1, type = "perc"),
      NonOptimal = boot.ci(bres, index = 2, type = "perc"))

    output <- cbind(output,
      LLPerc = c(Optimal = ci.perc[[1]]$percent[1, 4],
                 NonOptimal = ci.perc[[2]]$percent[1, 4]),
      ULPerc = c(Optimal = ci.perc[[1]]$percent[1, 5],
                 NonOptimal = ci.perc[[2]]$percent[1, 5]),
      LLBCa = c(Optimal = ci.bca[[1]]$bca[1, 4],
                NonOptimal = ci.bca[[2]]$bca[1, 4]),
      ULBCa = c(Optimal = ci.bca[[1]]$bca[1, 5],
                NonOptimal = ci.bca[[2]]$bca[1, 5])
    )
  }
  return(output)
}

