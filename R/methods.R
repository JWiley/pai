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

#' Summary Statistics in a Named Vector
#'
#' Simple utility wrapper to make it easy to get a variety of summary
#' statistics in a named vector, in whatever order desired.
#'
#' @param x A data vector
#' @param types The types of summaries to be computed.
#' @param na.rm Logical passed on to stat functions.
#' @return a named vector
#' @export
#' @importFrom psych skew kurtosi
#' @keywords utils
#' @examples
#'
#' statsummary(1:10, c("Mean", "Min", "Max", "IQR", "Kurtosis"))
statsummary <- function(x, types, na.rm = TRUE) {
  res <- sapply(types, function(f) {
    fun <- switch(f,
      Mean = mean,
      Median = median,
      SD = sd,
      IQR = function(x, na.rm) as.vector(abs(diff(quantile(x, probs = c(.25, .75), na.rm = na.rm)))),
      Min = min,
      Max = max,
      Skew = skew,
      Kurtosis = kurtosi,
      LL90 = function(x, na.rm) as.vector(quantile(x, probs = .1, na.rm = na.rm)),
      LL95 = function(x, na.rm) as.vector(quantile(x, probs = .05, na.rm = na.rm)),              UL90 = function(x, na.rm) as.vector(quantile(x, probs = .9, na.rm = na.rm)),
      UL95 = function(x, na.rm) as.vector(quantile(x, probs = .95, na.rm = na.rm))
    )
    fun(x, na.rm = na.rm)
  })
  names(res) <- types
  return(res)
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

#' Calculate the Residuals for Predictions
#'
#' @param object an mipai or mibootpai class object returned from
#'   \code{mi_lm_pai_cvboot}
#' @param type A character string indicating the type of residuals to be calculated.
#'   Currents mean or median of the absolute or root squared values.  Only applies
#'   when summarize is \code{FALSE}.
#' @param summarize Logical whether or not to summarize values for individuals
#'   across imputations or leave as a matrix of residuals.  Defaults to \code{TRUE}.
#' @param boot Not currently used.
#' @param dots Additional arguments to pass down.  Not currently used.
#' @return A vector of residuals for each person or matrix across the
#'   multiple imputations.
#' @export
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
#' residuals(m, "RootMeanSq")
#' }
residuals.mipai <- function(object, type = c("MeanAbs", "MedianAbs", "RootMeanSq", "RootMedianSq"), summarize = TRUE, boot = TRUE, ...) {
  stopifnot(inherits(object, "mipai"))
  type <- match.arg(type)

  if (summarize) {
    apply(
      object$nobootresults[, "Y", , ] - object$nobootresults[, "Yhat.reality", , ],
      1, switch(type,
        MeanAbs = function(x) mean(abs(x), na.rm = TRUE),
        MedianAbs = function(x) median(abs(x), na.rm = TRUE),
        RootMeanSq = function(x) sqrt(mean(x^2, na.rm = TRUE)),
        RootMedianSq = function(x) sqrt(median(x^2, na.rm = TRUE)))
      )
  } else {
    ## if (inherits(object, "mibootpai") && boot) {
    ##   object$bootresults[, "Y", , ] - object$bootresults[, "Yhat.reality", , ]
    ## } else {
      object$nobootresults[, "Y", , ] - object$nobootresults[, "Yhat.reality", , ]
    ## }
  }
}

#' Calculate Summary Statistics for PAI Model
#'
#' @param object an mipai or mibootpai class object returned from
#'   \code{mi_lm_pai_cvboot}.
#' @param plot logical whether to plot a histogram of the PAI scores.
#' @param dots Additional arguments, not currently used.
#'   a label for the results.
#' @param boot logical whether to use bootstrapped results. Defaults to
#'   \code{TRUE} if the object includes bootstrapped results.
#' @return A list of data frames with an element of summary statistics
#'   for the PAI (PAI) and observed results for patients assigned to
#'   their optimal and nonoptimal treatment (ObservedOutcome) and
#'   residuals (Residuals), and for the histogram (Graph).
#' @export
#' @import boot pander ggplot2
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
summary.mipai <- function(object, plot = TRUE, boot = TRUE, ...) {
  stopifnot(inherits(object, "mipai"))
  pais <- apply(object$nobootresults[, "Yhat.PAI", , ], 1, function(x) {
    mean(x, na.rm = TRUE)
  })

  output <- data.frame(
    Mean = mean(abs(pais), na.rm=TRUE),
    Median = median(abs(pais), na.rm=TRUE),
    Min = min(abs(pais), na.rm=TRUE),
    Max = max(abs(pais), na.rm=TRUE),
    SD = sd(abs(pais), na.rm=TRUE)
  )

  resoutput <- data.frame(AbsError =
    c(Mean = mean(residuals(object, "MeanAbs"), na.rm = TRUE),
      Median = median(residuals(object, "MedianAbs"), na.rm = TRUE)))

  if (inherits(object, "mibootpai") && boot) {
    pai.means <- matrix(NA_real_, nrow = dim(object$bootresults)[3],
      ncol = dim(object$bootresults)[4])

    for (i in 1:dim(object$bootresults)[3]) {
      for (j in 1:dim(object$bootresults)[4]) {
        x <- object$bootresults[, , i, j]
        pai.means[i, j] <- mean(x[, "Yhat.PAI"], na.rm=TRUE)
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

    if (FALSE) {
    e <- residuals(object, summarize = FALSE, boot = TRUE)

    e.medians <- e.means <- matrix(NA_real_, nrow = dim(object$bootresults)[3],
      ncol = dim(object$bootresults)[4])

    for (i in 1:dim(object$bootresults)[3]) {
      for (j in 1:dim(object$bootresults)[4]) {
        x <- e[, i, j]
        e.means[i, j] <- mean(abs(x), na.rm=TRUE)
        e.medians[i, j] <- median(abs(x), na.rm=TRUE)
      }
    }

    ebres <- boot(data.frame(X = 1:dim(object$bootresults)[1]),
      function(d, i) 1,
      R = prod(dim(object$bootresults)[3:4]))

    ebres$data <- rep(NA_real_, dim(object$bootresults)[1])

    ebres$t0 <- resoutput$AbsError
    ebres$t <- cbind(as.vector(e.means), as.vector(e.medians))

    e.ci.res <- do.call(rbind, lapply(1:2, function(i) {
      e.ci.bca <- tryCatch(boot.ci(ebres, index = i, type = "bca"),
                    error = function(x) x)
      e.ci.perc <- tryCatch(boot.ci(ebres, index = i, type = "perc"),
                     error = function(x) x)

      if (inherits(e.ci.bca, "error")) {
        e.LL.bca <- NA
        e.UL.bca <- NA
        message(e.ci.bca)
      } else {
        e.LL.bca <- e.ci.bca$bca[1, 4]
        e.UL.bca <- e.ci.bca$bca[1, 5]
      }

      if (inherits(e.ci.perc, "error")) {
        e.LL.perc <- NA
        e.UL.perc <- NA
        message(e.ci.perc)
      } else {
        e.LL.perc <- e.ci.perc$percent[1, 4]
        e.UL.perc <- e.ci.perc$percent[1, 5]
      }

      data.frame(LLPerc = e.LL.perc, ULPerc = e.UL.perc,
                 LLBCa = e.LL.bca, ULBCa = e.UL.bca)
    }))
    rownames(e.ci.res) <- c("Mean", "Median")

    resoutput <- cbind(resoutput, e.ci.res)
  }
  }

  obs.output <- observed_outcome(object)

  g <- ggplot(data.frame(PAI = pais), aes(PAI)) +
    geom_histogram() +
    theme_classic() +
    xlab("PAI Score") +
    ggtitle("Histogram of PAI Scores")

  if (plot) print(g)

  pander(output, digits = 3, split.tables = 200,
    caption = "Summary Statistics for the Average PAI")

  pander(obs.output, digits = 3, split.tables = 200,
    caption = "Observed Outcome for Patients Receiving Optimal and NonOptimal Treatment")

  pander(resoutput, digits = 3, split.tables = 200,
    caption = "Summary Statistics for the Residuals (Predicted Outcome - True Outcome)")


  return(invisible(list(PAI = output, ObservedOutcome = obs.output,
                        Residuals = resoutput, Graph = g)))
}


#' Calculate Individual PAIs
#'
#' @param object an mipai or mibootpai class object returned from
#'   \code{mi_lm_pai_cvboot}.
#' @param plot logical whether to plot a dotplot.
#' @param confint logical whether to calculate confidence intervals.
#'   Only applies to bootstrapped model results.  Defaults to
#'   \code{TRUE}.
#' @param cl An existing cluster to use (optional)
#' @param cores The number of cores to use when creating a cluster, if an existing
#'   cluster is not passed (optional).  If left blank and no cluster passed,
#'   defaults to the number of cores available, but is only used when calculating
#'   bootstrapped confidence intervals.
#' @param dots Additional arguments, not currently used.
#'   a label for the results.
#' @return An invisible list of two elements, the \dQuote{PAI} scores and
#'   the \dQuote{Graph}, which is returned even if not plotted.
#' @export
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
#' # use two cores to speed up when obtaining bootstrapped CIs
#' predict(m, cores = 2L)
#' }
predict.mipai <- function(object, plot = TRUE, confint = TRUE, cores, cl, ...) {
  stopifnot(inherits(object, "mipai"))

  if (inherits(object, "mibootpai")) {
    if (missing(cl)) {
      if (missing(cores)) {
        cores <- detectCores()
      }
      cl <- makeCluster(cores)
      on.exit(stopCluster(cl))
    }

    env <- environment()
    clusterEvalQ(cl, library(pai))
  }


  pais <- apply(object$nobootresults[, "Yhat.PAI", , ], 1, function(x) {
    mean(x, na.rm = TRUE)
  })

  finalout <- data.frame(PAI = pais)

  if (inherits(object, "mibootpai") && confint) {
    spais <- as.data.frame(t(apply(object$bootresults[, "Yhat.PAI", , ], 1, function(x) {
      statsummary(x, c("Mean", "Median", "SD", "IQR", "LL95", "UL95"))
    })))

    bres <- boot(data.frame(X = 1:dim(object$bootresults)[1]),
      function(d, i) 1,
      R = prod(dim(object$bootresults)[3:4]))

    bres$data <- rep(NA_real_, dim(object$bootresults)[1])

    bres$t0 <- pais
    bres$t <- apply(object$bootresults[, "Yhat.PAI", , ], 1, as.vector)

    clusterExport(cl, "bres", env)

    output <- do.call(rbind, parLapplyLB(cl, 1:dim(object$bootresults)[1], function(i) {
      ci.bca <- tryCatch(boot.ci(bres, index = i, type = "bca"),
                  error = function(x) x)
      ci.perc <- tryCatch(boot.ci(bres, index = i, type = "perc"),
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

      data.frame(LLPerc = LL.perc, ULPerc = UL.perc,
        LLBCa = LL.bca, ULBCa = UL.bca)
    }))

    finalout <- cbind(PAI = pais, output)
  }

  pdat <- finalout[order(finalout$PAI), , drop = FALSE]
  pdat$Index <- 1:nrow(pdat)

  if (inherits(object, "mibootpai") && confint) {
    g <- ggplot(pdat, aes(Index, PAI)) +
      geom_linerange(aes(ymin = LLPerc, ymax = ULPerc)) +
      geom_point() +
      theme_classic() +
      xlab("Index") + ylab("Predicted PAI Score") +
      ggtitle("PAI Scores with 95% Bootstrap Percentile CIs")
  } else {
    g <- ggplot(pdat, aes(Index, PAI)) +
      geom_point() +
      theme_classic() +
      xlab("Index") + ylab("Predicted PAI Score") +
      ggtitle("PAI Scores")
  }

  if (plot) {
    print(g)
  }

  return(invisible(list(PAI = finalout, Graph = g)))
}

