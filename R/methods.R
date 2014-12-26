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



#' Calculate Summary Statistics for PAI Scores
#'
#' @param object an mipai or mibootpai class object returned from
#'   \code{mi_lm_pai_cvboot}.
#' @param plot logical whether to plot a histogram of the PAI scores.
#' @param dots Additional arguments, not currently used.
#'   a label for the results.
#' @return A list of data frames with an element of summary statistics
#'   for the PAI (PAI) and observed results for patients assigned to
#'   their optimal and nonoptimal treatment (ObservedOutcome) and
#'   for the histogram (Graph).
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

  return(invisible(list(PAI = output, ObservedOutcome = obs.output, Graph = g)))
}


#' Calculate Individual PAIs
#'
#' @param object an mipai or mibootpai class object returned from
#'   \code{mi_lm_pai_cvboot}.
#' @param plot logical whether to plot a dotplot.
#' @param cores An integer, the number of cores to be used for a local cluster
#'   when calculating CIs on bootstrapped models.
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
#' predict(m, cores = 2L)
#' }
predict.mipai <- function(object, plot = TRUE, cores = 1L, ...) {
  stopifnot(inherits(object, "mipai"))

  if (inherits(object, "mibootpai")) {
    env <- environment()
    cl <- makeCluster(cores)
    clusterEvalQ(cl, library(pai))

    on.exit(stopCluster(cl))
  }


  pais <- apply(object$nobootresults[, "Yhat.PAI", , ], 1, function(x) {
    mean(x, na.rm = TRUE)
  })

  finalout <- data.frame(PAI = pais)

  if (inherits(object, "mibootpai")) {
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

  if (inherits(object, "mibootpai")) {
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

