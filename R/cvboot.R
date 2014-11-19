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
#' @param seed A seed to make the bootstrap reproducible
#' @return A three dimensional array with cases on the first dimension,
#'   multiple estimates calculated on the second dimension, and bootstrap resamples
#'   on the third dimension.
#' @import Rcpp RcppEigen
#' @export
#' @examples
#'
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
#' m <- lm_pai_cvboot(formula, DV, TreatVar, dat, boot = TRUE, k = k, holdouts = holdouts)
#'
#' # see result from first bootstrap resample
#' m[, , 1]
lm_pai_cvboot <- function(formula, DV, TreatVar, dat, boot = TRUE, k, holdouts, yhigherisbetter = TRUE, seed = 1234) {

  stopifnot(all(unique(dat[, TreatVar]) %in% c(0, 1)))

  res <- lapply(holdouts, function(n) {
    nholdout <- length(n)
    ntrain <- nrow(dat) - nholdout

    if (boot) {
      set.seed(seed)
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
      out <- cbind(out,
                   Yhat.PAI = out[, "Yhat.0"] - out[, "Yhat.1"],
                   Residual = out[, "Y"] - out[, "Yhat.reality"])

      if (yhigherisbetter) {
        out <- cbind(out, Lucky = ifelse(out[, "Yhat.0"] < out[, "Yhat.1"], 1, 0))
      } else {
        out <- cbind(out, Lucky = ifelse(out[, "Yhat.0"] < out[, "Yhat.1"], 0, 1))
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
#' @param imputedData An imputed dataset (list of datasets) to use
#' @param boot A logical value whether to bootstrap or not
#' @param k The number of bootstrap resamples to take
#' @param holdouts The holdouts to be used.  If a vector from 1 to the number of
#'   cases available in the dataset, it is equivalent to leave-one-out.
#'   Must be constructed manually.
#' @param yhigherisbetter A logical value whether higher scores indicate better
#'   treatment response on the dependent variable
#' @param seed A seed to make the bootstrap reproducible
#' @param cl An existing cluster to use (optional)
#' @param cores The number of cores to use when creating a cluster, if an existing
#'   cluster is not passed (optional).  If left blank and no cluster passed,
#'   defaults to the number of cores available.
#' @return A three dimensional array with cases on the first dimension,
#'   multiple estimates calculated on the second dimension, and bootstrap resamples
#'   on the third dimension.
#' @export
#' @import doParallel doRNG mice
#' @examples
#' # make me!!!
mi_lm_pai_cvboot <- function(formula, DV, TreatVar, imputedData, boot = TRUE, k = 500, holdouts, yhigherisbetter = TRUE, seed = 1234, cl, cores) {
  if (missing(cl)) {
    if (missing(cores)) {
      cores <- detectCores()
    }
    cl <- makeCluster(cores)
  }

  on.exit(stopCluster(cl))

  clusterEvalQ(cl, {
    library(pai)
  })

  Nimp <- length(imputedData)

  clusterExport(cl, c("indexer", "lm_pai_cvboot", "holdouts"))

  env <- environment()
  clusterExport(cl, c("imputedData"), envir = env)

  set.seed(seed)

  parLapplyLB(cl, 1:Nimp, function(i) {
    lm_pai_cvboot(formula = formula,
      DV = DV, TreatVar = TreatVar, dat = imputedData[[i]],
      boot = boot, k = k, holdouts = holdouts,
      yhigherisbetter = yhigherisbetter, seed = seed)
  })
}
