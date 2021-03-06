% Generated by roxygen2 (4.0.1): do not edit by hand
\name{observed_outcome}
\alias{observed_outcome}
\title{Calculate the Average Outcome for Those in Optimal and Non-optimal Treatment}
\usage{
observed_outcome(object)
}
\arguments{
\item{object}{an mipai or mibootpai class object returned from
\code{mi_lm_pai_cvboot}.}
}
\value{
A data frame with summary statistics.
}
\description{
This function claculates the average outcome for cases assigned to their optimal and nonoptimal treatment.  If bootstrapping was performed, confidence intervals based on the bootstrap.  This is done by calculating the optimal treatment based on predicted PAI scores from each bootstrap replication, and using these to select cases and calculate the average (observed) outcome.
}
\examples{
\dontrun{
# builtin dataset with 32 cases
dat <- mtcars
dat$cyl <- factor(dat$cyl)

m <- mi_lm_pai_cvboot(~ cyl + am * (mpg + hp + drat + wt),
  "disp", "am", list(dat, dat),
  nboot = 50, holdouts = "10", cores = 2)

observed_outcome(m)
}
}

