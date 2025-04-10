% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repsetup.R
\name{repsetup}
\alias{repsetup}
\title{Setup for Analysis with Replicate Weights}
\usage{
repsetup(
  repwt,
  wt,
  df,
  method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA", "TALIS"),
  group = NULL,
  exclude = NULL
)
}
\arguments{
\item{repwt}{a string indicating the common names for the replicate weights
columns (within \code{df}), or a data frame with the replicate weights.}

\item{wt}{a string specifying the name of the column (within \code{df}) with the
total weights.}

\item{df}{a data frame.}

\item{method}{a string indicating the name of the large-scale assessment
to determine the replication method to use. Available options are:
\code{"TIMSS"}, \code{"PIRLS"}, \code{"ICILS"}, \code{"ICCS"},
\code{"PISA"}, and \code{"TALIS"}. Note that \code{"TIMSS"} and \code{"PIRLS"}
refer to the method used from 2016 onwards.
Their method has not yet been implemented for previous cycles.}

\item{group}{a string specifying the variable name (within \code{df}) to be
used for grouping. Categories in \code{group} are treated as independent, e.g.,
countries.}

\item{exclude}{a vector indicating which groups
(in the same format as \code{group})
should be excluded from the pooled and composite estimates.}
}
\value{
a list to be used in other functions.
}
\description{
Creates a list with common arguments used for analysis with replicate
weights.
}
\examples{
# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

### No groups ----
stp1 <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "ICILS")
stp1

### Groups ----
stp2 <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "ICILS",
                 group = "GROUP", exclude = "GR2")
stp2


### repmean ----

repmean(x = "Math1",setup = stp1)

repmean(x = "Math1",setup = stp2)
}
