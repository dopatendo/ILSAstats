% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repprop.R
\name{repprop}
\alias{repprop}
\title{Proportions with Replicate Weights}
\usage{
repprop(
  x,
  categories = NULL,
  setup = NULL,
  repwt,
  wt,
  df,
  method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA", "TALIS"),
  group = NULL,
  exclude = NULL
)
}
\arguments{
\item{x}{a string vector specifying variable names (within \code{df}) for
analysis.}

\item{categories}{a vector indicating all possible response categories.
If \code{NULL}, categories will be derived from the data.}

\item{setup}{an optional list produced by \code{\link{repsetup}}.}

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
a list.
}
\description{
Estimates proportions using replicate weights
for a variable or a group of plausible values variables and for one or more
populations.
For a detailed explanation on how the standard errors are estimated
see \code{\link{repse}}.
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

# One variable - weights within df
repprop(x = c("item01"),
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
        method = "ICILS")

# One variable - weights weights as a separate data frame
repprop(x = c("item01"),
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

# Multiple variables - PVs are assumed
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

### Groups ----

# One variable - weights within df
repprop(x = c("item01"),
        group = "GROUP",
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
        method = "ICILS")

# One variable - weights weights as a separate data frame
repprop(x = c("item01"),
        group = "GROUP",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

# Multiple variables - PVs are assumed
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        group = "GROUP",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

# Multiple variables - excluding one group
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        group = "GROUP",
        exclude = "GR2",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

}
