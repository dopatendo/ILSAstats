% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reprho.R
\name{reprho}
\alias{reprho}
\title{Correlations with Replicate Weights}
\usage{
reprho(
  x = NULL,
  pv = NULL,
  pv2 = NULL,
  relatedpvs = TRUE,
  setup = NULL,
  repwt,
  wt,
  df,
  rho = c("pearson", "spearman", "polychoric"),
  method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA", "TALIS"),
  group = NULL,
  exclude = NULL
)
}
\arguments{
\item{x}{a string vector specifying variable names (within \code{df}) for analysis.
If \code{pv} is \code{NULL}, this function estimates correlations between all variables in the vector.
If \code{pv2} is NOT \code{NULL}, then \code{x} should be set to \code{NULL}.}

\item{pv}{a string vector indicating the variable names for all plausible values
of a construct. If not \code{NULL}, this function estimates correlations only between
\code{x} and the plausible values construct.}

\item{pv2}{a string vector indicating the variable names for all plausible values
of a second construct (distinct from \code{pv}).}

\item{relatedpvs}{a logical value indicating if \code{pv} and \code{pv2} are drawn
from the same model, and have the same number of plausible values.
If \code{TRUE} (default), a total of \eqn{n} estimations will be done,
where \eqn{n} is the number of plausible values of each.
If \code{FALSE}, a total of \eqn{n_1 \times n_2}
estimations will be done, where \eqn{n_1} is the number of plausible values in \code{pv}
and \eqn{n_2} is the number of plausible values in \code{pv2}.}

\item{setup}{an optional list produced by \code{\link{repsetup}}.}

\item{repwt}{a string indicating the common names for the replicate weights
columns within \code{df}, or a data frame with the replicate weights.}

\item{wt}{a string specifying the name of the column in \code{df} that contains the total weights.}

\item{df}{a data frame.}

\item{rho}{a string indicating the correlation coefficient to be computed:
\code{"pearson"}, \code{"polychoric"}, or \code{"spearman"} (lower or uppercase).}

\item{method}{a string indicating the name of the large-scale assessment to
determine the replication method to use. Available options are:
\code{"TIMSS"}, \code{"PIRLS"}, \code{"ICILS"}, \code{"ICCS"}, and \code{"PISA"}.}

\item{group}{a string specifying the variable name (within \code{df}) to be used for grouping.}

\item{exclude}{a vector indicating which groups (in the same format as \code{group})
should be excluded from the estimation of pooled and composite estimates.}
}
\value{
a data frame.
}
\description{
Estimates correlation coefficients using replicate weights.
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

# Non PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = NULL,
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")

# X var and PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = paste0("Reading",1:5),
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")

# PVs and PVs (related)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")

# PVs and PVs (UNrelated)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       relatedpvs = FALSE,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")


### Groups ----

# Non PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = NULL,
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")

# X var and PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = paste0("Reading",1:5),
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")

# PVs and PVs (related)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")

# PVs and PVs (UNrelated)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       relatedpvs = FALSE,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")
}
