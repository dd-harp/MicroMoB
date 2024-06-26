% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-basicM.R
\name{MBionomics.basicM}
\alias{MBionomics.basicM}
\title{Reset bloodfeeding and mortality rates to baseline}
\usage{
\method{MBionomics}{basicM}(t, y, pars, s)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{a \link{list}}

\item{s}{the species index}
}
\value{
the model as a \link{list}
}
\description{
Implements \link{MBionomics} for the basicM model
}
