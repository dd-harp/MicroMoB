% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-basicM.R
\name{dMYZdt.basicM}
\alias{dMYZdt.basicM}
\title{Derivatives for adult mosquitoes}
\usage{
\method{dMYZdt}{basicM}(t, y, pars, s)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{a \link{list}}

\item{s}{the species index}
}
\value{
a \link{numeric} vector
}
\description{
Implements \link{dMYZdt} for the basicM model.
}
