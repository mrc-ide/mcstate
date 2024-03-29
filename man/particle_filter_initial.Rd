% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/particle_filter_initial.R
\name{particle_filter_initial}
\alias{particle_filter_initial}
\title{Create restart initial state}
\usage{
particle_filter_initial(state)
}
\arguments{
\item{state}{A matrix of state (rows are different states, columns
are different realisations).  This is the form of a slice pulled
from a restart.}
}
\value{
A function with arguments \code{info}, \code{n_particles} and \code{pars}
that will sample, with replacement, a matrix of state suitable
as a starting point for a particle filter.  The \code{info} and
\code{pars} arguments are ignored.
}
\description{
Create a suitable initial condition function from a set of restart
state.  This takes care of a few bookkeping and serialisation
details and returns a function appropriate to pass to
\link{particle_filter} as \code{initial}.
}
