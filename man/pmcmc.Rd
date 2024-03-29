% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmcmc.R
\name{pmcmc}
\alias{pmcmc}
\title{Run a pmcmc sampler}
\usage{
pmcmc(pars, filter, initial = NULL, control = NULL)
}
\arguments{
\item{pars}{A \code{\link{pmcmc_parameters}} object containing
information about parameters (ranges, priors, proposal kernel,
translation functions for use with the particle filter).}

\item{filter}{A \code{\link{particle_filter}} object}

\item{initial}{Optional initial starting point. If given, it must
be compatible with the parameters given in \code{pars}, and must be
valid against your prior. You can use this to override the
initial conditions saved in your \code{pars} object. You can provide
either a vector of initial conditions, or a matrix with
\code{n_chains} columns to use a different starting point for each
chain.}

\item{control}{A \link{pmcmc_control} object which will
control how the MCMC runs, including the number of steps etc.}
}
\value{
A \code{mcstate_pmcmc} object containing \code{pars}
(sampled parameters) and \code{probabilities} (log prior, log
likelihood and log posterior values for these
probabilities). Two additional fields may be present:
\code{state} (if \code{return_state} was \code{TRUE}),
containing the final state of a randomly selected particle at
the end of the simulation, for each step (will be a matrix with
as many rows as your state has variables, and as \code{n_steps + 1} columns corresponding to each step). \code{trajectories} will
include a 3d array of particle trajectories through the
simulation (if \code{return_trajectories} was \code{TRUE}).
}
\description{
Run a pmcmc sampler
}
\details{
This is a basic Metropolis-Hastings MCMC sampler.  The
\code{filter} is run with a set of parameters to evaluate the
likelihood. A new set of parameters is proposed, and these
likelihoods are compared, jumping with probability equal to their
ratio. This is repeated for \code{n_steps} proposals.

While this function is called \code{pmcmc} and requires a particle
filter object, there's nothing special about it for particle
filtering. However, we may need to add things in the future that
make assumptions about the particle filter, so we have named it
with a "p".
}
