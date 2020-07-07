## mcstate <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![Build Status](https://travis-ci.com/mrc-ide/mcstate.svg?branch=master)](https://travis-ci.com/mrc-ide/mcstate)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/mcstate/badge)](https://www.codefactor.io/repository/github/mrc-ide/mcstate)
[![codecov.io](https://codecov.io/github/mrc-ide/mcstate/coverage.svg?branch=master)](https://codecov.io/github/mrc-ide/mcstate?branch=master)
<!-- badges: end -->

## Installation

Install from drat with

```
# install.packages("drat") # -- if you don't have drat installed
drat:::add("mrc-ide")
install.packages("mcstate")
```

You will need a compiler to install dependencies for the package, and to build any models.  Use `pkgbuild::check_build_tools()` to see if your system is ok to use.

The development version of the package can be installed directly from GitHub if you prefer with:

```r
remotes::install_github("mrc-ide/mcstate", upgrade = FALSE)
```

## License

MIT © Imperial College of Science, Technology and Medicine
