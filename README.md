# GaussCop: Utilities for the Gaussian Copula Distribution

*Martin Lysy, Jonathan Ramkissoon*

---

### Description

Provides functions for density evaluation and both joint and conditional random sampling from a Gaussian Copula distribution with arbitrary margins.  This is achieved by storing univariate marginal distributions as *extended density* (`xDensity`) objects which provide `d/p/q/r` methods, and which can be estimated from data using a variety of parametric, nonparametric, and semiparametric approaches.

### Installation

Install the R package [**devtools**](https://CRAN.R-project.org/package=devtools) and run
```r
devtools::install_github("mlysy/GaussCop")
```

### Usage

An overview of the `GaussCop::xDensity` class and other package utilities is provided in the package [vignette](http://htmlpreview.github.io/?https://github.com/mlysy/GaussCop/master/doc/xDensity-quicktut.html).
