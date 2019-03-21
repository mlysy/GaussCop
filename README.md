# GaussCop: Utilities for the Gaussian Copula Distribution

*Martin Lysy, Jonathan Ramkissoon*

*June 26, 2018*

---

### Description

Provides functions for density evaluation and both joint and conditional random sampling from a Gaussian Copula distribution with arbitrary margins.  This is achieved by storing univariate marginal distributions as *extended density* (`xDensity`) objects which provide `d/p/q/r` methods, and which can be estimated from data using a variety of parametric, nonparametric, and semiparametric approaches.

### Installation

Install the R package [devtools](https://CRAN.R-project.org/package=devtools) and run
```r
devtools::install_github("mlysy/GaussCop")
```

### Usage

An overview of the `xDensity` class and other package utilities is provided in this [vignette](http://htmlpreview.github.com/?https://github.com/mlysy/GaussCop/master/inst/doc/xDensity-quicktut.html).
