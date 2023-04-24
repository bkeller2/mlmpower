# mlmpower
An R Package for simulating and completing  power analysis for multilevel models

[![R-CMD-check](https://github.com/bkeller2/mlmpower/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/bkeller2/mlmpower/actions/workflows/check-standard.yaml)

## Installing Package

To install the package via CRAN run the following:
```r
install.packages('mlmpower')
```

To install the dev build via GitHub use the following command:

```r
remotes::install_github('bkeller2/mlmpower', build_vignettes = TRUE)
```

## Using the Package

A guide to using the package is provided through the vignette.

If you set `build_vignettes = TRUE` then the following command should open up said vignette:
```r
vignette('mlmpower')
```

Additional details can be obtained throughout the documentation, starting with:
```r
?mlmpower::mlmpower
```
or

```r
help(package = 'mlmpower')
```

## Bug Reporting and Feature Request
Please use the [Issues](https://github.com/bkeller2/mlmpower/issues) tab.

## To Do List
- [x] Complete and Test Implementation
- [x] Add Documentation
- [x] Add Vignettes
- [ ] Release 1.0 and accepted to CRAN
