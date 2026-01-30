---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# csbSTATS

<!-- badges: start -->
<!-- badges: end -->

**csbSTATS** is is an open-source project that provides an user-friendly,
interactive user interface for exploratory data analysis, visualization,
and statistical testing.

## Installation

**csbSTATS** can be installed by first installing R and typing the following
commands in R console:

``` r
#to allow installation of GitHub packages
install.packages("devtools")

#csbSTATS can then be installed
devtools::install_github("BonilhaCaio/csbSTATS")
```

## Usage

**csbSTATS** Shiny-based user interface with all its tools is generated through
the package's single function:

``` r
csbSTATS::runCsbSTATS()
```

## Credits and citation

A publication where the package dependencies and imports will be appropriately
credited is on the way. Meanwhile, you can cite **csbSTATS** as follows:

``` r
citation("csbSTATS")
```
