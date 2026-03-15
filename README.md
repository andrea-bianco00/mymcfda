# mymcfda

`mymcfda` is an independent R package for simulation and estimation in the context of functional snippet data.

The package was developed after studying the methodology introduced in:

Lin, Z. and Wang, J.-L. (2022). *Mean and Covariance Estimation for Functional Snippets*. Journal of the American Statistical Association, 117(537), 348-360.

## Purpose

This package was created for personal research, methodological study, and comparative evaluation.

Its goal is to provide a standalone implementation of simulation, mean estimation, variance estimation, and covariance estimation procedures for functional snippet data.

This package does **not** import, wrap, or reuse source code from the authors’ original `mcfda` package. It is an independent reimplementation developed from the paper and methodological analysis.

## Main components

The package currently includes tools for:

- simulation of functional snippet data
- mean estimation
- variance component estimation
- covariance estimation
- performance evaluation of estimates

## Installation

For now, the package is under active development and can be installed from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("andrea-bianco00/mymcfda")
