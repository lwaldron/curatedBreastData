# curatedBreastData

## Curated Breast Cancer Transcriptomic Data

https://doi.org/10.18129/B9.bioc.curatedBreastData

**Bioconductor:**
[![BioC status](https://www.bioconductor.org/shields/build/release/data-experiment/curatedBreastData.svg)](https://bioconductor.org/checkResults/release/data-experiment-LATEST/curatedBreastData)
[![Platforms](https://bioconductor.org/images/shields/availability/all.svg)](https://bioconductor.org/packages/release/data/experiment/html/curatedBreastData.html#archives)
[![Downloads](https://www.bioconductor.org/shields/downloads/release/curatedBreastData.svg)](https://bioconductor.org/packages/stats/data-experiment/curatedBreastData)

**GitHub Actions:**
[![R CMD check](https://github.com/waldronlab/curatedBreastData/actions/workflows/ci.yml/badge.svg)](https://github.com/waldronlab/curatedBreastData/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/waldronlab/curatedBreastData/graph/badge.svg)](https://codecov.io/gh/waldronlab/curatedBreastData)

## Installation

We recommend installing the stable release version of `curatedBreastData` in Bioconductor. This can be done using `BiocManager`:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("curatedBreastData")
```

## Overview

The `curatedBreastData` package provides a collection of curated breast cancer transcriptomic datasets. See the package vignettes for examples on how to retrieve and post-process these datasets.
