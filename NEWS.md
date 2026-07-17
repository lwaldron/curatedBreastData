# curatedBreastData 2.40.1

## SIGNIFICANT USER-VISIBLE CHANGES

- Package data has been successfully offloaded to Zenodo (https://zenodo.org/records/14842512) and is no longer directly embedded in the package. The data is now dynamically retrieved and cached via `BiocFileCache`.

## NEW FEATURES

- Added new retrieval functions `getCuratedBreastDataExprSetList()` and `getClinicalData()` to access the remote dataset from Zenodo.
- Embedded a miniature test dataset archive in `inst/extdata` for rapid local testing and offline execution without remote dependencies.

## BUG FIXES

- Fixed operator precedence bugs and replaced deprecated `1:...` loop sequence generators with type-safe `seq_along()` and `seq_len()`.
- Converted `sapply` calls to type-safe `vapply()` in `collapseDupProbes()`.
- Resolved redundant error/warning signal conditions in processing handlers.
- Updated `DESCRIPTION` metadata to include the latest R dependency (`>= 4.6.0`), missing `biocViews`, `URL`, and `BugReports` fields.
- Fixed `BiocFileCache` path query matching on character-coerced resource IDs.

## INTERNAL CHANGES

- Migrated all package documentation to `roxygen2` and removed legacy `.Rd` manual pages.
- Converted the package vignette from Sweave (`.Rnw`) to R Markdown (`.Rmd`) for improved readability and compilation performance.
- Added `README.md` and `_pkgdown.yml` mirroring the design of `MultiAssayExperiment`.
- Configured GitHub Actions CI workflows for automated `R CMD check` and Bioconductor compliance validation (`BiocCheck`).

# curatedBreastData 0.99.0

- Initial package release.
