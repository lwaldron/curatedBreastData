#!/usr/bin/env Rscript

# Script to extract a miniature subset of curatedBreastData objects and save them as an archive for local testing
# Run this from the package root directory.

# Load the original data (you must have access to the original large rda files or the downloaded archive)
message("Loading data...")
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to run inst/scripts/make-test-data.R", call. = FALSE)
}
devtools::load_all()
# Load the full dataset (remote/cache) before subsetting
esets <- getCuratedBreastDataExprSetList(test = FALSE)
clin <- getClinicalData(test = FALSE)

# Create a small subset
message("Subsetting data...")
# Subset to first 500 features and first 30 samples for both studies
eset1 <- esets[[1]][1:500, 1:30]
eset2 <- esets[[2]][1:500, 1:30]
curatedBreastDataExprSetList <- list(study_1379_GPL1223_all = eset1, study_2034_GPL96_all = eset2)

clinicalData <- clin
target_samples <- c(sampleNames(eset1), sampleNames(eset2))
clinicalData$clinicalTable <- clinicalData$clinicalTable[clinicalData$clinicalTable$GEO_GSMID %in% target_samples, ]

# Create a temporary directory to hold the raw files
tmp_dir <- tempfile("curatedBreastData_test_")
dir.create(tmp_dir)

# Define the paths for the extracted objects
exprset_path <- file.path(tmp_dir, "curatedBreastDataExprSetList.rda")
clinical_path <- file.path(tmp_dir, "clinicalData.rda")

# Save them into the temp directory
message("Saving objects to temporary directory...")
save(curatedBreastDataExprSetList, file = exprset_path, compress = "xz")
save(clinicalData, file = clinical_path, compress = "xz")

# Create the tar.gz archive in inst/extdata
extdata_dir <- "inst/extdata"
if (!dir.exists(extdata_dir)) {
  dir.create(extdata_dir, recursive = TRUE)
}
archive_name <- "curatedBreastData_test_archive.tar.gz"
dest_archive_path <- file.path(extdata_dir, archive_name)
message(sprintf("Creating archive at %s...", dest_archive_path))

# Tar it up from the temp dir so the paths inside the tarball are clean
old_dir <- setwd(tmp_dir)
tar(archive_name, files = c("curatedBreastDataExprSetList.rda", "clinicalData.rda"), compression = "gzip")
setwd(old_dir)

# Copy the archive to inst/extdata
file.copy(file.path(tmp_dir, archive_name), dest_archive_path, overwrite = TRUE)

message(sprintf("Done! The test archive has been saved to %s", dest_archive_path))
message(sprintf("Archive file size is: %d bytes", file.info(dest_archive_path)$size))

# Cleanup
unlink(tmp_dir, recursive = TRUE)
