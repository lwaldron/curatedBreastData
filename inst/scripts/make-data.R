#!/usr/bin/env Rscript

# Script to extract curatedBreastData objects and save them as an archive for Zenodo
# Run this from the package root directory.

# Load the original data
message("Loading original data...")
load("data/curatedBreastDataExprSetList.rda")
load("data/clinicalData.rda")

# Create a temporary directory to hold the raw files
tmp_dir <- tempfile("curatedBreastData_")
dir.create(tmp_dir)

# Define the paths for the extracted objects
exprset_path <- file.path(tmp_dir, "curatedBreastDataExprSetList.rda")
clinical_path <- file.path(tmp_dir, "clinicalData.rda")

# Save them into the temp directory
message("Saving objects to temporary directory...")
save(curatedBreastDataExprSetList, file = exprset_path, compress = "xz")
save(clinicalData, file = clinical_path, compress = "xz")

# Create the tar.gz archive
archive_name <- "curatedBreastData_archive.tar.gz"
message(sprintf("Creating archive %s...", archive_name))

# Tar it up from the temp dir so the paths inside the tarball are clean
old_dir <- setwd(tmp_dir)
tar(archive_name, files = c("curatedBreastDataExprSetList.rda", "clinicalData.rda"), compression = "gzip")
setwd(old_dir)

# Move the archive to the package root or current working dir
file.copy(file.path(tmp_dir, archive_name), archive_name, overwrite = TRUE)

message(sprintf("Done! The archive has been saved to %s", archive_name))
message("Please upload this file to Zenodo and record the URL.")

# Cleanup
unlink(tmp_dir, recursive = TRUE)
