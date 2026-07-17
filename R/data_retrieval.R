#' getCuratedBreastDataExprSetList
#' 
#' A list of ExpressionSet objects, one for each curated study, containing
#' study-specific gene expression and phenotype data. FeatureNames are gene
#' symbols. Data is already quantile normalized according to standard protocols
#' for 1 and 2-channel arrays, depending on the platform used for this study.
#' 
#' Batches from studies are treated as individual datasets, as the signal
#' can differ between batches. Thus, an expression object named using a GSE
#' study number followed by an underscore means this ExpressionSet contains
#' samples either from a distinct platform (and the study used >1 platforms),
#' or from a distinct batch or tissue site.  An "all" tag means that there were no
#' batches for this study. Raw data files downloaded from GEO oftentimes have clear
#' batch/site information appended to sample names; this was often the source of
#' batch identification and how the package developer chose to create the batch
#' name string.
#' 
#' @param test A boolean. If TRUE, downloads a miniature test dataset instead of the full dataset. Default is FALSE.
#' @return A list, with each index containing an ExpressionSet object from a specific study, and potentially a specific batch.
#' @references Planey, Butte. Database integration of 4923 publicly-available samples of breast cancer molecular and clinical data. AMIA Joint Summits Translational Science Proceedings. (2003) PMC3814460
#' @examples
#' curatedBreastDataExprSetList <- getCuratedBreastDataExprSetList(test=TRUE)
#' #what are all the names of the studies?
#' names(curatedBreastDataExprSetList)
#' #what is the dimension of the gene 
#' #expression matrix for study GSE2034?
#' dim(exprs(curatedBreastDataExprSetList$study_2034_GPL96_all))
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd bfcrpath
#' @importFrom utils untar
#' @export

getCuratedBreastDataExprSetList <- function(test = FALSE) {
  if (test) {
    zenodo_url <- "https://zenodo.org/records/21415886/files/curatedBreastData_test_archive.tar.gz?download=1"
    archive_name <- "curatedBreastData_test_archive.tar.gz"
  } else {
    zenodo_url <- "https://zenodo.org/records/21415886/files/curatedBreastData_archive.tar.gz?download=1"
    archive_name <- "curatedBreastData_archive.tar.gz"
  }
  
  bfc <- BiocFileCache()
  res <- bfcquery(bfc, archive_name, "rname", exact = TRUE)
  
  if (nrow(res) == 0) {
    message("Downloading curatedBreastData archive from Zenodo...")
    rpath <- bfcadd(bfc, archive_name, zenodo_url)
  } else {
    rpath <- bfcrpath(bfc, archive_name)
  }
  
  # The archive contains curatedBreastDataExprSetList.rda and clinicalData.rda
  # We extract them to a temporary directory
  tmp_dir <- tempdir()
  untar(rpath, exdir = tmp_dir)
  
  # Load the data
  rda_path <- file.path(tmp_dir, "curatedBreastDataExprSetList.rda")
  if (!file.exists(rda_path)) {
    stop("Archive does not contain curatedBreastDataExprSetList.rda")
  }
  
  # Load it into an environment so we can return it
  env <- new.env()
  load(rda_path, envir = env)
  
  return(env$curatedBreastDataExprSetList)
}

#' Clinical Data Table & Variable Definitions
#'
#' Clinical data for all samples across all studies, and corresponding variable 
#' definitions. Rownames are the GEO_GSMID feature, which corresponds to the sample 
#' names in the expression object for a certain study. Includes treatment 
#' information.
#' 
#' GEO study ID can be found form the study_ID variable. If site_ID is NA, it 
#' pertains to the batch ID, which may be due to different platforms being used 
#' in the same study or different tissue site collections. Columns 112-151 pertaint
#' to treatment information. radiotherapyClass,  chemotherapyClass, and 
#' hormone_therapyClass are indicator variables used to signal whether a patient 
#' had radiotherapy, chemotherapy, and/or some form of hormone therapy 
#' (usually an estrogen or aromatse inhibitor.) 
#' 
#' More granular information, when available, is provided: for example, whether 
#' the chemotherapy drug was capecitabine is coded as the indicator "capecitabine" 
#' variable. A value of 1 = yes, 0 = no, NA = not recorded/could not infer from 
#' publically available information. "Other" means that most likely, gleaned from 
#' the study's Pubmed publication, that the patient may have had other treatments 
#' that were not recorded (oftentimes radiotherapy, as this is not always recorded 
#' and up to a clinician's discretion in a clinical trial.)
#' 
#' Survival information, such as DFS, RFS, OS, and treatment response information, 
#' such pCR and RCB, is also recorded when available.
#'
#' @param test A boolean. If TRUE, downloads a miniature test dataset instead of the full dataset. Default is FALSE.
#' @return A list with the following two items:
#' \item{clinicalTable}{A data frame. Rownames are the GEO_GSMID feature, which corresponds to the sample names in the expression object for a certain study.}
#' \item{clinicalVarDef}{Character string descriptions of each variable.}
#' @references Planey, Butte. Database integration of 4923 publicly-available samples of breast cancer molecular and clinical data. AMIA Joint Summits Translational Science Proceedings. (2003) PMC3814460
#' @examples
#' clinicalData <- getClinicalData(test=TRUE)
#' #check out some of the variable name/definitions
#' clinicalData$clinicalVarDef[c(1:2),]
#' #Check out the treatment information. 
#' #look at first three patients
#' head(clinicalData$clinicalTable)[c(1:3),c(112:ncol(clinicalData$clinicalTable))]
#' #how many had chemotherapy?
#' numChemoPatients <- length(which(
#' clinicalData$clinicalTable$chemotherapyClass==1))
#' #how many patients have non-NA OS binary data?
#' length(which(!is.na(clinicalData$clinicalTable$OS)))
#' #how many have OS data in the more granular form of months until OS? 
#' #this variable includes studies that had a cieling for tracking OS
#' length(which(!is.na(clinicalData$clinicalTable$OS_months_or_MIN_months_of_OS)))
#' #how many patients have OS information that is definitively 
#' #followed up until their death
#' #(details on how studies collect OS data can be surprising!)
#' length(which(!is.na(clinicalData$clinicalTable$OS_up_until_death)))
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd bfcrpath
#' @importFrom utils untar
#' @export

getClinicalData <- function(test = FALSE) {
  if (test) {
    zenodo_url <- "https://zenodo.org/records/21415886/files/curatedBreastData_test_archive.tar.gz?download=1"
    archive_name <- "curatedBreastData_test_archive.tar.gz"
  } else {
    zenodo_url <- "https://zenodo.org/records/21415886/files/curatedBreastData_archive.tar.gz?download=1"
    archive_name <- "curatedBreastData_archive.tar.gz"
  }
  
  bfc <- BiocFileCache()
  res <- bfcquery(bfc, archive_name, "rname", exact = TRUE)
  
  if (nrow(res) == 0) {
    message("Downloading curatedBreastData archive from Zenodo...")
    rpath <- bfcadd(bfc, archive_name, zenodo_url)
  } else {
    rpath <- bfcrpath(bfc, archive_name)
  }
  
  # The archive contains curatedBreastDataExprSetList.rda and clinicalData.rda
  # We extract them to a temporary directory
  tmp_dir <- tempdir()
  untar(rpath, exdir = tmp_dir)
  
  # Load the data
  rda_path <- file.path(tmp_dir, "clinicalData.rda")
  if (!file.exists(rda_path)) {
    stop("Archive does not contain clinicalData.rda")
  }
  
  # Load it into an environment so we can return it
  env <- new.env()
  load(rda_path, envir = env)
  
  return(env$clinicalData)
}
