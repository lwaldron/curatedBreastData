#' Curated breast gene expression data with survival and treatment information
#'
#' 34 manually curated high-quality gene expression microarray datasets with 
#' advanced beast cancer samples collected from GEO. All datasets provided have 
#' some form of survival and treatment information, and all such clinical variables
#' are semantically normalized across all datasest for easy analyses across 
#' datasets.  Authors of the Pubmed article linked to each GEO dataset was
#' contacted in an effort to collect as much extra clinical data as possible.
#' See vignette and publication reference from AMIA Translational Science Joint
#' Summits presentation in 2013 for more details on how this data was curated.  
#' 
#' Functions are provided to post-process standard S4 ExpressionSet objects 
#' to remove samples with high NA rates, impute missing values, collapse duplicated
#' gene symbols or probes, remove duplicated samples that share the same patient
#' ID, and filter genes by variance magnitude or percentile.
#'
#'
#' @author Katie Planey
#' @import XML
#' @import Biobase
#' @import BiocStyle
#' @importFrom ggplot2 ggplot ggsave geom_histogram aes labs
#' @importFrom impute impute.knn
#' @importFrom methods validObject new
#' @importFrom stats cov median na.omit quantile var
#' 
#' @keywords internal
"_PACKAGE"
#' @note Suggestions for new datasets to add are always welcome; the maintainer does aim
#' to only include datasets that have minimal treatment and some form of survival
#' (and/or treatment response) to allow for richer analyses. Raw data is always
#' preferred in order to control normalization schemes. Normalization details for
#' each dataset can be found the the Github repo in the References section.
#' @references Planey, Butte. Database integration of 4923 publicly-available samples of breast
#' cancer molecular and clinical data. AMIA Joint Summits Translational
#' Science Proceedings. (2003) PMC3814460
#' 
#' Github repo with code, further documentation on datasets and baseline
#' normalization schemes, and database quality checks: 
#' https://github.com/kplaney/curatedBreastCancer
#' @keywords package
#' @examples
#' \dontrun{
#' #load up datasets that are in S4 expressionSet format.
#' curatedBreastDataExprSetList <- getCuratedBreastDataExprSetList(test=TRUE)
#' 
#' #process only the first two datasets to avoid a long-running example:
#' #take top 5000 genes by variance from each dataset.
#' proc_curatedBreastDataExprSetList <- processExpressionSetList(exprSetList=
#' curatedBreastDataExprSetList[1:2], 
#' outputFileDirectory = "./", numTopVarGenes=5000)
#' 
#' #now we have processed expression matrices,
#' #each with the top 5000 genes by variance 
#' }


#' @title curatedBreastData_pData
#'
#' @description A single \code{data.frame} containing the concatenated \code{pData} from all datasets in the \code{curatedBreastDataExprSetList} list. Rownames correspond to sample IDs. 
#'
#' @format A \code{data.frame} with patient samples as rows and phenotype variables as columns.
#' @details This data object was created by combining the \code{pData} slot of all ExpressionSets provided in the \code{curatedBreastDataExprSetList}.
#' @examples
#' data(curatedBreastData_pData)
#' dim(curatedBreastData_pData)
#' head(curatedBreastData_pData)
#' @keywords datasets
"curatedBreastData_pData"
