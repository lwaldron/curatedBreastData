#' Process a list of S4 expressionSet objects.
#'
#' A wrapper function for the post-processing function processExpressionSet() on
#' a list of S4 expressionSet objects. This function is run after initial dataset
#' normalization, such as quantile normalization on microarray datasets.
#'
#' @param exprSetList List of S4 expression sets.
#' @param outputFileDirectory Output file directory for messages that print status of post-processing the ExpressionSets.
#' @param numTopVarGenes A numeric value indicating the number of genes (features) to select; the function will only take this number of genes that have the highest variance across all genes.
#' @param minVarPercentile Minimum variance percentile. Must be provided in conjunction with maxVarPercentile to use percentiles to threshold genes.
#' @param maxVarPercentile Maximum variance percentile. Defaul is 1, i.e. 1\%. Must be provided in conjunction with minVarPercentile to use percentiles to threshold genes.
#' @param minVar If maxVar is provided, as opposed to minVarPercentile and maxVarPercentile, genes are removed that are below a certain variance magnitude. This is helpful before running certain algorithms, such as the popular Combat batch normalization technique, that can throw errors if genes with extremely low variances are in the data matrix. May be used in conjunction with maxVar or in isolation.
#'
#' @return A list of processed S4 ExpressionSet objects.
#'
#' @author Katie Planey <katie.planey@@gmail.com>
#' @seealso \code{\link{processExpressionSet}} 
#' @examples
#' #load up our datasets
#' curatedBreastDataExprSetList <- getCuratedBreastDataExprSetList(test=TRUE)
#' 
#' #just take top 100 genes by variance
#' #this will post-process every dataset in the package
#' #to make them ready for downstream analyses.
#' proc_curatedBreastDataExprSetList <- processExpressionSetList(
#' exprSetList=curatedBreastDataExprSetList, 
#' outputFileDirectory = tempdir(), numTopVarGenes=100)
#' @export
processExpressionSetList <- function(exprSetList,outputFileDirectory="./",
                                     numTopVarGenes,minVarPercentile,maxVarPercentile=1,minVar){
  
  outputFile = paste0(outputFileDirectory,
                      "/curatedBreastData_processExpressionSetMessages.txt")

  for(e in 1:length(exprSetList)){
    message("\nAnalyzing dataset ",e, " or dataset named ",
            names(exprSetList)[e]);
    
    if(!missing(numTopVarGenes)){
      
     exprSetList[[e]] <- processExpressionSet( exprSetList[[e]],
                        outputFileDirectory=outputFileDirectory,
                        numTopVarGenes=numTopVarGenes)
    
    }else if(!missing(minVarPercentile) && !missing(maxVarPercentile) 
             && missing(minVar)){
      
      exprSetList[[e]] <- processExpressionSet( exprSetList[[e]],
                          outputFileDirectory=outputFileDirectory,
                          maxVarPercentile=maxVarPercentile,
                          minVarPercentile=minVarPercentile)
      
    }else if(!missing(minVar)){
      
      exprSetList[[e]] <- processExpressionSet( exprSetList[[e]],
                          outputFileDirectory=outputFileDirectory,minVar=minVar)
      
      
   }else{
      
      exprSetList[[e]] <- processExpressionSet(exprSetList[[e]],
                          outputFileDirectory=outputFileDirectory)
    
    }
   
  }
  
  return(exprSetList)
  
}

#assumption:  your feature data has a column that labeled gene_symbol.
#' Post-process a normalized assayData in an ExpressionSet object
#'
#' This function Post-processes a normalized assayData in an ExpressionSet object.
#' Is is assumed the assay data is already baseline normalized (for example, for
#' microarray data, this could mean quantile normalized and then logged.)
#'
#' This function performs several post-processing tasks: filtering out genes and
#' samples with high NA rates, imputing missing values, collapsing duplicated
#' features/genes to make a unique feature list, removing any samples for which
#' there is already a sample with the sample patient ID (duplicated samples),
#' and filtering genes by variance.  This function is a wrapper for the functions:
#' filterAndImputeSamples(), collapseDupProbes(), removeDuplicatedPatients(), and
#' filterGenesByVariance(). It is is run after initial dataset normalization, such
#' as quantile normalization on microarray datasets.
#'
#' @param exprSet expressionSet S4 object with expression (assay) data, featureData and phenoData.
#' @param outputFileDirectory Output file directory for messages that print status of post-processing the ExpressionSet.  
#' @param numTopVarGenes A numeric value indicating the number of genes (features) to select; the function will only take this number of genes that have the highest variance across all genes.
#' @param minVarPercentile Minimum variance percentile. Must be provided in conjunction with maxVarPercentile to use percentiles to threshold genes.
#' @param maxVarPercentile Maximum variance percentile. Defaul is 1, i.e. 1\%. Must be provided in conjunction with minVarPercentile to use percentiles to threshold genes.
#' @param minVar If maxVar is provided, as opposed to minVarPercentile and maxVarPercentile, genes are removed that are below a certain variance magnitude. This is helpful before running certain algorithms, such as the popular Combat batch normalization technique, that can throw errors if genes with extremely low variances are in the data matrix. May be used in conjunction with maxVar or in isolation.
#'
#' @return A post-processed S4 expressionSet.  Tests are run to confirm the final S4 object is a valid ExpressionObject before it is returned.
#'
#' @author Katie Planey <katie.planey@@gmail.com>
#' @examples
#' #load up our datasets
#' curatedBreastDataExprSetList <- getCuratedBreastDataExprSetList(test=TRUE)
#' 
#' #just perform on one dataset as an example, GSE1379. 
#' #take only genes that fall in 
#' #the variance percentiles between .75 and 1 
#' #(i.e. top 75th percentile genes by variance.)
#' 
#' post_procExprSet <- processExpressionSet(exprSet=
#' curatedBreastDataExprSetList[[1]], 
#' outputFileDirectory = tempdir(),
#' minVarPercentile=.75, maxVarPercentile = 1)
#' @export
processExpressionSet <- function(exprSet,outputFileDirectory="./",
                                 numTopVarGenes,minVarPercentile,
                                 maxVarPercentile=1,minVar){
  
  study <- list(expr=exprSet@assayData$exprs,
                keys=exprSet@featureData$gene_symbol,phenoData=pData(exprSet))
  tmp <- filterAndImputeSamples(study,studyName = "study",
         outputFile = paste0(outputFileDirectory,
         "/curatedBreastData_processExpressionSetMessages.txt"),impute=TRUE, 
          knnFractionSize=.01,fractionSampleNAcutoff=.005,
         fractionGeneNAcutoff = .01,exprIndex="expr",classIndex="phenoData",
         sampleCol=TRUE,returnErrorRate=FALSE)
    
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))
    
    tmpExpr <- data.matrix(tmp$exprFilterImpute)
    exprSet@assayData <- assayDataNew(exprs =  tmpExpr)
    featureNames(exprSet@assayData) <- c(1:length(tmp$keysFilterImpute))
    #need a data.frame and not cbind: otherwise coerces 
    #gene symbols into numeric data type.
    fData <- data.frame(c(1:length(tmp$keysFilterImpute)), tmp$keysFilterImpute)
    colnames(fData) <- c("number","gene_symbol")
    fData(exprSet) <- fData
  #only one set of classes data here - the phenoData, 
  #so will be in first list index.
    rownames(tmp$classesFilter[[1]]) <- colnames(exprs(exprSet))
    pData(exprSet) <- data.frame(tmp$classesFilter[[1]])
  #must also set global featureNames to updated keys, not just featureData.
  #otheriwse won't return a valid ExpressionObject.
  featureNames(exprSet@featureData) <- c(1:length(tmp$keysFilterImpute))
    # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
  if(!validObject(exprSet)){
    
    warning("\nYour expression set is now not valid.
This is usually due to sample or feature names not matching up across
assay, pheno and feature slots.\n Proceed through data analysis with caution!")
    
  }
  #collapse duplicated probes (really gene symbols here)
   tmp <-  collapseDupProbes(expr=exprSet@assayData$exprs,
            sampleColNames=colnames(exprSet@assayData$exprs),
            keys=fData(exprSet)[[2]], method="highestVariance",
            debug=TRUE,removeNA_keys=TRUE,varMetric = "everything")


    phenoData <- pData(exprSet)
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))
    
    tmpExpr <- data.matrix(tmp$expr)
    fData <- data.frame(tmp$keys,stringsAsFactors=FALSE)
    colnames(fData) <- "gene_symbol"
    fData(exprSet) <- fData
    pData(exprSet) <- phenoData
    #must also set global featureNames to updated keys, not just featureData.
    #otheriwse won't return a valid ExpressionObject.
    featureNames(exprSet@featureData) <- tmp$keys
    exprSet@assayData <- assayDataNew(exprs =  tmpExpr)
    featureNames(exprSet@assayData) <- tmp$keys
    # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
  
  tmp <-  removeDuplicatedPatients(exprMatrix=exprSet@assayData$exprs, 
          outputFile=paste0(outputFileDirectory,
          "/curatedBreastData_processExpressionSetMessages.txt"), 
          varMetric = "everything")

  if(ncol(exprs(exprSet))!=ncol(tmp)){
  
    exprs(exprSet) <- tmp
    featureN <- featureNames(exprSet@featureData)
    #lost some samples. need to subset the matrix
    phenoData <- pData(exprSet)[ ,na.omit(match(colnames(tmp),
                                    rownames(pData(exprSet))))]
    fData <- fData(exprSet)
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))

    fData(exprSet) <- fData
    pData(exprSet) <- phenoData
    #must also set global featureNames to updated keys, not just featureData.
    #otheriwse won't return a valid ExpressionObject.
    featureNames(exprSet@featureData) <- featureN
    exprSet@assayData <- assayDataNew(exprs =  tmp)
    featureNames(exprSet@assayData) <- featureN
    # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
        # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
    
  }

  study <- list(expr=exprSet@assayData$exprs,
                keys=exprSet@featureData$gene_symbol)
  
  if(!missing(numTopVarGenes)){
    
    tmp <- filterGenesByVariance(study, plotSaveDir="~/",
          exprIndex = "expr", keysIndex = "keys", 
          outputFile=paste0(outputFileDirectory,
          "/curatedBreastData_processExpressionSetMessages.txt"),
          varMetric = c("everything"),sampleCol=TRUE,
          numTopVarGenes=numTopVarGenes,plotVarianceHist=FALSE);
    

       tmpExpr <- data.matrix(tmp$filteredStudy$expr)
    #if lost some samples. need to subset the matrix
    phenoData <- pData(exprSet)[na.omit(match(colnames(tmpExpr),
                                    rownames(pData(exprSet)))), ]
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))

    fData <- data.frame(tmp$filteredStudy$keys,stringsAsFactors=FALSE)
    colnames(fData) <- "gene_symbol"
    fData(exprSet) <- fData
    pData(exprSet) <- phenoData
    #must also set global featureNames to updated keys, not just featureData.
    #otheriwse won't return a valid ExpressionObject.
    featureNames(exprSet@featureData) <- tmp$filteredStudy$keys
    exprSet@assayData <- assayDataNew(exprs =  tmpExpr)
    featureNames(exprSet@assayData) <- tmp$filteredStudy$keys
   

    # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
        # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
    #do quantiles (percentage cut-offs.)
  }else if(!missing(minVarPercentile) && !missing(maxVarPercentile)
           && missing(minVar)){
    
    tmp <- filterGenesByVariance(study, plotSaveDir="~/",
          minVarPercentile=minVarPercentile,maxVarPercentile=maxVarPercentile,
          exprIndex = "expr", keysIndex = "keys", 
          outputFile=paste0(outputFileDirectory,
          "/curatedBreastData_processExpressionSetMessages.txt"),
          varMetric = c("everything"),sampleCol=TRUE,plotVarianceHist=FALSE);
    
     tmpExpr <- data.matrix(tmp$filteredStudy$expr)
    #if lost some samples. need to subset the matrix
    phenoData <- pData(exprSet)[na.omit(match(colnames(tmpExpr),
                                    rownames(pData(exprSet)))), ]
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))

    fData <- data.frame(tmp$filteredStudy$keys,stringsAsFactors=FALSE)
    colnames(fData) <- "gene_symbol"
    fData(exprSet) <- fData
    pData(exprSet) <- phenoData
    #must also set global featureNames to updated keys, not just featureData.
    #otheriwse won't return a valid ExpressionObject.
    featureNames(exprSet@featureData) <- tmp$filteredStudy$keys
    exprSet@assayData <- assayDataNew(exprs =  tmpExpr)
    featureNames(exprSet@assayData) <- tmp$filteredStudy$keys

    # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
        # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
    
  }else if(!missing(minVar)){
    
    tmp <- filterGenesByVariance(study, plotSaveDir="~/",minVar=minVar,
           exprIndex = "expr", keysIndex = "keys", 
           outputFile=paste0(outputFileDirectory,
            "/curatedBreastData_processExpressionSetMessages.txt"),
            varMetric = c("everything"),sampleCol=TRUE,plotVarianceHist=FALSE);
    
     tmpExpr <- data.matrix(tmp$filteredStudy$expr)
    #if lost some samples. need to subset the matrix
    phenoData <- pData(exprSet)[na.omit(match(colnames(tmpExpr),
                                    rownames(pData(exprSet)))), ]
    exprSet <- new("ExpressionSet", assayData = assayDataNew(exprs=new("matrix")), 
                   phenoData = new("AnnotatedDataFrame"), 
                   featureData = new("AnnotatedDataFrame"), 
                   experimentData = new("MIAME"), annotation = character(0))

    fData <- data.frame(tmp$filteredStudy$keys,stringsAsFactors=FALSE)
    colnames(fData) <- "gene_symbol"
    fData(exprSet) <- fData
    pData(exprSet) <- phenoData
    #must also set global featureNames to updated keys, not just featureData.
    #otheriwse won't return a valid ExpressionObject.
    featureNames(exprSet@assayData) <- tmp$filteredStudy$keys
    featureNames(exprSet@featureData) <- tmp$filteredStudy$keys
    exprSet@assayData <- assayDataNew(exprs =  tmpExpr)

    # #may need to re-set protocolData field so dimensions match
   labelDescription <- rep("Breast cancer human tumor tissue sample", nrow(pData(exprSet)))
   labelDescription <- data.frame(labelDescription) 
   rownames(labelDescription) <- rownames(pData(exprSet))
   protocolData(exprSet) <- AnnotatedDataFrame(labelDescription)
    
  }
    
  #test: still a valid expression set??
    if(!validObject(exprSet)){
      
      warning("\nYour expression set is now not valid.
This is usually due to sample or feature names not matching up across assay, 
pheno and feature slots.\n Proceed through data analysis with caution!")
      
    }

    return(exprSet)

}
#fractionGeneNAcutoff:  max fraction of NAs allowed for a certain gene across all samples.
#fractionSampleNAcutoff: max fraction of NAs allowed for a certain sample across all genes.
#' Filter and Impute Samples
#'
#' A method that removes samples or genes with high NA rates and then KNN imputes
#' remaining missing values.
#'
#' @param study A list, of minimally the gene expression or some molecular data matrix with keys (molecular features, such as genes) in the rows and patient samples in the columns and a keys list. It is assumed that the keys entity is named "keys", but in line with using this function for any type of molecular data, the exprIndex name in the list can be altered.
#' @param studyName Character string to name the study. Useful in cases where looping over multiple datasets; output messages printed to the output file can then be identified by each individual study name. 
#' @param outputFile Output File for printing progress and stats on gene/sample filtering and data imputation.  Include full directory if file should not be printed to current working directory.
#' @param impute Impute data? A boolean TRUE or FALSE value. If FALSE, only genes and samples with high NA rates are removed, and the rest of the data is not imputed.
#' @param knnFractionSize What is the fraction of neighbors out of the total dataset to be used for knn impute nearest neighbor? This is translated into the "k" numeric magnitude in impute.knn() from the impute package. Default is .01, or 1\% of the data.
#' @param fractionSampleNAcutoff Max fraction of NAs allowed for a certain sample across all genes. Default is .005 (.005, or .5\%, still captures a large number of gnees for a sample if there are tens of thousands of genes in the data matrix.)
#' @param fractionGeneNAcutoff Max fraction of NAs allowed for a certain gene across all samples. Default is .01. Thus, a certain gene cannot be missing in greater than 1\% of patients. It is recommended that this threshold be increased for smaller datasets unless a user wants a gene to be removed that is missing in only 1 sample.
#' @param exprIndex Character string. List slot name for the data matrix, presumably an expression matrix.
#' @param classIndex Optional character string giving the list slot name for a phenotype vector or matrix if available. If phenotype/class data such as survival is already in the list, filtering out samples with high NA rates will result in the need to remove these samples from the phenotype data matrix; filterAndImputeSamples will appropriately filter out these samples from the phenoteyp data.
#' @param sampleCol Are samples in the columns of the expression matrix? If not, this function will first transpose the matrix to make sure impute.knn is running properly.
#' @param returnErrorRate Boolean TRUE or FALSE. If TRUE, a small amount of real expression data points are held out, and knn.impute is performed. The accuracy rate of the imputed values vs. the real values is returned. THis is helpful in early data analysis stages to determine whether KNN imputation is appropriate for your type of data. Default is FALSE to reduce computation time.
#'
#' @return A list containing the following objects:
#' \item{expr}{original expression matrix}
#' \item{exprFilterImputed}{final filtered and imputed expression matrix}
#' \item{keys}{original keys}
#' \item{keys}{final filtered and imputed keys}
#' \item{classes}{original classes/phenotype data}
#' \item{classes}{final classes/phenotype data, removing any sample rows that were removed from the expression matrix after filtering.}
#'
#' @author Katie Planey <katie.planey@@gmail.com>
#' @examples
#' #load up our datasets
#' curatedBreastDataExprSetList <- getCuratedBreastDataExprSetList(test=TRUE)
#' 
#' #just perform on one dataset as an example, GSE1379. 
#' #create study list object. 
#' study <- list(expr=exprs(curatedBreastDataExprSetList[[1]]),
#' keys=fData(curatedBreastDataExprSetList[[1]])[, "gene_symbol"],
#' phenoData=pData(curatedBreastDataExprSetList[[1]]))
#' 
#' filteredStudy <- filterAndImputeSamples(study, studyName = "study", 
#' outputFile = tempfile(), impute = TRUE, 
#' knnFractionSize = 0.01, fractionSampleNAcutoff = 0.005, 
#' fractionGeneNAcutoff = 0.01, exprIndex = "expr", classIndex="phenoData",
#' sampleCol = TRUE, returnErrorRate = TRUE)
#' 
#' #see output list names 
#' names(filteredStudy)
#' #what is the imputation error fraction (rate)?
#' filteredStudy$errorRate
#' @export
filterAndImputeSamples <- function(study,studyName = "study",
                                  outputFile = "createTestTrainSetsOutput.txt",
                                  impute=TRUE, knnFractionSize=.01,
                                  fractionSampleNAcutoff=.005,
                                  fractionGeneNAcutoff = .01,
                                  exprIndex="expr",classIndex,sampleCol=TRUE,
                                  returnErrorRate=TRUE){

  #   #knn impute format: An expression matrix with genes 
  #in the rows, samples in the columns
  if(sampleCol){
    #case were samples are the columns, not the rows.
    
    expr <- study[[exprIndex]]
    exprOrig <- study[[exprIndex]]
    
  }else{
    
    expr <- t(study[[exprIndex]])
    exprOrig <- t(study[[exprIndex]])
    message("dimensions of expression study will be returned transposed: 
samples are now the columns for an pxn matrix.")
    
  }  

  #NOTE: must have a "keys" list in here!
  keysOrig <-study$keys
  numPatients <- dim(expr)[2]
  
  if(is.null(expr) || is.null(keysOrig) || all(is.na(expr)) || all(is.na(keysOrig))){
    
    stop("you didn't provide the right expr or keys keyword 
         for the filter genes by variance function.")
    
  }
  
  totalGen <- dim(expr)[1]
  
  message("\nNote: this function assumes your missing values
  are proper NAs, not \"null\",etc.\n")
  
  gene_fractionNAsamples <- apply(expr,MARGIN=1, 
                                  FUN = function(studyRow,numPatients){
    
    numNA <- length(which(is.na(studyRow)))
    fractionNAsamples <- numNA/numPatients
    return(fractionNAsamples)
    
    #pass in extra arguments to the function here
  },numPatients)
  
  #chances are, more cases where a ton of patients missing a 
  #certain gene than a ton of really bad patients with tons of NAs 
  #across lots of genes. so filter out these certain NA genes first.    
  goodGeneIndices <- which(gene_fractionNAsamples < fractionGeneNAcutoff)
  cat("there were ",length(goodGeneIndices), "genes that passed your max", 
      fractionGeneNAcutoff, "NA filter out of a total of ",totalGen, 
      "genes","\n",append=TRUE,file=outputFile)
  
  expr <- expr[goodGeneIndices,]
  #also need to update keys
  keysFilterImpute <- study$keys[goodGeneIndices]
  
  sample_fractionNAgenes <- apply(expr,MARGIN=2, 
                                  FUN = function(studyCol,totalGenes){
    
    numNA <- length(which(is.na(studyCol)))
    fractionNAgenes <- numNA/totalGenes
    return(fractionNAgenes)
    
    #pass in extra arguments to the function here
  },totalGen)
  
  goodIndices <- which(sample_fractionNAgenes < fractionSampleNAcutoff)   
  #now look at "bad" patient indices for training set
  
  cat("there were ",length(goodIndices), "patients that passed your max", 
      fractionSampleNAcutoff, "NA filter out of a total of ",numPatients, 
      "patients","\n",append=TRUE,file=outputFile)
  #cut out these patients from your already filtered gene list.
  expr <- expr[,goodIndices]
  
  #NOW run KNN impute
  #come back...what if it's a small sample size?
  knnFractionSize <- .1
  knnImpute <- function(knnFractionSize, expr,fractionSampleNAcutoff, 
                        fractionGeneNAcutoff,maxp=1500,rng.seed= 362436069){
    
    #k is the # samples use for neighborhood - try about 5% of patients?
    knnSize <- round(knnFractionSize*(dim(expr)[2]))
    #function below returns a few things - just keep them all, 
    #put out $data after function is run.
    exprFilterImpute <- impute.knn(expr,k=knnSize,colmax=fractionSampleNAcutoff,
                        rowmax=fractionGeneNAcutoff,maxp=maxp,rng.seed=rng.seed)
    return(exprFilterImpute)
    
  }
  
  if(impute){
    
    #what genes have at least 1 NA value?
    
    genesNA <- apply(expr,MARGIN=1, FUN = function(studyRow){
      
      return(any(is.na(studyRow)))
      
      #pass in extra arguments to the function here
    })
    
    #odd... if just try to index with -genesNA get nowhere...
    numNAgenes <- length(which(genesNA==TRUE))
    
    if(numNAgenes==0){
      
      cat("there are zero NA genes in this dataset 
          so returning original study \n",file=outputFile,append=TRUE)
      
      exprFilterImpute <- expr
      
      #there are some missing genes, continue on to impute!       
    }else{
      
      cat("there are ",numNAgenes, "genes with at least one NA value \n",
          file=outputFile,append=TRUE)
      
      #expr must NOT be a list!
      expr <- data.matrix(expr)
      exprFilterImpute <- knnImpute(k=knnFractionSize, expr=expr,
                          fractionSampleNAcutoff, fractionGeneNAcutoff,
                          maxp=1500,rng.seed= 362436069)
      exprFilterImpute <- exprFilterImpute$data
      
    }
    
    
    #don't bother with this if there aren't any NA genes.
    if(numNAgenes > 0 && returnErrorRate){
      
      noNAexpr <- expr[which(genesNA == FALSE),]
      
      #now try making some of these NA values, for each sample (randomly.)
      #mimic mean of actual NA genes missing in each sample.
      #average NA fraction
      fractionNA <- mean(sample_fractionNAgenes[
        which(sample_fractionNAgenes >0)])
      #how many NA genes per sample (column)?
      NAperCol <- round(dim(noNAexpr)[1]*fractionNA)
      
      
      fakeNAexpr <- apply(noNAexpr,MARGIN=2, FUN = function(studyCol, NAperCol){
        
        #make random genes (rows) for this patient NA
        NArows <- sample(1:length(studyCol),NAperCol,replace=FALSE)
        
        studyCol[NArows] <- NA
        
        return(studyCol)
        
        #pass in extra arguments to the function here.
        #match the limit for sample of fractions, or the average?
      },NAperCol)
      
      #now re-run imputation.
      #possible that randomly, for each sample, the same gene
      #was removed...so let that cutoff be 1.
      fakeNAindices <- which(is.na(fakeNAexpr))
      #add a .01 "buffer" to the fractionNA so the function doesn't 
      #throw an error, as our % column NA matches exactly fractionNA.
      fakeexprFilterImpute <- knnImpute(knnFractionSize, expr=fakeNAexpr,
                              fractionSampleNAcutoff=(fractionNA+.01), 
                              fractionGeneNAcutoff=1,maxp=1500,
                              rng.seed= 362436069)
      
      imputedNAs <- fakeexprFilterImpute$data[fakeNAindices]
      trueValues <- noNAexpr[fakeNAindices]
      errorRate <- sum(abs(trueValues-imputedNAs)/trueValues)/
        length(fakeNAindices)
      meanAbsDiff <- sum(abs(trueValues-imputedNAs))/length(fakeNAindices)
      #COME BACK: make more flexible - not just knn impute?
      #get NAs by comparing fakeNAexpr and noNAexpr.  compare output 
      #of imputation and noNAexpr.
      cat("the average absolute error for imputation for study ", studyName,
" is ", errorRate, " with an average absolute difference of log-fold
          change expression of ", meanAbsDiff, "\n","\n", 
file=outputFile,append=TRUE)
      
      #now calculate error. take only non-NA genes. 
      
      
      
    }else{
      
      cat("error rate not specified to be calculated or 
      there are no missing genes in study ", studyName, 
      " so imputation error rate not calculated. \n",
      file=outputFile,append=TRUE)
      errorRate <- NA
      meanAbsDiff <- NA
      
    }
    
    #end of gene filtering statements.
  }
  
  if(length(na.omit(match(colnames(exprFilterImpute),colnames(expr))))==0){
    
    stop("\nEither all of your samples were removed, or 
    somehow the sample names got changed after imputing data.\n")
    
  }
  #was a class index supplied? (i.e. is 
  #there a class list in this study object?)
  #then we need to also filter out the removed patients here.
  if(missing(classIndex)){
    
    cat("finished imputing study ",studyName, "\n",file=outputFile,append=TRUE)
    
    message("no list index name for a class/outcomes given,
    so this will not be returned")
    
    study <- list(expr=exprOrig,exprFilterImpute = exprFilterImpute,
    keysFilterImpute=keysFilterImpute,keys=keysOrig,meanAbsDiff=meanAbsDiff,
    errorRate=errorRate)
    return(study)
    
  }else{
    
    classesFilter <- list()
    
    
    
    for(c in 1:length(classIndex)){
      
      if( !(classIndex[c] %in% names(study)) || 
            (length(which(names(study) ==classIndex[c])) >1) ){
        print(classIndex[c])
        stop("you either supplied an incorrect class/outcomes index name or
             there are duplicated names in your study list.")
        
      }
      #get the correctly subsetted classes/outcomes variables 
      #for all designated variables.
      classes <- study[[classIndex[c]]]
      cat("getting the re-ordered and trimmed down class labels.",
          append=TRUE,file=outputFile)
      
      classesFilter[[c]] <- study[[classIndex[c]]][
        na.omit(match(colnames(exprFilterImpute),colnames(study[[exprIndex]]))),
        ,drop=FALSE]
      
      
    }
    
    
    study <- list(expr=exprOrig,exprFilterImpute = exprFilterImpute, 
                  class=classes,classesFilter=classesFilter,
                  keysFilterImpute=keysFilterImpute,keys=keysOrig,
                  meanAbsDiff=meanAbsDiff,errorRate=errorRate)
    
    cat("finished imputing study ",studyName, "\n",file=outputFile,append=TRUE)
    return(study)
    
  }
  
  
  #end of function
}


#remove NAs is remove NA keys.
#' Collapse/handle duplicated probes (genes) in a dataset
#'
#' Used internally by processExpressionSet.  Code to either take the average across
#' a set of duplicated "keys" (can be probes or genes, which correspond to the rows
#' in the expression matrix "expr"), or take the keys that has the highest variance 
#' across the set of duplicated keys.
#'
#' @param expr An expression matrix with genes in the rows and samples in the columns.
#' @param sampleColNames Sample column names. Needed for internal debugging; usually the default colnames(expr) is appropriate.
#' @param keys Generally the list of gene symbols, or some molecular key, that needs to be "collapsed" because it contains duplicated names.
#' @param method Method used to collapse probes: take the mean across all duplicated keys, or just pick the key with the highest variance?
#' @param debug Use internal unit tests that will stop the code if it detects a bug?
#' @param removeNA_keys Remove any NA keys?
#' @param varMetric Standard options taken from the base var() function. May be important if you have NA values in your data matrix; otherwise, "everything" is usually fine.
#'
#' @return Returns a processed list with the items "expr" and "keys", the expression matrix and final keys list.
#'
#' @author Katie Planey <katie.planey@@gmail.com>
#' @examples
#' #load up our datasets
#' curatedBreastDataExprSetList <- getCuratedBreastDataExprSetList(test=TRUE)
#' 
#' #just perform on second dataset, GSE2034, as an example.
#' #This dataset has no NAs already but does have duplicated genes
#' #highestVariance calculation make take a minute to run.
#' collapsedData <- collapseDupProbes(expr=exprs(curatedBreastDataExprSetList[[2]]),  
#' keys=fData(curatedBreastDataExprSetList[[2]])[, "gene_symbol"], 
#' method = c("highestVariance"), debug = TRUE, removeNA_keys = TRUE, 
#' varMetric = c("everything"))
#' #look at names of outputs
#' names(collapsedData)
#' @export
collapseDupProbes <- function(expr,sampleColNames=colnames(expr),keys, 
                              method=c("average","highestVariance"),debug=TRUE,
                              removeNA_keys=TRUE,
                              varMetric = c("everything", "all.obs", 
                                            "complete.obs", "na.or.complete", 
                                            "pairwise.complete.obs")){
  
  message("It's best to impute NA values before running this function
otherwise it may set averages to NA if there is 1 NA present.
This function just removes any genes whose key is NA.")
  
  #make sure it's an expression matrix to start with
  #(otherwise with 1 sample won't work).
  #and make sure it's numeric! the data.matrix solves that...
  #but have to make it a data.frame first. weird!
  #if don't as stringsAsFactors=FALSE...get totally wrong values.
  expr <- data.matrix(as.data.frame(expr,stringsAsFactors=FALSE))
  
  if(length(colnames(expr))==0 && missing(sampleColNames)){
    
    stop("error: must have a sample (column) ID specified for this sample.")
    
  }
  
  if(dim(expr)[1] != length(keys)) 
    stop("error: length of keys doesn't match 
         number of rows in expression matrix.")
  
  if(ncol(expr)==1){
    
    message("\nOnly 1 sample. May encounter edge cases 
            when collapsing duplicated probes.")

  }
  
  if(removeNA_keys){
    
    #catch both NA and "NA"
    #need to alternate between filtering expr and keys 
    #so that they keep the same indices
    #as.matrix() prevents it from turning into a vector after the subsetting. 
    #if you only have a 1-column matrix 
    expr <- as.matrix(expr[which(!is.na(keys)),]) 
    keys <- keys[which(!is.na(keys))]
    expr <- as.matrix(expr[which(keys != "NA"),])
    keys <- keys[which(keys != "NA")]
    #also remove "" keys just in case those 
    #somehow got coded up as "" instead of NA.
    expr <- as.matrix(expr[which(keys != ""),])
    keys <- keys[which(keys != "")]
    
  }
  
  singles.keys <- names(which(table(keys) == 1))
  singles.ind <- which(keys %in% singles.keys)
  
  message("\nYou may get a warning here because key (usually gene) names are
  duplicated so it can't use them as row names. 
  That's OK, because we are immediately collapsing them into one feature.\n")
  gems <- data.frame(expr = expr, keys = keys,stringsAsFactors = FALSE)
  #hmm this isn't working...but aren't all the lengths the same now???
  #add data.frame() again to force it into dataframe unless it's one column.
  out <- data.frame(gems[ singles.ind, ])
  singleKeys <- out$keys
  rownames(out) <- singleKeys
  #now remove the (last) keys column so only have expression values
  #to convince yourself of this, just type in names(out) 
  #and you'll see what I mean.
  out <- out[,1:(dim(out)[2]-1)]
  
  #re-set to patient names (data.frame action put "expr." as prefix)
  #only 1 sample? that means it'll pop up as a numeric class, 
  #not a matrix. re-set indices.
  
  if(is.numeric(out)){
    
    #uggh so much workaround to keep row names for single columns!
    out <- as.matrix(out)
    rownames(out) <- singleKeys
    
  }
  
  #make sure can add sample ID for later reference.
  if(!missing(sampleColNames)){
    
    colnames(expr) <- sampleColNames
    
  }
  
  if(ncol(expr)==0){
    
    stop("\nBug in your code that removes all samples.")
  }
  
  colnames(out) <-colnames(expr)
  
  ## Next, deal with multiple keys within a study ##
  #don't include probes for right now - makes it tricky.
  multis <- gems[ -singles.ind, ]
  
  if(nrow(multis) > 0){  # skip if no multiple ID found
    
    #split will give a list of the length of all the duplicated keys 
    #(1 index per key name - THIS part isn't duplicated now!
    #duplicated values will be in the element's list for that specific key name.)
    multis$keys <- as.factor(multis$keys)
    #levels may be holding on to keys you removed earlier, 
    #so drop these that have no data attached to them now.
    tmp <- split(multis, multis$keys,drop=TRUE)
    
    if (method == "average") {
      
      out2 <- sapply(tmp, function(m) {
        
        #odd bc if you run this code outside of sapply: 
        #it will have flipped indices for m.
        #inside sapply: keys are in the ROWS. 
        #so colMean takes mean INTRA patient across the duplicated keys 
        #- what we want!
        
        #REMOVE LAST COLUMN: this is the key string character, 
        #not an expression value!! 
        m <- m[, (1:dim(m)[2]-1) ]
        #just average across all duplicated probes (for each patient/column: 
        #want a new collapsed value for all duplicated keys, for each patient.)
        
        #if only one sample: set as a matrix. if there's only 1 col but > 1 row:
        #means 1 GENE, not 1 sample.
        #we want to remedy the 1 sample, duplicated genes scenario. 
        #there should be no just 1 gene (row) scenario anyways:
        #we're only looking at multis here...
        
        #well have to re-set every time as a matrix anyways to make all numeric 
        #again so numbers don't end up weird!!
        m <- data.matrix(as.data.frame(m,stringsAsFactors=FALSE))
        #this colMeans is just the duplicates across each SEPARATE patient 
        #(not lumping patient expression values together)
        out2 <- colMeans(m,na.rm=TRUE)
        
      })
      
      #sapply flips it - puts the keys as columns.
      
      #but if it's one PATIENT (sample): that means we need to make it a matrix,
      #and don't transpose it after
      #...but then lose the row names...so make sure to add them back!
      
      
      if(is.numeric(out2)){
        
        out2 <- as.matrix(out2)
        
        
      }else{
        #sapply flips it- puts the keys as columns.
        out2 <- t(out2)
        
      }
      
      colnames(out2) <-colnames(expr)
      rownames(out2) <- names(tmp)
      
      
      if(debug){
        #randomly find a key that was duplicated, check its means.
        if(length(unique(multis$keys)) != dim(out2)[1]) stop ("function is not 
        collapsing duplicated keys to the expected number of unique keys")
        
        #try just the first duplicated key - we don't know if 
        #there will be more duplicated ones.
        keyInd <- grep(rownames(out2)[1],multis$keys)
        #is mean of these rows (excluding last column - 
        #the key characters) as expected?
        exprRows <- as.matrix(multis[keyInd,1:(dim(multis)[2]-1)])
        
        
        #if(dim(as.matrix(exprRows))[2]==1 && dim(as.matrix(out2))[1]>1){
        
        #if only 1 sample: have to make it a 2D matrix first.
        #means <- colMeans(as.matrix(exprRows)) 
        
        # }else{
        
        means <- as.numeric(colMeans(exprRows,na.rm=TRUE))
        
        #  }
        
        #ignore samples who actually just have an NA value for this key.
        #once I got a bug in this becuase I hadn't filtered out a few weird
        #keys with "" as names and not NA.
        if(any ((means==out2[1,]) == FALSE,na.rm=TRUE)) stop("function is not 
returning the expected mean of rows of expression values for duplicated keys ")
        
      }
      
      #end of average calculations.
    }
    
    
    ######
    
    if (method == "highestVariance") {
      
      #we need lapply here bc tmp is a list - can't just use apply.
      #sapply was too buggy with this variance function option - 
      #never quite figured out why.
      out2 <- lapply(tmp, function(m) {
        
        #just average across all duplicated probes (for each patient.)
        #var() computes the variance across columns...what we want - 
        #1 gene's variance across all patients
        #REMOVE LAST COLUMN: this is the key string character, 
        #not an expression value!!
        m <- m[, (1:dim(m)[2]-1) ]
        #make all numeric again so numbers don't end up weird!!
        m <- data.matrix(as.data.frame(m,stringsAsFactors=FALSE))
        #here, we want to look at the variance (diag of cov) of each 
        #row (i.e. across all the columns (patients))
        #remember covariance is t(X)X - want numrows*numPatients = X
        #m will have multiple rows...so use apply.
        geneVar <- apply(m, MARGIN=1, FUN = function(geneRow,varMetric){
          
          
          varGene <- var(geneRow,use=varMetric)
          
          
          return(varGene)  
        }, varMetric)
        
        #if there's a tie: just randomly pick the first one - 
        #otherwise will return double the patients!       
        rowIndex <- which(geneVar == max(geneVar))[1]
        
        #need to make this a vector - otherwise hard to format into a data.frame
        #when writing up to SQL table later...
        #and sapply tries to make m a list - not good here!
        #cat("\n",dim(m))
        out2 <- m[rowIndex, ]
        
      })
      
      
      #lapply flips it - puts the keys as columns.
      #also need to convert lists to 1 matrix
      out2 <- do.call(rbind,out2)
      colnames(out2) <-colnames(expr)
      rownames(out2) <- names(tmp)
      
      if(debug){
        #randomly find a key that was duplicated, check its means.
        if(length(unique(multis$keys)) != dim(out2)[1]) stop ("function is not 
        collapsing duplicated keys to the expected number of unique keys")
        
        #try just the first duplicated key - we don't know if there will
        #be more duplicated ones.
        keyInd <- grep(rownames(out2)[1],multis$keys)
        #is the row chosen out of this subset as expected?
        exprRows <- multis[keyInd,1:(dim(multis)[2]-1)]
        
        var <- diag(cov(t(exprRows),use=varMetric))
        
        rowIndex <- which(var == max(var))[1]          
        bestRow <- exprRows[rowIndex, ]
        
        #so do the values for each patient match up?
        #ignore samples who actually just have an NA value for this key.
        if( any((bestRow==out2[1,]) == FALSE ,na.rm=TRUE)) stop("function 
      is not returning the expected highest row variance for duplicated keys ")
        
      }
      
      #end of variance option. 
    }
    
    #save some original values to test this function does what we want.
    exprOut <- rbind(out, out2)
    #backstory: the m[row, ] part in method = "highVariance" makes like nested 
    #lists...and RmySQL goes crazy when you do this
    #tne entire expr table can be a list, but it shouldn't be a list of lists..
    #..found this trick on stack overflow:
    #http://stackoverflow.com/questions/8305735/error-writing-to-csv
    #print(dim(exprOut))
    #exprOutFirst <- exprOut
    #careful: below does weird stuff...also remove row names!
    #exprOut <- as.data.frame(lapply(exprOut,unlist))
    #add row names afterwards?
    
    #save(exprOutFirst, exprOut,file="exprOut.RData")
    #keys are in row names of out,out2.
    keysOut <- append(rownames(out),rownames(out2))
    #hmm getting matching probes may be harder...probesOut <- 
    #end of multis processing 
    
    #no duplicates found.
  }else{
    
    exprOut <- out
    keysOut <- rownames(out)
    
  }

  if(debug==TRUE){

    
    if(length(keysOut) != dim(exprOut)[1]) stop("function is not returning 
    the same number of expression rows as keys")
    
    if(dim(expr)[2] != dim(exprOut)[2]) stop("function is not returing the 
    orginal number of patients/samples (columns in expr matrix)")
    
    if(length(unique(keysOut)) != length(keysOut)) stop("function is not 
    returning the correct number of unique keys in final step.")
    
    
  }
  
  
  
  gems <- list(expr = exprOut, keys = keysOut)
  
  return(gems)
  
}

#picks patient's sample that has the highest variance across all features.
#' Remove duplicated patient samples (samples from the same patient/column ID)
#'
#' Function to keep only 1 sample per patient (column ID) in the data matrix.
#' Keeps the sample that has the overall highest variance.
#'
#' @param exprMatrix Expression matrix, with features in the rows and samples in the columns.
#' @param outputFile Output file for messages that print status of removing duplicated samples. Include full directory if file should not be printed to current working directory.
#' @param varMetric Standard options taken from the base var() function. May be important if you have NA values in your data matrix; otherwise, "everything" is usually fine.
#'
#' @return exprMatrix:  the final data matrix with only 1 sample per patient ID.
#'
#' @author Katie Planey <katie.planey@@gmail.com>
#' @note Suggestions are welcome for further ways to pick the best sample from samples
#' from the same patient. No curatedBreastData matrices currently have samples
#' that share the same patient ID, but this function is especially useful for
#' say TCGA data, where this is often the case. 
#' 
#' It is suggseted one imputes missing values using the filterAndImpute function
#' before running this function to avoid -Inf and NA values in the variance
#' calculations.
#' @examples
#' #No curatedBreastData has duplicated samples, 
#' #but we can still run this function on one of the datasets:
#' #load up our datasets
#' curatedBreastDataExprSetList <- getCuratedBreastDataExprSetList(test=TRUE)
#' 
#' #This dataset does not have NA values, which makes for a good example without
#' #extra pre-processing.
#' outputMatrix <- removeDuplicatedPatients(exprMatrix=
#' exprs(curatedBreastDataExprSetList[[1]]), 
#' outputFile = tempfile(), varMetric = c("everything"))
#' #final dimensions - unchanged in this case with 
#' #no samples sharing the same patient ID.
#' dim(outputMatrix)
#' @export
removeDuplicatedPatients <- function(exprMatrix, 
                                     outputFile="duplicatedPatientsOutput.txt", 
                                    varMetric = c("everything", "all.obs",
                                                  "complete.obs", 
                                                  "na.or.complete", 
                                                  "pairwise.complete.obs")){
  
  message("\nStarting with  ", ncol(exprMatrix) ,"patients.")
  cat("\nStarting with ", ncol(exprMatrix) ,"patients. \n\n",
      append=TRUE,file=outputFile)
  
  
  if(length(varMetric)>1){
    
    message("defaulting the everything variance metric.")
    varMetric = c("everything")
    
  }
  nSamples <- ncol(exprMatrix)
  
  if(all(is.na(exprMatrix))){
    
    stop("your counts matrix  is empty.")
  }
  #all.obs this means that missing observations will produce an error. 
  #we want that! there should be NO NAs.
  
  #take only the UNIQUE set of these (if there's 3 samples from 1 patient, 
  #will show up twice in duplicated list)
  duplicated_patients <- unique(colnames(exprMatrix)[
    which(duplicated(colnames(exprMatrix)))])
  
  
  
  if(length(duplicated_patients)>0){
    
    message("Found duplicated samples for ",length(duplicated_patients),
            "patients. There may be more than 2 samples per patient.")
    cat("Found duplicated samples for patients ",duplicated_patients,
        "patients. There may be more than 2 samples per patient.\n",
        append=TRUE,file=outputFile)
    
    #remember: duplicated only returns the index of the n-1 duplicated items, 
    #NOT the first time the item appeared!
    #why need the "which" in the for loop to grap all samples 
    #related to 1 patient.
    
    for(d in 1:length(duplicated_patients)){
      
      
      #how many NAs in each?
      dupIndexes <- which(colnames(exprMatrix)== duplicated_patients[d])
      
      sampleVars <- apply(exprMatrix[,dupIndexes], MARGIN=2, 
                          FUN = function(sampleCol,varMetric){
        
        
        varSample <- var(sampleCol,use=varMetric)
        
        
        return(varSample)  
      }, varMetric)
      
      
      dup_sample_highestVar <- which(sampleVars==max(sampleVars))
      
      if(length(dup_sample_highestVar)>1){
        
        #a tie. ha probably means a bug...that would be rare! 
        #or you literally just duplicated a data column but tacked
        #on a different aliquot name...
        #just take arbitrary first sample then.
        dup_sample_highestVar <- dup_sample_highestVar[1]
        
      }
      
      
      #remove the other duplicated samples besides the ith one.
      dupSampleRemove <- dupIndexes[-dup_sample_highestVar]
      
      if(length(dupSampleRemove)==0){
        
        stop("bug in your remove duplicated patients code...no indices 
        left to remove the duplicated patient sample.")
        
      }
      
      exprMatrix <- exprMatrix[,-dupSampleRemove]
      
      
      
      #end of loop
    }
    
    
    #how many samples did we lose?
    numSamplesRemove <- nSamples - ncol(exprMatrix)
    message("removed ",numSamplesRemove, " samples from patients with 
            multiple samples. \n There are ", ncol(exprMatrix) ,"left.")
    cat("removed ",numSamplesRemove, " samples from patients with multiple 
        samples. \n There are ", ncol(exprMatrix) ,"left. \n\n",
        append=TRUE,file=outputFile)
    
    #end of if statement.
  }else{
    
    message("found no multiple samples from the same patient(s)")
    cat("found no multiple samples from the same patient(s) \n\n",
        append=TRUE,file=outputFile)
    
  }
  
  #return study (processed or not depending on if found duplicated)
  return(exprMatrix)
}


#note: this can handle duplicated keys: it never actually uses the key names, 
#so they don't need to be unique.
#this filters by variance, not coefficient of variance.
#if not using logged values, sometimes filtering by coeff of var is better.
#to-do: replace missing() with each variance option as NULL to make running code 
#in another function less cumbersome
#' Filter genes by variance
#'
#' A function that filters genes by variance; it can simply threshold out genes
#' that are above or below a certain magnitude of variance, filter out genes that
#' fall outside of a minimum and maximum percentile, or simply select the top
#' N varying genes.
#'
#' @param study A list, of minimally the gene expression or some molecular data matrix with keys (molecular features, such as genes) in the rows and patient samples in the columns and a keys list. In line with using this function for any type of molecular data, the exprIndex name, and also the keysIndex name, in the list can be altered.
#' @param plotSaveDir If plotVarianceHist is TRUE, then the plotSaveDir is a character string specifying where this histogram plot should be saved.
#' @param minVarPercentile Minimum variance percentile. Must be provided in conjunction with maxVarPercentile to use percentiles to threshold genes.
#' @param maxVarPercentile Maximum variance percentile. Defaul is 1, i.e. 1\%. Must be provided in conjunction with minVarPercentile to use percentiles to threshold genes.
#' @param maxVar If maxVar is provided, as opposed to minVarPercentile and maxVarPercentile, genes are removed that are above a certain variance magnitude. This may be useful if a user suspects very highly varying genes are actually technical noise/outliers. May be used in conjunction with minVar or in isolation.
#' @param minVar If maxVar is provided, as opposed to minVarPercentile and maxVarPercentile, genes are removed that are below a certain variance magnitude. This is helpful before running certain algorithms, such as the popular Combat batch normalization technique, that can throw errors if genes with extremely low variances are in the data matrix. May be used in conjunction with maxVar or in isolation.
#' @param exprIndex Character string. List slot name for the data matrix, presumably an expression matrix.
#' @param keysIndex Character string. List slot name for the feature names, presumably probes or gene names.
#' @param outputFile Output file for messages that print status of the filtering.  Include full directory if file should not be printed to current working directory.
#' @param plotVarianceHist Plot the histogram of variances overall? Good for exploratory analyses to understand the distribution of variance across all data points. Default is FALSE to avoid saving a ggplot image for every function run.
#' @param varMetric Standard options taken from the base var() function. May be important if you have NA values in your data matrix; otherwise, "everything" is usually fine.
#' @param sampleCol Are samples in the columns of the expression matrix? If not, this function will first transpose the matrix, as the function assumes samples are in the columns features are in the rows.
#' @param numTopVarGenes A numeric value indicating the number of genes (features) to select; the function will only take this number of genes that have the highest variance across all genes.
#'
#' @return A list:  output <- list(study=study,filteredStudy=filteredStudy,p=p);
#' \item{study}{Original study list object}
#' \item{filteredStudy}{filteredStudy object, i.e. the gene expression and keys only for the desired filtered keys/features.}
#'
#' @author Katie Planey <katie.planey@@gmail.com>
#' @note Filtering by variance is equivalent to filtering on the coefficient of variation
#' if data is logged.  Further work includes automatically allowing the user to use
#' the coefficient of variation as opposed to baseline variation for a threshold.
#' 
#' It is highly suggested you use filterAndImputeSamples() beforehand to remove any
#' NA values, to avoid -Inf or NA variance calculations.
#' @examples
#' #load up our datasets
#' curatedBreastDataExprSetList <- getCuratedBreastDataExprSetList(test=TRUE)
#' 
#' #just perform on one dataset as an example, GSE1379. 
#' #This dataset does not have NA values, which makes for a
#' #good example without extra preprocessing.
#' #highestVariance calculation make take a minute to run.
#' #create study list object. 
#' study <- list(expr=exprs(curatedBreastDataExprSetList[[1]]),
#' keys=fData(curatedBreastDataExprSetList[[1]])[, "gene_symbol"])
#' #take top 100 varying genes
#' 
#' filterGeneStudy <- filterGenesByVariance(study, exprIndex = "expr", 
#' keysIndex = "keys", outputFile = tempfile(), 
#' plotVarianceHist = FALSE,
#' varMetric = c("everything"), sampleCol = TRUE, numTopVarGenes=100)
#' 
#' #names of output
#' names(filterGeneStudy)
#' @export
filterGenesByVariance <- function(study, plotSaveDir="~/",minVarPercentile,
                                  maxVarPercentile=1,maxVar,minVar,
                                  exprIndex = "expr", keysIndex = "keys", 
                                  outputFile="varCal.txt",
                                  plotVarianceHist=FALSE,
                                  varMetric = c("everything", "all.obs", 
                                                "complete.obs", 
                                                "na.or.complete", 
                                                "pairwise.complete.obs"),
                                  sampleCol=TRUE,
                                  numTopVarGenes){
  
  
  if(sampleCol){
    #case were samples are the columns, not the rows.
    
    expr <- study[[exprIndex]];
    
  }else{
    
    expr <- t(study[[exprIndex]]);
    
    message("dimensions of expression data will be returned transposed: 
    samples are now the columns for a pxn matrix.");
    
  }
  
  
  keys <- study[[keysIndex]];
  
  if(all(is.na(expr)) || all(is.na(keys))){
    
    stop("you didn't provide the right expr and/or keys keyword for the 
         filter genes by variance function.")
    
  }
  
  cat("calculating variance across all genes for plotting","\n",
      file=outputFile,append=TRUE);
  
  
  #we want a variance calculation for each gene (row here)
  geneVar <- apply(expr, MARGIN=1, FUN = function(geneRow,varMetric){
    
    
    varGene <- var(geneRow,use=varMetric);
    
    
    return(varGene);  
  }, varMetric);
  
  names(geneVar) <- keys;
  cat("finished calculating variance","\n",file=outputFile,append=TRUE);
  
  #remove any genes with NA variance - cannot estimate these...
  
  geneVarNA <- which(is.na(geneVar));
  
  if(length(geneVarNA)>0){
    
    expr <- expr[-geneVarNA,];
    
    keys <- keys[-geneVarNA];
    
  }
  
  cat("found ",length(geneVarNA), "genes with NA variance. Removing them. \n",
      file=outputFile,append=TRUE);
  
  if(plotVarianceHist){
    
    cat("plotting a histogram of the variance for each gene 
    across all patients.","\n",file=outputFile,append=TRUE)
    p <- ggplot(as.data.frame(geneVar),aes(x=geneVar)) + 
      geom_histogram(binwidth=round((max(geneVar,na.rm=TRUE)-
                                       min(geneVar,na.rm=TRUE))/10)) +
      labs(title = "Histogram of gene expression variances",
           y ="number of genes", x="variance of expression");
    
    fileSave <- paste(plotSaveDir,"geneVarianceHist.png",sep="")
    ggsave(filename=fileSave);
    
  }else{
    
    p <- "no plot requested from function inputs";
    
  }
  
  
  if(missing(numTopVarGenes)){
    
    if(!missing(minVar)){
      #use an arbitrary variance cutoff
      
      medianGeneVar <- median(geneVar);
      
      
      geneVarMedianThresholded <- geneVar[which(geneVar > medianGeneVar)];
      cat("there were ",length(geneVar)," genes with no NA values 
     and the median variance of these genes was ", medianGeneVar, 
      "\n","there were", length(geneVarMedianThresholded) ,
      "genes above this median","\n" ,append=TRUE,file=outputFile);
      
      
      geneVarThreshholded <- geneVar[which(geneVar > minVar)];
      cat("there were ",length(geneVar),
          " genes with no NA values and there were ", 
          length(geneVarThreshholded) ,"genes above the mininum variance of ",
          minVar,"\n" ,append=TRUE,file=outputFile);
      
      expr <- expr[which(geneVar > minVar),];
      
      keys <- keys[which(geneVar > minVar)];
      
      filteredStudy <- list(expr=expr,keys=keys);
      
    }
    
    if(!missing(maxVar)){
      
      #only take genes UNDER a certain variance
      medianGeneVar <- median(geneVar);
      
      
      geneVarMedianThresholded <- geneVar[which(geneVar < medianGeneVar)];
      cat("there were ",length(geneVar)," genes with no NA values 
      and the median variance of these genes was ", medianGeneVar, "\n",
          "there were", length(geneVarMedianThresholded) ,
          "genes below this median","\n" ,append=TRUE,file=outputFile);
      
      
      geneVarThreshholded <- geneVar[which(geneVar > maxVar)];
      cat("there were ",length(geneVar),
          " genes with no NA values and there were ", 
          length(geneVarThreshholded) ,"genes below the maximum variance of ",
          minVar,"\n" ,append=TRUE,file=outputFile);
      
      expr <- expr[which(geneVar < maxVar),];
      
      keys <- keys[which(geneVar < maxVar)];
      
      filteredStudy <- list(expr=expr,keys=keys);
      
      
    }
    
    #percentiles will get messed up if don't do this together
    #default to min/max Var if that was also given.
    if(!missing(minVarPercentile) && !missing(maxVarPercentile)
       && missing(minVar) && missing(maxVar)){
      
      #do quantiles instead
      #if want .001: then times *1000+1 indices.
      varQuantiles <- quantile(geneVar,probs=seq(0,1,.01),na.rm=TRUE);
      #want to match up with varQuantiles index (.95 is 96th index)
      minVarPercentile <- minVarPercentile*100+1;
      minVar <- varQuantiles[minVarPercentile];
      
      
      filteredStudy <- list(expr=expr,keys=keys);
      #want to match up with varQuantiles index (.95 is 96th index)
      maxVarPercentile <- maxVarPercentile*100+1;
      maxVar <- varQuantiles[maxVarPercentile];
      
      #now filter out these genes. must do intersect: often, 
      #a lot of the genes will overlap!
      
      expr <- expr[intersect(which(geneVar > minVar),which(geneVar < maxVar)),]; 
      keys <- keys[intersect(which(geneVar > minVar),which(geneVar < maxVar))];
      filteredStudy <- list(expr=expr,keys=keys);
      
    }else if(missing(minVar) && missing(maxVar)){
      
      stop("you didn't specify the number of genes to threshold but didn't 
           #give a variance threshold metric either.")
    }
    
  }else{
    
    #just take top X varying genes.
    if(any(is.na(geneVar))){
      
      message("You have NA values in your gene variances. 
      variances with NA values will not be considered.")
      
    }
    #na.last=TRUE not allowed if returning an index like we are.
    geneVarSorted <- sort.int(geneVar,na.last=NA,decreasing=TRUE,
                              index.return=TRUE);
    #this way, can actually handle duplicated gene names. 
    #never sees the gene names, just the indices.
    if(length(geneVarSorted$ix)>= numTopVarGenes){
      
      topGeneIndices <- geneVarSorted$ix[1:numTopVarGenes];
      
    }else{
      #just take the rest of the genes
      cat("there weren't ",numTopVarGenes, " to choose from, 
      so just took all of the remaining genes.\n",append=TRUE,file=outputFile)
      topGeneIndices <- geneVarSorted$ix;
      
    }
    
    expr <- expr[topGeneIndices, ];
    keys <- keys[topGeneIndices];
    
    #come back and TEST:
    #take first index - there may be a tie
    topGeneIndex <- which(geneVar ==max(geneVar,na.rm=TRUE))[1];
    
    if(! (topGeneIndex %in% topGeneIndices) ){
      
      stop("looks like there's a bug in how you're 
    sorting the top gene indices.")
    }
    
    filteredStudy <- list(expr=expr,keys=keys);
    
    output <- list(study=study,filteredStudy=filteredStudy,p=p);
    
    return(output);
    
  }
  
  output <- list(study=study,filteredStudy=filteredStudy,p=p);
  
  return(output);
  
}

