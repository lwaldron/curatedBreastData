test_that("Data retrieval functions return correct types and sizes", {
  esets <- getCuratedBreastDataExprSetList(test=TRUE)
  expect_is(esets, "list")
  expect_is(esets[[1]], "ExpressionSet")
  expect_equal(length(esets), 2)
  
  clin <- getClinicalData(test=TRUE)
  expect_is(clin, "list")
  expect_true("clinicalTable" %in% names(clin))
  expect_true("clinicalVarDef" %in% names(clin))
  expect_is(clin$clinicalTable, "data.frame")
  expect_is(clin$clinicalVarDef, "data.frame")
  
  data(curatedBreastData_pData, package = "curatedBreastData")
  expect_is(curatedBreastData_pData, "data.frame")
})
