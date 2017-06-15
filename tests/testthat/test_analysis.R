#==============================================================================
context("Reading in data")

#==============================================================================
test_that("Ratios are calculated correctly", {
  
  file1 = "../testdata/experiment/condition1/rep2/spikein-col1-fwd.csv"
  file2 = "../testdata/experiment/condition1/rep2/spikein-col1-rev.csv"
  result1  = load_counts_fwdandrev(file1, file2)
  
  file1 = "../testdata/experiment/condition1/rep2/spikein-col2-fwd.csv"
  file2 = "../testdata/experiment/condition1/rep2/spikein-col2-rev.csv"
  result2  = load_counts_fwdandrev(file1, file2)
  
  s1 = subset(result1,select=c(Geneid,sense))
  s2 = subset(result2,select=c(Geneid,sense))
  a1 = subset(result1,select=c(Geneid,anti))
  a2 = subset(result2,select=c(Geneid,anti))
  
  s = merge(s1,s2,by="Geneid")
  a = merge(a1,a2,by="Geneid")
  g = c("Col","Col")
  
  fullresult = analysenew(data.matrix(s[,2:3]),data.matrix(a[,2:3]),data.matrix(s[,1]), g)
  result = fullresult[[1]]
  # overall ratio calc
  expect_equal(nrow(result), 1)
  expect_equal(result$ratio[[1]], 0.0003733957, tolerance=1e-6)
  # check sum of individual ratio calcs - note these can contain NaN
  expect_equal(sum(fullresult[[2]][[1]]$ratio,na.rm=TRUE), 0.02704632, tolerance=1e-6)
  
  g = c(1,2)
  fullresult = analysenew(data.matrix(s[,2:3]),data.matrix(a[,2:3]),data.matrix(s[,1]), g)
  result = fullresult[[1]]
  # overall ratio calc
  expect_equal(nrow(result),2)
  expect_equal(result$ratio[[1]], 0.0007461418, tolerance=1e-6)
  expect_equal(result$ratio[[2]], 0.0007817986, tolerance=1e-6)
  # check sum of individual ratio calcs - note these can contain NaN
  expect_equal(sum(fullresult[[2]][[1]]$ratio,na.rm=TRUE), 0.05486877, tolerance=1e-6)
  expect_equal(sum(fullresult[[2]][[2]]$ratio,na.rm=TRUE), 0.02777858, tolerance=1e-6)
})