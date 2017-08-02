# This file is part of RoSA.
# 
# RoSA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# RoSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with RoSA  If not, see <http://www.gnu.org/licenses/>.

#==============================================================================
context("Reading in data")

##############################################################################################################
# load_counts_fwdandrev: load spike-in data from fwd and rev files
# file1: forward counts file, stripped of initial lines, tab delimited
# file2: reverse counts file, ditto
load_counts_fwdandrev <- function(file1, file2)
{
  f <- read.csv(file1, sep="\t")
  r <- read.csv(file2, sep="\t")
  
  names(f)[7]<- "sense"
  names(r)[8]<- "anti"
  
  l <- merge(f,r,by="Geneid")
  
  return(l)
}

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
  
  #TODO
})

#==============================================================================
test_that("Calculate ratios works correctly", {
  
  ids = data.frame(c("id1"))
  lengths = data.frame(c(100))
  sensecounts = array(c(10,20,30,40,50,60),dim=c(2,3))
  anticounts = array(c(2,5,4,1,10,8), dim=c(2,3))
  groups = c("A", "B", "C")
  totalcounts = array(c(50,60,70), dim=c(1,3))
  
  result = calculate_ratios(sensecounts, anticounts, ids, lengths, groups, totalcounts)
  
  expect_equal(result[[1]][1,1],"A")
  expect_equal(result[[1]][[1,2]],0.24)
  expect_equal(result[[1]][3,1],"C")
  expect_equal(result[[1]][[3,2]],0.1606557, tolerance=1e-5)
})